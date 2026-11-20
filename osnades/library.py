"""Compound library: species with their σ-surfaces (.orcacosmo files, one per conformer).

Layout::

    <home>/library/index.json
    <home>/library/<slug>/conf0.orcacosmo, conf1.orcacosmo, ...
"""
from __future__ import annotations

import json
import re
import shutil
import threading
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Callable

import numpy as np

from . import orca
from .config import Settings
from .structures import ATOMIC_MASS, Geometry, parse_xyz, smiles_info, smiles_to_conformers

Log = Callable[[str], None]
_LOCK = threading.RLock()


def slugify(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", name.strip()).strip("_")
    return s or f"species_{int(time.time())}"


@dataclass
class Species:
    name: str
    conformers: list[str]  # paths relative to library dir
    smiles: str | None = None
    charge: int = 0
    multiplicity: int = 1
    formula: str = ""
    molar_mass: float = 0.0
    source: str = "orca"  # orca | imported
    method: str = ""
    optimizer: str = ""
    created: str = ""
    notes: str = ""
    # Optional experimental solid–liquid data (for solutes)
    dh_fus: float | None = None  # J/mol
    t_fus: float | None = None  # K
    dcp_fus: float | None = None  # J/(mol K)
    fusion_source: str = ""  # where dh_fus / t_fus come from (lookup reference, "user", preset note)

    def conformer_paths(self, library_dir: Path) -> list[Path]:
        return [library_dir / c for c in self.conformers]


class Library:
    def __init__(self, settings: Settings):
        self.settings = settings
        self.dir = settings.library_dir
        self.dir.mkdir(parents=True, exist_ok=True)
        self.index_path = self.dir / "index.json"

    # ------------------------------------------------------------------ storage
    def _load(self) -> dict[str, dict]:
        if not self.index_path.exists():
            return {}
        return json.loads(self.index_path.read_text())

    def _save(self, data: dict[str, dict]) -> None:
        tmp = self.index_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=2, sort_keys=True))
        tmp.replace(self.index_path)

    def all(self) -> list[Species]:
        with _LOCK:
            return [Species(**v) for _, v in sorted(self._load().items(), key=lambda kv: kv[0].lower())]

    def names(self) -> list[str]:
        return [s.name for s in self.all()]

    def get(self, name: str) -> Species:
        with _LOCK:
            data = self._load()
        if name not in data:
            raise KeyError(f"Species {name!r} not in library. Available: {', '.join(sorted(data)) or '(none)'}")
        return Species(**data[name])

    def has(self, name: str) -> bool:
        return name in self._load()

    @staticmethod
    def _folder_key(name: str) -> str:
        # macOS/Windows file systems are case-insensitive: "Glucose/" and "glucose/" are the same folder
        return slugify(name).lower()

    def clash(self, name: str) -> str | None:
        """Existing *other* species whose storage folder would collide with ``name`` (case-insensitive)."""
        key = self._folder_key(name)
        for other, v in self._load().items():
            if other == name:
                continue
            folders = {Path(c).parts[0].lower() for c in v.get("conformers", []) if Path(c).parts}
            if self._folder_key(other) == key or key in folders:
                return other
        return None

    def _check_new_name(self, name: str) -> None:
        other = self.clash(name)
        if other:
            raise ValueError(f"{name!r} would share the storage folder of existing species {other!r} "
                             f"(names differ only in case/spaces/punctuation). Use {other!r} with overwrite, "
                             "or choose a distinct name.")

    def put(self, sp: Species) -> None:
        with _LOCK:
            data = self._load()
            data[sp.name] = asdict(sp)
            self._save(data)

    def update(self, name: str, **fields) -> Species:
        sp = self.get(name)
        for k, v in fields.items():
            setattr(sp, k, v)
        self.put(sp)
        return sp

    def remove(self, name: str) -> None:
        with _LOCK:
            data = self._load()
            sp = data.pop(name, None)
            self._save(data)
        if sp:
            if self.is_solvent(name):
                self.remove_solvent(name)
            folder = self.dir / slugify(name)
            # never delete files another species still uses
            used = {Path(c).parts[0].lower() for v in data.values() for c in v.get("conformers", [])}
            if folder.is_dir() and folder.name.lower() not in used:
                shutil.rmtree(folder)

    # ------------------------------------------------------------------ user solvents
    @property
    def solvents_path(self) -> Path:
        return self.dir / "solvents.json"

    def user_solvents(self) -> dict[str, dict[str, float]]:
        """Solvent systems the user defined (name -> {species: molar ratio})."""
        if not self.solvents_path.exists():
            return {}
        return json.loads(self.solvents_path.read_text())

    def _save_solvents(self, data: dict) -> None:
        tmp = self.solvents_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=2))  # keep component order (HBA:HBD as entered)
        tmp.replace(self.solvents_path)

    def add_solvent(self, name: str, components: dict[str, float] | None = None, overwrite: bool = False) -> None:
        """Register a solvent (a pure species when ``components`` is None, or a molar-ratio mixture)."""
        from .thermo import presets

        name = name.strip()
        comps = {k: float(v) for k, v in (components or {name: 1.0}).items() if float(v) > 0}
        if not name or not comps:
            raise ValueError("A solvent needs a name and at least one component")
        with _LOCK:
            data = self.user_solvents()
            taken = {k.lower(): k for k in list(presets()["solvent_systems"]) + list(data)}
            other = taken.get(name.lower())
            if other in presets()["solvent_systems"]:
                raise ValueError(f"{name!r} is already a preset solvent ({other!r})")
            if other and not overwrite:
                raise ValueError(f"Solvent {other!r} already exists")
            if other and other != name:
                data.pop(other)
            idx = self._load()
            q = sum(idx[s]["charge"] * r for s, r in comps.items() if s in idx) / sum(comps.values())
            if abs(q) > 1e-9:
                raise ValueError(f"{name!r} is not electroneutral (net charge {q:+.3f} per mole)")
            data[name] = comps
            self._save_solvents(data)

    def remove_solvent(self, name: str) -> None:
        with _LOCK:
            data = self.user_solvents()
            if data.pop(name, None) is not None:
                self._save_solvents(data)

    def is_solvent(self, species: str) -> bool:
        return self.user_solvents().get(species) == {species: 1.0}

    def paths(self, name: str) -> list[Path]:
        return self.get(name).conformer_paths(self.dir)

    # ------------------------------------------------------------------ import
    def import_orcacosmo(
        self, name: str, files: list[Path], smiles: str | None = None, overwrite: bool = False, **meta
    ) -> Species:
        from opencosmorspy.input_parsers import SigmaProfileParser

        if self.has(name) and not overwrite:
            raise ValueError(f"{name!r} already exists")
        self._check_new_name(name)
        folder = self.dir / slugify(name)
        folder.mkdir(parents=True, exist_ok=True)
        rel = []
        spp = None
        for i, f in enumerate(files):
            dest = folder / f"conf{i}.orcacosmo"
            shutil.copyfile(f, dest)
            spp = spp or SigmaProfileParser(str(dest))
            rel.append(str(dest.relative_to(self.dir)))
        charge = int(round(float(np.sum(spp["seg_charge"])))) if spp is not None else 0
        symbols = list(spp["atm_elmnt"]) if spp is not None else []
        geom = Geometry(symbols, np.asarray(spp["atm_pos"])) if symbols else None
        info = {}
        if smiles:
            try:
                info = smiles_info(smiles)
            except Exception:  # noqa: BLE001
                info = {}
        sp = Species(
            name=name, conformers=rel, smiles=info.get("smiles", smiles), charge=charge,
            formula=info.get("formula") or (geom.formula if geom else ""),
            molar_mass=info.get("molar_mass") or (geom.molar_mass if geom else 0.0),
            source=meta.pop("source", "imported"), method=spp["method"] if spp is not None else "",
            created=time.strftime("%Y-%m-%d %H:%M"), **meta,
        )
        self.put(sp)
        return sp

    # ------------------------------------------------------------------ compute
    def compute(
        self,
        name: str,
        smiles: str | None = None,
        xyz: str | None = None,
        charge: int | None = None,
        multiplicity: int | None = None,
        n_conformers: int = 1,
        optimizer: str = "xtb",  # xtb | mmff | none
        overwrite: bool = False,
        lookup_fusion: bool = True,
        log: Log = print,
        **meta,
    ) -> Species:
        """Generate geometry → (optional) xTB optimisation → ORCA COSMO-RS → add to library."""
        if self.has(name) and not overwrite:
            raise ValueError(f"{name!r} already exists in the library")
        self._check_new_name(name)
        if not smiles and not xyz:
            raise ValueError("Provide a SMILES string or an XYZ geometry")

        info: dict = {}
        if smiles:
            info = smiles_info(smiles)
            if info["n_fragments"] > 1:
                raise ValueError(
                    "SMILES has several fragments (salt?). Add each ion as its own species "
                    "(COSMO-RS treats ions separately) or use `compute_salt`."
                )
            geoms = smiles_to_conformers(smiles, n_conformers=n_conformers)
            if optimizer == "none":
                optimizer = "mmff"
        else:
            geoms = [parse_xyz(xyz)]
            if optimizer == "mmff":
                optimizer = "none"
        charge = info.get("charge", 0) if charge is None else charge
        multiplicity = info.get("multiplicity", 1) if multiplicity is None else multiplicity

        stamp = time.strftime("%Y%m%d-%H%M%S")
        workroot = self.settings.jobs_dir / f"{slugify(name)}_{stamp}"
        files: list[Path] = []
        log(f"{name}: {len(geoms)} conformer(s), charge {charge}, optimiser {optimizer}")
        for i, g in enumerate(geoms):
            wd = workroot / f"conf{i}"
            if optimizer == "xtb":
                g = orca.xtb_optimize(self.settings, g, wd / "opt", "opt", charge, multiplicity, log=log)
            files.append(orca.cosmors_orcacosmo(self.settings, g, wd / "cosmo", "cosmo", charge, multiplicity, log=log))

        if overwrite and self.has(name):
            self.remove(name)
        sp = self.import_orcacosmo(
            name, files, smiles=smiles, overwrite=True, source="orca",
            optimizer=optimizer, multiplicity=multiplicity, **meta,
        )
        if sp.charge != charge:
            sp = self.update(name, charge=charge)
        log(f"{name}: added to library ({len(files)} σ-surface(s))")
        if lookup_fusion and sp.charge == 0 and not (sp.dh_fus and sp.t_fus):
            sp = self.lookup_fusion(name, log=log)
        return sp

    def lookup_fusion(self, name: str, overwrite: bool = False, log: Log = print, timeout: float = 15.0):
        """Fill in ΔH_fus / T_m of a species from published data (see :mod:`osnades.fusion`).
        Only empty fields are filled unless ``overwrite``. Network problems are logged, not raised."""
        from . import fusion

        sp = self.get(name)
        try:
            rec = fusion.find(sp.smiles, name, timeout=timeout)
        except Exception as exc:  # noqa: BLE001
            log(f"{name}: fusion-data lookup failed ({exc}); enter ΔH_fus and T_m by hand")
            return sp
        if rec is None:
            log(f"{name}: no published fusion data found; enter ΔH_fus and T_m by hand")
            return sp
        fields = {}
        if rec.dh_fus and (overwrite or not sp.dh_fus):
            fields["dh_fus"] = rec.dh_fus
        if rec.t_fus and (overwrite or not sp.t_fus):
            fields["t_fus"] = rec.t_fus
        if fields:
            partial = len(fields) == 1 and sp.fusion_source
            what = "ΔH_fus" if "dh_fus" in fields else "T_m"
            fields["fusion_source"] = f"{sp.fusion_source}; {what}: {rec.note}" if partial else rec.note
            sp = self.update(name, **fields)
            log(f"{name}: fusion data {rec.label}")
        return sp

    def ensure(self, name: str, smiles: str, log: Log = print, **kw) -> Species:
        """Compute a species if it is not in the library yet."""
        if self.has(name):
            return self.get(name)
        return self.compute(name, smiles=smiles, log=log, **kw)

    def xyz(self, name: str) -> str:
        """XYZ block of the first conformer (from the .orcacosmo file)."""
        from opencosmorspy.input_parsers import SigmaProfileParser

        spp = SigmaProfileParser(str(self.paths(name)[0]))
        g = Geometry(list(spp["atm_elmnt"]), np.asarray(spp["atm_pos"]), name)
        return g.to_xyz()
