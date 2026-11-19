"""Thermodynamics with openCOSMO-RS: activity coefficients, solid–liquid solubility,
solvent / NADES screening.

Solid solute dissolving in a (multi-component) solvent, ideal solid phase::

    ln(x_s γ_s) = -ΔH_fus/R (1/T - 1/T_m) + ΔCp_fus/R (T_m/T - 1 - ln(T_m/T))

γ_s is referenced to the pure (sub-cooled) liquid solute.  ``iterative`` solves the
equation self-consistently (γ_s evaluated at the saturated composition, bracketed
root search); ``infinite_dilution`` uses γ_s^∞ (the tutorial's non-iterative variant).

Internally every job is evaluated with the "cosmo" reference state and converted:
ln γ_s(pure ref) = ln γ_s^cosmo(mix) − ln γ_s^cosmo(pure s).  This avoids building
unphysical "pure ion" reference liquids when the solvent contains ions (ChCl NADES).
"""
from __future__ import annotations

import json
import threading
from dataclasses import dataclass, field, asdict
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from .library import Library

R = 8.314462618  # J/(mol K)
PRESETS_PATH = Path(__file__).resolve().parent / "data" / "presets.json"


@lru_cache(maxsize=1)
def presets() -> dict:
    return json.loads(PRESETS_PATH.read_text())


# ---------------------------------------------------------------------------- systems
@dataclass
class SolventSystem:
    name: str
    components: dict[str, float]  # species name -> molar ratio

    @property
    def species(self) -> list[str]:
        return list(self.components)

    @property
    def fractions(self) -> np.ndarray:
        r = np.array([float(v) for v in self.components.values()])
        if np.any(r < 0) or r.sum() <= 0:
            raise ValueError(f"Invalid ratios in {self.name}")
        return r / r.sum()

    def with_water(self, water_wt_percent: float, library: Library) -> "SolventSystem":
        """Return a copy with water added at the given mass percent of the solvent."""
        if water_wt_percent <= 0:
            return self
        if water_wt_percent >= 100:
            return SolventSystem("Water", {"water": 1.0})
        m_mix = sum(f * library.get(s).molar_mass for s, f in zip(self.species, self.fractions))
        m_w = library.get("water").molar_mass or 18.015
        w = water_wt_percent / 100.0
        # moles of water per mole of (dry) mixture
        n_w = (w / (1 - w)) * m_mix / m_w
        comps = {s: float(f) for s, f in zip(self.species, self.fractions)}
        comps["water"] = comps.get("water", 0.0) + n_w
        return SolventSystem(f"{self.name} + {water_wt_percent:g} wt% water", comps)

    def charge_balance(self, library: Library) -> float:
        return float(sum(library.get(s).charge * f for s, f in zip(self.species, self.fractions)))

    @staticmethod
    def parse(text: str) -> "SolventSystem":
        """Parse ``"choline:chloride:urea=1:1:2"``, ``"water"`` or a preset name."""
        text = text.strip()
        systems = presets()["solvent_systems"]
        if text in systems:
            return SolventSystem(text, dict(systems[text]))
        for k in systems:  # case-insensitive preset lookup
            if k.lower() == text.lower():
                return SolventSystem(k, dict(systems[k]))
        if "=" in text:
            names, ratios = text.split("=", 1)
            names = [n.strip() for n in names.split(":")]
            ratios = [float(r) for r in ratios.split(":")]
            if len(names) != len(ratios):
                raise ValueError(f"Number of species and ratios differ in {text!r}")
            return SolventSystem(text, dict(zip(names, ratios)))
        return SolventSystem(text, {text: 1.0})


def preset_systems() -> list[SolventSystem]:
    return [SolventSystem(k, dict(v)) for k, v in presets()["solvent_systems"].items()]


# ---------------------------------------------------------------------------- results
@dataclass
class SolubilityResult:
    solute: str
    solvent: str
    T: float
    x: float  # mole fraction of solute (solvent species counted individually)
    ln_gamma: float  # at the solution composition
    ln_gamma_inf: float
    x_ideal: float
    w: float  # mass fraction of solute
    g_per_kg: float  # g solute per kg solvent
    method: str
    note: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def ideal_rhs(T: float, dh_fus: float, t_fus: float, dcp_fus: float = 0.0) -> float:
    rhs = -dh_fus / R * (1.0 / T - 1.0 / t_fus)
    if dcp_fus:
        rhs += dcp_fus / R * (t_fus / T - 1.0 - np.log(t_fus / T))
    return float(rhs)


# ---------------------------------------------------------------------------- engine
class Engine:
    """openCOSMO-RS engine with a cache of loaded σ-surfaces per species set."""

    def __init__(self, library: Library, parameterization: str = "openCOSMORS24a"):
        self.library = library
        self.parameterization = parameterization
        self._cache: dict[tuple, object] = {}
        self._lock = threading.RLock()

    def _crs(self, species: tuple[str, ...]):
        from opencosmorspy import parameterization as P
        from opencosmorspy.cosmors import COSMORS

        key = tuple(
            (s, tuple(str(p) for p in self.library.paths(s))) for s in species
        ) + (self.parameterization,)
        if key not in self._cache:
            crs = COSMORS(par=getattr(P, self.parameterization)())
            for s in species:
                crs.add_molecule([str(p) for p in self.library.paths(s)])
            self._cache[key] = crs
            if len(self._cache) > 64:
                self._cache.pop(next(iter(self._cache)))
        return self._cache[key]

    # ------------------------------------------------------------- activity coeffs
    def ln_gamma(self, species: list[str], X: np.ndarray, T: float | np.ndarray,
                 refst: str = "cosmo", part: str = "tot") -> np.ndarray:
        """ln γ for each composition row of X (shape (n_jobs, n_species)).

        refst: "cosmo" (ideal conductor) or "pure_component" (pure liquid of each species).
        """
        X = np.atleast_2d(np.asarray(X, dtype=float))
        Ts = np.broadcast_to(np.asarray(T, dtype=float), (X.shape[0],))
        with self._lock:
            crs = self._crs(tuple(species))
            crs.clear_jobs()
            for x, t in zip(X, Ts):
                crs.add_job(x / x.sum(), float(t), refst=refst)
            res = crs.calculate()
            crs.clear_jobs()
        return np.asarray(res[part]["lng"])

    def ln_gamma_pure_ref(self, species: list[str], X: np.ndarray, T: float, idx: int = 0) -> np.ndarray:
        """ln γ of species ``idx`` referenced to its own pure liquid, for each row of X."""
        X = np.atleast_2d(np.asarray(X, dtype=float))
        pure = np.zeros(len(species))
        pure[idx] = 1.0
        lg = self.ln_gamma(species, np.vstack([X, pure]), T, refst="cosmo")
        return lg[:-1, idx] - lg[-1, idx]

    # ------------------------------------------------------------------ solubility
    def _composition(self, system: SolventSystem, xs: np.ndarray) -> np.ndarray:
        xs = np.atleast_1d(xs)[:, None]
        return np.hstack([xs, (1 - xs) * system.fractions[None, :]])

    def solubility(
        self, solute: str, system: SolventSystem, T: float,
        dh_fus: float | None = None, t_fus: float | None = None, dcp_fus: float | None = None,
        method: str = "iterative",
    ) -> SolubilityResult:
        sp = self.library.get(solute)
        dh_fus = dh_fus if dh_fus is not None else sp.dh_fus
        t_fus = t_fus if t_fus is not None else sp.t_fus
        dcp_fus = dcp_fus if dcp_fus is not None else (sp.dcp_fus or 0.0)
        if dh_fus is None or t_fus is None:
            raise ValueError(f"ΔH_fus and T_m are required for {solute} (set them in the library or pass them)")
        if solute in system.components:
            raise ValueError("Solute cannot also be a solvent component")
        species = [solute] + system.species
        rhs = ideal_rhs(T, dh_fus, t_fus, dcp_fus)
        x_ideal = float(min(1.0, np.exp(rhs)))

        lng_inf = float(self.ln_gamma_pure_ref(species, self._composition(system, np.array([0.0])), T)[0])
        note = ""
        if rhs >= 0:
            x, lng = 1.0, 0.0
            note = "T ≥ T_m: solute is liquid; SLE does not apply"
        elif method == "infinite_dilution":
            x, lng = float(min(1.0, np.exp(rhs - lng_inf))), lng_inf
        else:
            def g(lnx: float) -> float:
                lg = self.ln_gamma_pure_ref(species, self._composition(system, np.array([np.exp(lnx)])), T)[0]
                return lg + lnx - rhs

            grid = np.linspace(-46.0, 0.0, 47)
            lgs = self.ln_gamma_pure_ref(species, self._composition(system, np.exp(grid)), T)
            vals = lgs + grid - rhs
            idx = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
            if len(idx) == 0:
                x = float(np.exp(grid[0])) if vals[0] > 0 else 1.0
                lng = float(lgs[0])
                note = "no root found in [1e-20, 1]"
            else:
                i = idx[0]
                if len(idx) > 1:
                    note = "several roots (possible liquid–liquid demixing); lowest-x root reported"
                lnx = brentq(g, grid[i], grid[i + 1], xtol=1e-10)
                x = float(np.exp(lnx))
                lng = float(rhs - lnx)

        # mass-based units
        m_s = sp.molar_mass
        m_mix = float(sum(f * self.library.get(s).molar_mass for s, f in zip(system.species, system.fractions)))
        if m_s and m_mix and x < 1:
            w = x * m_s / (x * m_s + (1 - x) * m_mix)
            gkg = 1000.0 * x * m_s / ((1 - x) * m_mix)
        else:
            w, gkg = (1.0 if x >= 1 else float("nan")), float("inf" if x >= 1 else "nan")
        return SolubilityResult(solute, system.name, T, x, lng, lng_inf, x_ideal, w, gkg, method, note)

    def screen(self, solute: str, systems: list[SolventSystem], T: float, method: str = "iterative",
               **fus) -> pd.DataFrame:
        rows = []
        for sys_ in systems:
            try:
                r = self.solubility(solute, sys_, T, method=method, **fus)
                rows.append(r.as_dict())
            except Exception as exc:  # noqa: BLE001
                rows.append({"solute": solute, "solvent": sys_.name, "T": T, "note": f"error: {exc}"})
        df = pd.DataFrame(rows)
        if "x" in df:
            df = df.sort_values("x", ascending=False, na_position="last").reset_index(drop=True)
        return df

    def solubility_curve(self, solute: str, system: SolventSystem, temps: np.ndarray, method: str = "iterative",
                         **fus) -> pd.DataFrame:
        return pd.DataFrame([self.solubility(solute, system, float(T), method=method, **fus).as_dict() for T in temps])

    # ------------------------------------------------------------------ binaries
    def binary(self, a: str, b: str, T: float, n: int = 41) -> pd.DataFrame:
        """ln γ (pure-liquid reference) and G^E/RT across a binary at temperature T."""
        x1 = np.linspace(0.0, 1.0, n)
        X = np.column_stack([x1, 1 - x1])
        pure = np.array([[1.0, 0.0], [0.0, 1.0]])
        lg = self.ln_gamma([a, b], np.vstack([X, pure]), T, refst="cosmo")
        lg1 = lg[:n, 0] - lg[n, 0]
        lg2 = lg[:n, 1] - lg[n + 1, 1]
        ge = x1 * lg1 + (1 - x1) * lg2
        return pd.DataFrame({"x1": x1, f"ln_gamma_{a}": lg1, f"ln_gamma_{b}": lg2, "GE_RT": ge})

    def mixture_activity(self, system: SolventSystem, T: float) -> pd.DataFrame:
        """ln γ of every species in the mixture (cosmo reference) plus residual/combinatorial split."""
        X = system.fractions[None, :]
        with self._lock:
            crs = self._crs(tuple(system.species))
            crs.clear_jobs()
            crs.add_job(X[0], T, refst="cosmo")
            res = crs.calculate()
            crs.clear_jobs()
        out = {"species": system.species, "x": system.fractions}
        for part in ("tot", "res", "comb"):
            if part in res:
                out[f"ln_gamma_{part}"] = np.asarray(res[part]["lng"])[0]
        return pd.DataFrame(out)
