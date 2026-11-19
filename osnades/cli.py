"""Command line interface: ``osnades <command>``."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .config import load_settings
from .library import Library
from .thermo import Engine, SolventSystem, presets, preset_systems


def _lib():
    s = load_settings()
    return s, Library(s)


def cmd_setup(a):
    from . import orca

    s = load_settings()
    if a.orca_home:
        s.orca_home = Path(a.orca_home).expanduser()
    if a.mode:
        s.orca_mode = a.mode
    if a.ncores:
        s.ncores = a.ncores
    print(f"Saved settings to {s.save()}")
    if s.resolved_mode == "docker" and not a.skip_build:
        orca.build_docker_image(s)
    st = orca.check_backend(s, quick_run=not a.no_test)
    print("\n".join(st.messages))
    print("Backend OK" if st.ok else "Backend NOT ready")
    return 0 if st.ok else 1


def cmd_status(a):
    from . import orca

    s = load_settings()
    for k, v in s.to_dict().items():
        print(f"{k:16s} {v}")
    print(f"{'resolved_mode':16s} {s.resolved_mode}")
    st = orca.check_backend(s, quick_run=a.test)
    print("\n".join(st.messages))
    return 0 if st.ok else 1


def cmd_add(a):
    s, lib = _lib()
    pre = presets()["species"].get(a.name, {})
    smiles = a.smiles or (None if a.xyz else pre.get("smiles"))
    xyz = Path(a.xyz).read_text() if a.xyz else None
    meta = {k: pre[k] for k in ("dh_fus", "t_fus", "notes") if k in pre}
    if a.dh_fus is not None:
        meta["dh_fus"] = a.dh_fus
    if a.t_fus is not None:
        meta["t_fus"] = a.t_fus
    lib.compute(a.name, smiles=smiles, xyz=xyz, charge=a.charge, multiplicity=a.mult,
                n_conformers=a.conformers, optimizer=a.optimizer, overwrite=a.overwrite, **meta)
    return 0


def cmd_import(a):
    s, lib = _lib()
    meta = {}
    if a.dh_fus is not None:
        meta["dh_fus"] = a.dh_fus
    if a.t_fus is not None:
        meta["t_fus"] = a.t_fus
    sp = lib.import_orcacosmo(a.name, [Path(f) for f in a.files], smiles=a.smiles, overwrite=a.overwrite, **meta)
    print(f"Imported {sp.name}: {len(sp.conformers)} conformer(s), charge {sp.charge}, M={sp.molar_mass:.2f}")
    return 0


def cmd_list(a):
    s, lib = _lib()
    rows = [
        {"name": sp.name, "smiles": sp.smiles, "charge": sp.charge, "M": round(sp.molar_mass, 2),
         "conf": len(sp.conformers), "dH_fus": sp.dh_fus, "T_m": sp.t_fus, "source": sp.source}
        for sp in lib.all()
    ]
    print(pd.DataFrame(rows).to_string(index=False) if rows else "Library is empty.")
    return 0


def cmd_set(a):
    s, lib = _lib()
    fields = {k: v for k, v in (("dh_fus", a.dh_fus), ("t_fus", a.t_fus), ("dcp_fus", a.dcp_fus)) if v is not None}
    sp = lib.update(a.name, **fields)
    print(f"{sp.name}: dH_fus={sp.dh_fus} J/mol, T_m={sp.t_fus} K, dCp_fus={sp.dcp_fus}")
    return 0


def cmd_remove(a):
    s, lib = _lib()
    lib.remove(a.name)
    print(f"Removed {a.name}")
    return 0


def _resolve_systems(lib, specs, water):
    systems = [SolventSystem.parse(t) for t in specs] if specs else preset_systems()
    out = []
    for sy in systems:
        missing = [c for c in sy.species if not lib.has(c)]
        if missing:
            print(f"  skip {sy.name}: missing species {missing} (osnades add <name>)", file=sys.stderr)
            continue
        out.append(sy.with_water(water, lib) if water else sy)
    return out


def cmd_solubility(a):
    s, lib = _lib()
    eng = Engine(lib, a.parameterization)
    systems = _resolve_systems(lib, a.solvent, a.water)
    fus = {k: v for k, v in (("dh_fus", a.dh_fus), ("t_fus", a.t_fus), ("dcp_fus", a.dcp_fus)) if v is not None}
    frames = [eng.screen(a.solute, systems, T, method=a.method, **fus) for T in a.T]
    df = pd.concat(frames, ignore_index=True)
    cols = [c for c in ["solvent", "T", "x", "w", "g_per_kg", "ln_gamma", "ln_gamma_inf", "x_ideal", "note"] if c in df]
    with pd.option_context("display.float_format", "{:.5g}".format, "display.width", 200):
        print(f"Solute: {a.solute}   method: {a.method}\n")
        print(df[cols].to_string(index=False))
    if a.csv:
        df.to_csv(a.csv, index=False)
        print(f"\nSaved {a.csv}")
    return 0


def cmd_curves(a):
    from .plots import compute_curves, matplotlib_png

    s, lib = _lib()
    eng = Engine(lib, a.parameterization)
    systems = _resolve_systems(lib, a.solvent, a.water)
    sp = lib.get(a.solute)
    dh = a.dh_fus if a.dh_fus is not None else sp.dh_fus
    tm = a.t_fus if a.t_fus is not None else sp.t_fus
    dcp = a.dcp_fus if a.dcp_fus is not None else (sp.dcp_fus or 0.0)
    if dh is None or tm is None:
        raise ValueError(f"Set dH_fus and T_m for {a.solute} (osnades set) or pass --dh-fus/--t-fus")
    temps = np.linspace(a.tmin, a.tmax, a.n) + 273.15
    df = compute_curves(eng, a.solute, systems, temps, method=a.method, dh_fus=dh, t_fus=tm, dcp_fus=dcp,
                        progress=lambda f, n: print(f"  {f:5.0%} {n}", file=sys.stderr))
    out = Path(a.out or f"{a.solute}_solubility_vs_T")
    df.to_csv(out.with_suffix(".csv"), index=False)
    out.with_suffix(".png").write_bytes(matplotlib_png(df, a.solute, dh, tm, dcp))
    print(f"Saved {out.with_suffix('.csv')} and {out.with_suffix('.png')}")
    return 0


def cmd_binary(a):
    s, lib = _lib()
    df = Engine(lib, a.parameterization).binary(a.a, a.b, a.T, a.n)
    print(df.to_string(index=False, float_format="{:.5f}".format))
    return 0


def cmd_bootstrap(a):
    """Compute all preset species needed for the given solvent systems."""
    s, lib = _lib()
    systems = [SolventSystem.parse(t) for t in a.systems] if a.systems else preset_systems()
    needed = sorted({c for sy in systems for c in sy.species} | set(a.species or []))
    pre = presets()["species"]
    for name in needed:
        if lib.has(name):
            print(f"✓ {name}")
            continue
        if name not in pre:
            print(f"? {name}: no preset SMILES, add it manually", file=sys.stderr)
            continue
        meta = {k: pre[k] for k in ("dh_fus", "t_fus", "notes") if k in pre}
        lib.compute(name, smiles=pre[name]["smiles"], n_conformers=a.conformers, optimizer=a.optimizer, **meta)
    return 0


def cmd_app(a):
    app = Path(__file__).resolve().parent / "app.py"
    return subprocess.call([sys.executable, "-m", "streamlit", "run", str(app), "--server.port", str(a.port)])


def main(argv=None):
    p = argparse.ArgumentParser(prog="osnades", description="ORCA + openCOSMO-RS solubility and NADES screening")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("setup", help="configure ORCA backend, build Docker runtime, run a test")
    q.add_argument("--orca-home"); q.add_argument("--mode", choices=["auto", "native", "docker"])
    q.add_argument("--ncores", type=int); q.add_argument("--skip-build", action="store_true")
    q.add_argument("--no-test", action="store_true"); q.set_defaults(fn=cmd_setup)

    q = sub.add_parser("status", help="show settings and backend status")
    q.add_argument("--test", action="store_true", help="run a tiny ORCA job"); q.set_defaults(fn=cmd_status)

    q = sub.add_parser("add", help="compute a species with ORCA (xTB opt + COSMO-RS) and add it")
    q.add_argument("name"); q.add_argument("--smiles"); q.add_argument("--xyz")
    q.add_argument("--charge", type=int); q.add_argument("--mult", type=int)
    q.add_argument("--conformers", type=int, default=1)
    q.add_argument("--optimizer", choices=["xtb", "mmff", "none"], default="xtb")
    q.add_argument("--dh-fus", type=float, help="J/mol"); q.add_argument("--t-fus", type=float, help="K")
    q.add_argument("--overwrite", action="store_true"); q.set_defaults(fn=cmd_add)

    q = sub.add_parser("import", help="import existing .orcacosmo file(s) (one per conformer)")
    q.add_argument("name"); q.add_argument("files", nargs="+"); q.add_argument("--smiles")
    q.add_argument("--dh-fus", type=float); q.add_argument("--t-fus", type=float)
    q.add_argument("--overwrite", action="store_true"); q.set_defaults(fn=cmd_import)

    q = sub.add_parser("list", help="list library"); q.set_defaults(fn=cmd_list)
    q = sub.add_parser("set", help="set fusion data of a species")
    q.add_argument("name"); q.add_argument("--dh-fus", type=float); q.add_argument("--t-fus", type=float)
    q.add_argument("--dcp-fus", type=float); q.set_defaults(fn=cmd_set)
    q = sub.add_parser("remove", help="remove a species"); q.add_argument("name"); q.set_defaults(fn=cmd_remove)

    q = sub.add_parser("bootstrap", help="compute preset species for solvent systems (default: all presets)")
    q.add_argument("systems", nargs="*"); q.add_argument("--species", nargs="*")
    q.add_argument("--conformers", type=int, default=1)
    q.add_argument("--optimizer", choices=["xtb", "mmff", "none"], default="xtb"); q.set_defaults(fn=cmd_bootstrap)

    q = sub.add_parser("solubility", help="solid solubility of a solute in solvents / NADES")
    q.add_argument("solute")
    q.add_argument("-s", "--solvent", action="append",
                   help='preset name, species name, or "a:b:c=1:1:2" (repeatable; default: all presets available)')
    q.add_argument("-T", type=float, nargs="+", default=[298.15])
    q.add_argument("--method", choices=["iterative", "infinite_dilution"], default="iterative")
    q.add_argument("--water", type=float, default=0.0, help="add water, wt%% of solvent")
    q.add_argument("--dh-fus", type=float); q.add_argument("--t-fus", type=float); q.add_argument("--dcp-fus", type=float)
    q.add_argument("--parameterization", default="openCOSMORS24a"); q.add_argument("--csv")
    q.set_defaults(fn=cmd_solubility)

    q = sub.add_parser("curves", help="solubility vs temperature curves (CSV + grouped PNG)")
    q.add_argument("solute"); q.add_argument("-s", "--solvent", action="append")
    q.add_argument("--tmin", type=float, default=0.0, help="°C"); q.add_argument("--tmax", type=float, default=80.0, help="°C")
    q.add_argument("-n", type=int, default=17)
    q.add_argument("--method", choices=["iterative", "infinite_dilution"], default="iterative")
    q.add_argument("--water", type=float, default=0.0)
    q.add_argument("--dh-fus", type=float); q.add_argument("--t-fus", type=float); q.add_argument("--dcp-fus", type=float)
    q.add_argument("--parameterization", default="openCOSMORS24a"); q.add_argument("-o", "--out")
    q.set_defaults(fn=cmd_curves)

    q = sub.add_parser("binary", help="activity coefficients across a binary mixture")
    q.add_argument("a"); q.add_argument("b"); q.add_argument("-T", type=float, default=298.15)
    q.add_argument("-n", type=int, default=21); q.add_argument("--parameterization", default="openCOSMORS24a")
    q.set_defaults(fn=cmd_binary)

    q = sub.add_parser("app", help="launch the web application"); q.add_argument("--port", type=int, default=8501)
    q.set_defaults(fn=cmd_app)

    a = p.parse_args(argv)
    try:
        return a.fn(a) or 0
    except (KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
