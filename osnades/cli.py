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
    if "dh_fus" in meta or "t_fus" in meta:
        meta["fusion_source"] = "entered by user" if (a.dh_fus or a.t_fus) else pre.get("notes", "preset")
    lib.compute(a.name, smiles=smiles, xyz=xyz, charge=a.charge, multiplicity=a.mult,
                n_conformers=a.conformers, optimizer=a.optimizer, overwrite=a.overwrite,
                lookup_fusion=not a.no_fusion_lookup, **meta)
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
    if a.dh_fus is not None or a.t_fus is not None:
        fields["fusion_source"] = "entered by user"
    sp = lib.update(a.name, **fields)
    print(f"{sp.name}: dH_fus={sp.dh_fus} J/mol, T_m={sp.t_fus} K, dCp_fus={sp.dcp_fus}")
    return 0


def cmd_fusion(a):
    """Look up published ΔH_fus / T_m and (with --apply) store the suggested pair."""
    from . import fusion

    s, lib = _lib()
    def wanted(sp):
        src = sp.fusion_source or ""
        return not (sp.dh_fus and sp.t_fus) or (a.estimates and ("Walden" in src or "estimate" in src))

    names = [sp.name for sp in lib.all() if sp.charge == 0 and wanted(sp)] if a.all else a.names
    if not names:
        print("Nothing to look up (give species names or --all).")
        return 0
    for n in names:
        sp = lib.get(n)
        recs = fusion.lookup(sp.smiles, n)
        pick = fusion.choose(recs, n)
        print(f"== {n}: {len(recs)} published value(s)")
        for r in recs:
            print(("  * " if pick and r.source == pick.source and r.url == pick.url else "    ") + r.label)
        if pick is None:
            print("  no melting data found – set it with `osnades set`")
            continue
        print(f"  suggested: {pick.label}")
        if pick.estimated:
            print("  ⚠ ΔH_fus is a Walden's-rule estimate – replace it with a measured value if you find one")
        if a.apply:
            lib.update(n, dh_fus=pick.dh_fus, t_fus=pick.t_fus, fusion_source=pick.note)
            print("  saved to library")
    return 0


def cmd_remove(a):
    s, lib = _lib()
    lib.remove(a.name)
    print(f"Removed {a.name}")
    return 0


def _resolve_systems(lib, specs, water):
    systems = [SolventSystem.parse(t, lib) for t in specs] if specs else preset_systems(lib)
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


def cmd_design(a):
    """Find the component mole fractions of a multicomponent NADES that dissolve the most solute."""
    from . import designer as D

    s, lib = _lib()
    comps = []
    for spec in a.component:
        parts = spec.rsplit(":", 2) if spec.count(":") >= 2 else [spec]
        name = parts[0]
        lo, hi = (float(parts[1]), float(parts[2])) if len(parts) == 3 else \
            ((0.1, 0.5) if len(D.species_of(name)) > 1 else (0.0, 0.3 if name == "water" else 1.0))
        comps.append(D.make_component(name, lib, lo, hi))
    miss = D.missing_species(a.solute, comps, lib)
    if miss:
        print("Missing σ-surfaces – compute them first: " + ", ".join(f"osnades add {m}" for m in miss))
        return 1
    sp = lib.get(a.solute)
    dh, tm = a.dh_fus or sp.dh_fus, a.t_fus or sp.t_fus
    if not (dh and tm):
        print(f"ΔH_fus / T_m of {a.solute} unknown: run `osnades fusion {a.solute} --apply` or pass --dh-fus/--t-fus")
        return 1
    temps = np.linspace(a.tmin, a.tmax, a.nT) + 273.15 if a.nT > 1 else np.array([a.tmin + 273.15])
    print(f"{a.solute} in " + ", ".join(f"{c.name} [{c.lo:g}–{c.hi:g}]" for c in comps)
          + f" · {a.tmin:g}–{a.tmax:g} °C · goal: {D.OBJECTIVES[a.objective]}")
    des = D.Designer(Engine(lib, a.parameterization), a.solute, comps, dh, tm, sp.dcp_fus or 0.0)
    res = des.run(temps, step=a.step, objective=a.objective, require_liquid=not a.allow_solid, top_n=a.top,
                  local_rounds=a.rounds, progress=lambda f, m: print(f"  [{100 * f:5.1f}%] {m}", flush=True))
    for n in res.notes:
        print("note:", n)
    cols = ["rank", "name", "liquid_all_T"] + [c for c in res.best if c.startswith("x@") or c.startswith("g/kg@")]
    print(res.best[cols].to_string(index=False))
    if a.out:
        res.best.to_csv(f"{a.out}_best.csv", index=False)
        res.screen.to_csv(f"{a.out}_screen.csv", index=False)
        print(f"wrote {a.out}_best.csv and {a.out}_screen.csv")
    if a.save_as:
        lib.add_solvent(a.save_as, res.best_system(0).components, overwrite=True)
        print(f"saved the best composition as solvent {a.save_as!r}")
    return 0


def cmd_ternary(a, lib):
    from .eutectic import ratio_label, resolve_component, ternary_diagram, ternary_png

    cs = [resolve_component(n, lib, getattr(a, f"dh_{k}"), getattr(a, f"tm_{k}"), getattr(a, f"dcp_{k}"))
          for n, k in zip((a.a, a.b, a.c), "abc")]
    for c in cs:
        print(f"{c.name}: ΔH_fus = {c.dh_fus / 1000:.2f} kJ/mol, T_m = {c.t_fus:.1f} K ({c.source})")
    res = ternary_diagram(Engine(lib, a.parameterization), cs, step=a.step, T_lo=a.T_lo)
    names = ":".join(res.names)
    print()
    if np.isfinite(res.T_e):
        print(f"Ternary eutectic (openCOSMO-RS): {names} = {ratio_label(res.z_e)} "
              f"(x = {', '.join(f'{z:.3f}' for z in res.z_e)}), T = {res.T_e:.1f} K ({res.T_e - 273.15:.1f} °C)")
    elif np.isfinite(res.T_e_upper):
        print(f"Ternary eutectic (openCOSMO-RS): below {res.T_e_upper:.0f} K around {names} = {ratio_label(res.z_e)}")
    print(f"Ternary eutectic (ideal):        {names} = {ratio_label(res.z_e_ideal)}, T = {res.T_e_ideal:.1f} K "
          f"({res.T_e_ideal - 273.15:.1f} °C)")
    if a.T is not None:
        print(f"Liquid at {a.T:.2f} K: {res.liquid_fraction(a.T):.0%} of the composition triangle")
    for n in res.notes:
        print(f"note: {n}")
    if a.out:
        out = Path(a.out)
        res.data.to_csv(out.with_suffix(".csv"), index=False)
        out.with_suffix(".png").write_bytes(ternary_png(res, a.T))
        print(f"Saved {out.with_suffix('.csv')} and {out.with_suffix('.png')}")
    return 0


def cmd_eutectic(a):
    from .eutectic import matplotlib_png, phase_diagram, resolve_component

    s, lib = _lib()
    if a.c:
        return cmd_ternary(a, lib)
    A = resolve_component(a.a, lib, a.dh_a, a.tm_a, a.dcp_a)
    B = resolve_component(a.b, lib, a.dh_b, a.tm_b, a.dcp_b)
    for c in (A, B):
        print(f"{c.name}: ΔH_fus = {c.dh_fus / 1000:.2f} kJ/mol, T_m = {c.t_fus:.1f} K ({c.source})")
    res = phase_diagram(Engine(lib, a.parameterization), A, B, n=a.n)
    print()
    if np.isfinite(res.T_e):
        print(f"Eutectic (openCOSMO-RS): x_{A.name} = {res.x_e:.3f}, T = {res.T_e:.1f} K ({res.T_e - 273.15:.1f} °C)")
    elif np.isfinite(res.T_e_upper):
        print(f"Eutectic (openCOSMO-RS): below {res.T_e_upper:.0f} K near x_{A.name} ≈ {res.x_e:.2f}")
    print(f"Eutectic (ideal):        x_{A.name} = {res.x_e_ideal:.3f}, T = {res.T_e_ideal:.1f} K "
          f"({res.T_e_ideal - 273.15:.1f} °C)")
    if np.isfinite(res.depression):
        print(f"Depression below ideal: {res.depression:.1f} K")
    if a.T is not None:
        w = res.liquid_window(a.T)
        print(f"Liquid at {a.T:.2f} K for x_{A.name} = " + (f"{w[0]:.3f}–{w[1]:.3f}" if w else "none"))
    for n in res.notes:
        print(f"note: {n}")
    if a.out:
        out = Path(a.out)
        res.data.to_csv(out.with_suffix(".csv"), index=False)
        out.with_suffix(".png").write_bytes(matplotlib_png(res, a.T))
        print(f"Saved {out.with_suffix('.csv')} and {out.with_suffix('.png')}")
    return 0


def cmd_solvent(a):
    s, lib = _lib()
    if a.action == "list":
        user = lib.user_solvents()
        if not user:
            print("No user solvents. Add one with: osnades solvent add <name> [\"a:b=1:2\"]")
        for k, v in user.items():
            missing = [c for c in v if not lib.has(c)]
            print(f"{k}: " + ", ".join(f"{c} {r:g}" for c, r in v.items()) + (f"   (missing: {missing})" if missing else ""))
        return 0
    if not a.name:
        raise ValueError("Give a solvent name")
    if a.action == "remove":
        lib.remove_solvent(a.name)
        print(f"Removed solvent {a.name}")
        return 0
    comps = SolventSystem.parse(a.spec).components if a.spec else None
    lib.add_solvent(a.name, comps, overwrite=a.overwrite)
    missing = [c for c in (comps or {a.name: 1}) if not lib.has(c)]
    print(f"Added solvent {a.name}" + (f" (compute first: {', '.join(missing)})" if missing else ""))
    return 0


def cmd_bootstrap(a):
    """Compute all preset species needed for the given solvent systems."""
    s, lib = _lib()
    systems = [SolventSystem.parse(t, lib) for t in a.systems] if a.systems else preset_systems()
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
    q.add_argument("--no-fusion-lookup", action="store_true",
                   help="do not fetch published ΔH_fus/T_m for fields left empty")
    q.add_argument("--overwrite", action="store_true"); q.set_defaults(fn=cmd_add)

    q = sub.add_parser("import", help="import existing .orcacosmo file(s) (one per conformer)")
    q.add_argument("name"); q.add_argument("files", nargs="+"); q.add_argument("--smiles")
    q.add_argument("--dh-fus", type=float); q.add_argument("--t-fus", type=float)
    q.add_argument("--overwrite", action="store_true"); q.set_defaults(fn=cmd_import)

    q = sub.add_parser("list", help="list library"); q.set_defaults(fn=cmd_list)
    q = sub.add_parser("set", help="set fusion data of a species")
    q.add_argument("name"); q.add_argument("--dh-fus", type=float); q.add_argument("--t-fus", type=float)
    q.add_argument("--dcp-fus", type=float); q.set_defaults(fn=cmd_set)
    q = sub.add_parser("fusion", help="look up published fusion data (ΔH_fus, T_m) of species")
    q.add_argument("names", nargs="*"); q.add_argument("--all", action="store_true", help="all species without data")
    q.add_argument("--estimates", action="store_true", help="with --all: also re-check estimated values")
    q.add_argument("--apply", action="store_true", help="store the suggested values"); q.set_defaults(fn=cmd_fusion)
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

    q = sub.add_parser("solvent", help="manage your own solvents (pure species or mixtures)")
    q.add_argument("action", choices=["add", "list", "remove"]); q.add_argument("name", nargs="?")
    q.add_argument("spec", nargs="?", help='species name or "a:b:c=1:1:2" (default: the species called NAME)')
    q.add_argument("--overwrite", action="store_true"); q.set_defaults(fn=cmd_solvent)

    q = sub.add_parser("eutectic", help="solid–liquid phase diagram / eutectic point of 2 or 3 components")
    q.add_argument("a", help='component A: species, preset salt ("choline chloride") or "cation+anion"')
    q.add_argument("b", help="component B")
    q.add_argument("c", nargs="?", help="optional component C: computes the ternary liquidus surface")
    q.add_argument("-n", type=int, default=41, help="composition grid points (binary)")
    q.add_argument("--step", type=float, default=0.1, help="composition grid step (ternary)")
    q.add_argument("--T-lo", type=float, help="lowest temperature checked (K, ternary)")
    q.add_argument("-T", type=float, help="operating temperature (K): report the liquid composition window")
    for c in ("a", "b", "c"):
        q.add_argument(f"--dh-{c}", type=float, help=f"ΔH_fus of {c.upper()} (J/mol)")
        q.add_argument(f"--tm-{c}", type=float, help=f"T_m of {c.upper()} (K)")
        q.add_argument(f"--dcp-{c}", type=float, help=f"ΔCp_fus of {c.upper()} (J/mol/K)")
    q.add_argument("--parameterization", default="openCOSMORS24a"); q.add_argument("-o", "--out")
    q.set_defaults(fn=cmd_eutectic)

    q = sub.add_parser("design", help="optimise the mole fractions of a multicomponent NADES for a solute")
    q.add_argument("solute")
    q.add_argument("-c", "--component", action="append", required=True,
                   help='component[:min:max], e.g. "choline chloride:0.1:0.5", glycerol, "water:0:0.3" (repeat)')
    q.add_argument("--tmin", type=float, default=15.0, help="°C"); q.add_argument("--tmax", type=float, default=45.0)
    q.add_argument("--nT", type=int, default=3, help="temperature points")
    q.add_argument("--objective", choices=["mean_x", "min_T", "g_per_kg", "cooling"], default="mean_x")
    q.add_argument("--step", type=float, default=0.1, help="composition grid step")
    q.add_argument("--rounds", type=int, default=2, help="local refinement rounds")
    q.add_argument("--top", type=int, default=6); q.add_argument("--allow-solid", action="store_true",
                   help="do not require the solvent to be liquid over the range")
    q.add_argument("--dh-fus", type=float); q.add_argument("--t-fus", type=float)
    q.add_argument("--save-as", help="register the best composition as a solvent")
    q.add_argument("--parameterization", default="openCOSMORS24a"); q.add_argument("-o", "--out")
    q.set_defaults(fn=cmd_design)

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
