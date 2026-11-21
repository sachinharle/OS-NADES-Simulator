"""Accuracy benchmark of OS-NADES against experimental solubility data (NIST/TRC ThermoML).

    .venv/bin/python -m benchmarks.accuracy harvest   # find + parse the data → benchmarks/data/
    .venv/bin/python -m benchmarks.accuracy predict   # openCOSMO-RS + ideal predictions for every point
    .venv/bin/python -m benchmarks.accuracy report    # statistics, plots and benchmarks/REPORT.md

Predictions use exactly what the app uses: the library σ-surfaces, the library fusion data
(ΔH_fus, T_m, ΔCp_fus) and ``Engine.solubility`` (iterative SLE, openCOSMORS24a). The tool counts
solvent species individually (choline⁺ and Cl⁻ separately); predictions are converted to formula
units so that they compare with the published mole fractions.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from benchmarks import thermoml as TM

HERE = Path(__file__).parent
DATA = HERE / "data"
RAW = DATA / "thermoml_solubility.csv"
PRED = DATA / "predictions.csv"
SALTS = {"choline chloride": ("C[N+](C)(C)CCO.[Cl-]", {"choline": 1, "chloride": 1})}


def _lib():
    from osnades.cli import _lib as L

    return L()[1]


def _key(smiles: str) -> str:
    from rdkit import Chem

    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(m) if m is not None else ""


def mapping_table(lib) -> dict[str, tuple[str, float]]:
    """InChIKey first block → (library name, molar mass of one formula unit)."""
    out = {}
    for sp in lib.all():
        if sp.smiles and sp.charge == 0:
            out[TM.first_block(_key(sp.smiles))] = (sp.name, sp.molar_mass)
    for salt, (smi, ions) in SALTS.items():
        if all(lib.has(i) for i in ions):
            out[TM.first_block(_key(smi))] = (salt, sum(lib.get(i).molar_mass * n for i, n in ions.items()))
    return out


def solutes(lib) -> dict[str, str]:
    """Library compounds that are solid near room temperature (T_m ≥ 285 K) and have fusion data.

    Low-melting liquids (water, alcohols…) are excluded: they appear in thousands of ThermoML records
    as solvents, and their freezing-point data are not what the tool is for."""
    return {sp.name: _key(sp.smiles) for sp in lib.all()
            if sp.charge == 0 and sp.smiles and sp.dh_fus and sp.t_fus and sp.t_fus >= 285.0}


# --------------------------------------------------------------------------- harvest
def cmd_harvest():
    lib = _lib()
    table = mapping_table(lib)
    print(f"{len(table)} library compounds can be matched; searching ThermoML for {len(solutes(lib))} solutes")

    def mapping(c):
        return table.get(TM.first_block(c.inchikey)) if c is not None else None

    df, info = TM.harvest(solutes(lib), mapping)
    df, qc = quality_control(df)
    info["quality_control"] = qc
    DATA.mkdir(exist_ok=True)
    df.to_csv(RAW, index=False)
    (DATA / "harvest_info.json").write_text(json.dumps(info, indent=1))
    print(f"{len(df)} points from {df['doi'].nunique() if len(df) else 0} papers "
          f"({info['records']} ThermoML records searched); skipped: {info['skipped']}")
    print(f"quality control: {qc['duplicates_removed']} duplicates, {qc['points_removed_inconsistent']} points in "
          f"{len(qc['inconsistent_sets'])} inconsistent sets removed")
    for line in qc["inconsistent_sets"]:
        print("   ", line)
    if len(df):
        print(df.groupby("solute").agg(points=("x_exp", "size"), papers=("doi", "nunique")).to_string())


def quality_control(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Drop republished duplicates and mixed-solvent series that contradict pure-solvent data.

    A mixed-solvent data set is inconsistent when its end point with only one solvent component
    differs from measured pure-solvent solubility (any paper, same T ± 0.6 K) by more than ×1.5."""
    n0 = len(df)
    df = df.drop_duplicates(["solute", "composition", "T", "x_exp"]).reset_index(drop=True)
    qc = {"duplicates_removed": n0 - len(df), "inconsistent_sets": []}
    comp = [json.loads(c) for c in df.composition]
    solv = [tuple(sorted(k for k, v in c.items() if k != s and v > 1e-9)) for c, s in zip(comp, df.solute)]
    df["solvset"] = solv
    pure = df[[len(x) == 1 and n == 2 for x, n in zip(solv, df.n_components)]]
    drop = set()
    for (doi, ds), g in df[df.n_components >= 3].groupby(["doi", "dataset"]):
        ends = g[[len(x) == 1 for x in g.solvset]]
        for r in ends.itertuples():
            ref = pure[(pure.solute == r.solute) & (pure.solvset == r.solvset) & ((pure["T"] - r.T).abs() < 0.6)]
            if len(ref) and abs(np.log10(r.x_exp / ref.x_exp.median())) > np.log10(1.5):
                drop.add((doi, ds))
                qc["inconsistent_sets"].append(f"{doi} set {ds} ({r.solute} in {'+'.join(r.solvset)} end point "
                                               f"{r.x_exp:.3g} vs pure-solvent {ref.x_exp.median():.3g} at {r.T:.1f} K)")
                break
    keep = [(d, s) not in drop for d, s in zip(df.doi, df.dataset)]
    qc["points_removed_inconsistent"] = int(len(df) - sum(keep))
    return df[keep].drop(columns="solvset").reset_index(drop=True), qc


# --------------------------------------------------------------------------- predict
def solvent_system(comp: dict[str, float], solute: str):
    from osnades.thermo import SolventSystem

    sp: dict[str, float] = {}
    for name, x in comp.items():
        if name == solute:
            continue
        for s, nu in SALTS.get(name, (None, {name: 1}))[1].items():
            sp[s] = sp.get(s, 0.0) + x * nu
    tot = sum(sp.values())
    return SolventSystem("+".join(sorted(sp)), {k: v / tot for k, v in sp.items() if v > 0})


def species_per_unit(comp: dict[str, float], solute: str) -> float:
    """Moles of COSMO species per mole of solvent formula units (1 + salt fraction)."""
    xs = {k: v for k, v in comp.items() if k != solute}
    tot = sum(xs.values())
    return sum(v / tot * sum(SALTS.get(k, (None, {k: 1}))[1].values()) for k, v in xs.items())


def category(row, lib) -> str:
    comp = json.loads(row["composition"])
    solv = [k for k in comp if k != row["solute"]]
    if any(k in SALTS for k in solv):
        return "NADES / ChCl mixtures"
    if len(solv) == 1:
        s = solv[0]
        sp = lib.get(s)
        if sp.t_fus and sp.t_fus > row["T"] + 0.5:
            return "binary SLE (both solid-formers)"
        return "water" if s == "water" else "pure organic solvent"
    return "aqueous mixtures" if "water" in solv else "organic mixtures"


def cmd_predict():
    from osnades.thermo import Engine

    lib = _lib()
    eng = Engine(lib, "openCOSMORS24a")
    df = pd.read_csv(RAW)
    out, cache = [], {}
    t0 = time.time()
    for i, r in df.iterrows():
        comp = json.loads(r["composition"])
        sysm = solvent_system(comp, r["solute"])
        key = (r["solute"], tuple(sorted((k, round(v, 6)) for k, v in sysm.components.items())), round(r["T"], 2))
        if key not in cache:
            t = time.time()
            try:
                res = eng.solubility(r["solute"], sysm, float(r["T"]))
                cache[key] = (res.x, res.x_ideal, res.note, time.time() - t)
            except Exception as exc:  # noqa: BLE001
                cache[key] = (np.nan, np.nan, f"error: {exc}", time.time() - t)
        x_sp, x_id_sp, note, dt = cache[key]
        k = species_per_unit(comp, r["solute"])
        to_fu = lambda x: x / (x + (1 - x) / k) if np.isfinite(x) else np.nan  # noqa: E731
        out.append({"x_pred": to_fu(x_sp), "x_ideal": to_fu(x_id_sp) if k != 1 else x_id_sp,
                    "note": note, "seconds": dt, "category": category(r, lib)})
        if (i + 1) % 200 == 0:
            print(f"  {i + 1}/{len(df)} points ({time.time() - t0:.0f} s)", flush=True)
    res = pd.concat([df, pd.DataFrame(out)], axis=1)
    res.to_csv(PRED, index=False)
    print(f"predicted {len(res)} points ({len(cache)} unique) in {time.time() - t0:.0f} s → {PRED}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd in ("harvest", "all"):
        cmd_harvest()
    if cmd in ("predict", "all"):
        cmd_predict()
    if cmd in ("report", "all"):
        from benchmarks.report import cmd_report

        cmd_report()
