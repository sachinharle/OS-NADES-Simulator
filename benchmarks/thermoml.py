"""Harvest experimental solid–liquid solubility data from the NIST/TRC ThermoML Archive.

Every ThermoML record that contains one of the library's solid compounds is found through the
ThermoML-API (search by standard InChIKey), downloaded once into ``benchmarks/cache/`` and parsed.
A data point is kept when

* the property is the composition of a component in the liquid phase (mole fraction, mass fraction
  or molality) that coexists with that same component's crystal (= its solubility / liquidus), and
* every other component is in the library (choline chloride maps to choline⁺ + Cl⁻), and
* the solvent composition is fully specified (pure solvent, or all but one solvent component given
  as a mole or mass fraction, on a solute-free or total basis).

The result is one row per point with the full mole-fraction composition on a formula-unit basis
(a salt counts once), so it compares directly with the published numbers.
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

API = "https://trc.nist.gov/ThermoML-API/objects"
XML = "https://trc.nist.gov/ThermoML/{doi}.xml"
NS = {"t": "http://www.iupac.org/namespaces/ThermoML"}
CACHE = Path(__file__).parent / "cache"
UA = {"User-Agent": "osnades-benchmark/0.1 (research use)"}


# --------------------------------------------------------------------------- download
def _get(url: str, timeout: float = 60.0) -> bytes:
    from osnades.fusion import _ssl_context

    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout, context=_ssl_context()) as r:
        return r.read()


def search_dois(inchikey: str) -> list[str]:
    """DOIs of all ThermoML records that contain the compound."""
    q = f'type:TRCTml4 AND /Compound/_/sStandardInChIKey:"{inchikey}"'
    out, page = [], 0
    while True:
        url = API + "?" + urllib.parse.urlencode({"query": q, "ids": "true", "pageSize": 200, "pageNum": page})
        d = json.loads(_get(url))
        out += [r.split("thermoml/", 1)[1] for r in d["results"]]
        if len(out) >= d["size"] or not d["results"]:
            return out
        page += 1


def fetch(doi: str) -> Path:
    """Download (once) the ThermoML file of a record."""
    p = CACHE / (doi.replace("/", "_") + ".xml")
    if not p.exists() or p.stat().st_size == 0:
        CACHE.mkdir(parents=True, exist_ok=True)
        p.write_bytes(_get(XML.format(doi=doi)))
        time.sleep(0.3)  # be gentle with the server
    return p


# --------------------------------------------------------------------------- parsing
def _t(el, path: str, default=None):
    x = el.find(path, NS)
    return x.text.strip() if x is not None and x.text else default


def _num(el, path: str) -> float | None:
    v = _t(el, path)
    return float(v) if v is not None else None


@dataclass
class Compound:
    org: int
    inchikey: str
    names: list[str]
    formula: str


def _compounds(root) -> dict[int, Compound]:
    out = {}
    for c in root.findall("t:Compound", NS):
        org = int(_t(c, "t:RegNum/t:nOrgNum", "0"))
        out[org] = Compound(org, _t(c, "t:sStandardInChIKey", "") or "",
                            [n.text for n in c.findall("t:sCommonName", NS) if n.text],
                            _t(c, "t:sFormulaMolec", "") or "")
    return out


def _comp_spec(vid) -> tuple[str, int | None] | None:
    """(kind, component) of a composition variable/constraint, e.g. ('solvent:x', 3) or ('x', 1)."""
    ct = vid.find("t:VariableType", NS) if vid.tag.endswith("VariableID") else vid.find("t:ConstraintType", NS)
    if ct is None:
        return None
    reg = vid.find("t:RegNum/t:nOrgNum", NS)
    org = int(reg.text) if reg is not None else None
    comp = ct.find("t:eComponentComposition", NS)
    solv = ct.find("t:eSolventComposition", NS)
    if comp is not None:
        kind = {"Mole fraction": "x", "Mass fraction": "w", "Molality, mol/kg": "m"}.get(comp.text.strip())
        return (kind, org) if kind else ("unsupported", org)
    if solv is not None:
        kind = {"Solvent: Mole fraction": "solvent:x", "Solvent: Mass fraction": "solvent:w"}.get(solv.text.strip())
        return (kind, org) if kind else ("unsupported", org)
    return None


def parse(path: Path, mapping) -> tuple[list[dict], dict]:
    """Solubility points of one ThermoML file. ``mapping(compound) -> (name, molar_mass) | None``."""
    root = ET.parse(path).getroot()
    cit = root.find("t:Citation", NS)
    doi = _t(cit, "t:sDOI", path.stem)
    title = _t(cit, "t:sTitle", "")
    year = _t(cit, "t:yrPubYr", "")
    comps = _compounds(root)
    rows, skipped = [], {}

    def skip(why):
        skipped[why] = skipped.get(why, 0) + 1

    for ds in root.findall("t:PureOrMixtureData", NS):
        members = [int(c.text) for c in ds.findall("t:Component/t:RegNum/t:nOrgNum", NS)]
        if len(members) < 2:
            continue
        crystals = {int(r.text) for ph in ds.findall("t:PhaseID", NS)
                    if (_t(ph, "t:ePhase", "") or "").startswith("Crystal")
                    for r in ph.findall("t:RegNum/t:nOrgNum", NS)}
        phases = [(_t(ph, "t:ePhase", "") or "") for ph in ds.findall("t:PhaseID", NS)]
        props = {}
        for pr in ds.findall("t:Property", NS):
            name = _t(pr, ".//t:ePropName", "")
            org = pr.find("t:Property-MethodID/t:RegNum/t:nOrgNum", NS)
            phase = _t(pr, "t:PropPhaseID/t:ePropPhase", "")
            kind = {"Mole fraction": "x", "Mass fraction": "w", "Molality, mol/kg": "m"}.get(name)
            if kind and org is not None and phase == "Liquid" and int(org.text) in crystals:
                props[int(_t(pr, "t:nPropNumber"))] = (kind, int(org.text))
        if not props:
            continue
        # every crystal phase must be a single pure component (no co-crystals, solvates, hydrates as
        # separate compounds are fine – they are other compounds and will not map)
        if sum(p.startswith("Crystal") for p in phases) != 1:
            skip("several solid phases in the data set")
            continue
        names = {}
        for o in members:
            m = mapping(comps.get(o)) if o in comps else None
            if m is None:
                break
            names[o] = m
        else:
            m = None
        if len(names) != len(members):
            skip("component not in the library")
            continue
        # composition specs: constraints (fixed) and variables (per point)
        fixed, var_kind, T_fixed = {}, {}, None
        bad = False
        for c in ds.findall("t:Constraint", NS):
            cid = c.find("t:ConstraintID", NS)
            if cid.find(".//t:eTemperature", NS) is not None:
                T_fixed = _num(c, "t:nConstraintValue")
                continue
            sp = _comp_spec(cid)
            if sp is None:
                continue
            if sp[0] == "unsupported":
                bad = True
            fixed[sp] = _num(c, "t:nConstraintValue")
        T_var = None
        for v in ds.findall("t:Variable", NS):
            vid = v.find("t:VariableID", NS)
            n = int(_t(v, "t:nVarNumber"))
            if vid.find(".//t:eTemperature", NS) is not None:
                T_var = n
                continue
            sp = _comp_spec(vid)
            if sp is None:
                continue
            if sp[0] == "unsupported":
                bad = True
            var_kind[n] = sp
        if bad:
            skip("unsupported composition variable")
            continue
        for pn, (kind, solute) in props.items():
            solvent = [o for o in members if o != solute]
            for nv in ds.findall("t:NumValues", NS):
                vals = {int(_t(x, "t:nVarNumber")): float(_t(x, "t:nVarValue"))
                        for x in nv.findall("t:VariableValue", NS)}
                pv = {int(_t(x, "t:nPropNumber")): (float(_t(x, "t:nPropValue")),
                                                    _num(x, ".//t:nCombExpandUncertValue"))
                      for x in nv.findall("t:PropertyValue", NS)}
                if pn not in pv:
                    continue
                T = vals.get(T_var) if T_var else T_fixed
                if T is None:
                    skip("no temperature")
                    continue
                spec = dict(fixed)
                spec.update({var_kind[k]: v for k, v in vals.items() if k in var_kind})
                x = _resolve(kind, pv[pn][0], solute, solvent, spec, names)
                if x is None:
                    skip("solvent composition not fully specified")
                    continue
                rows.append({"doi": doi, "year": year, "title": title, "dataset": int(_t(ds, "t:nPureOrMixtureDataNumber")),
                             "solute": names[solute][0], "T": T, "x_exp": x[names[solute][0]],
                             "u_exp": pv[pn][1], "reported_as": kind, "n_components": len(members),
                             "composition": json.dumps({k: round(v, 8) for k, v in x.items()})})
    return rows, skipped


def _resolve(kind, value, solute, solvent, spec, names) -> dict[str, float] | None:
    """Full mole-fraction composition (formula units) from the reported solute amount + solvent specs."""
    M = {o: names[o][1] for o in [solute] + solvent}
    # solvent composition on a solute-free mole basis
    xs = {}
    if len(solvent) == 1:
        xs = {solvent[0]: 1.0}
    else:
        known_x = {o: v for (k, o), v in spec.items() if k == "solvent:x" and o in solvent}
        known_w = {o: v for (k, o), v in spec.items() if k == "solvent:w" and o in solvent}
        tot_x = {o: v for (k, o), v in spec.items() if k == "x" and o in solvent}
        tot_w = {o: v for (k, o), v in spec.items() if k == "w" and o in solvent}
        if len(known_x) == len(solvent) - 1:
            rest = [o for o in solvent if o not in known_x][0]
            xs = {**known_x, rest: 1 - sum(known_x.values())}
        elif len(known_w) == len(solvent) - 1:
            rest = [o for o in solvent if o not in known_w][0]
            w = {**known_w, rest: 1 - sum(known_w.values())}
            n = {o: w[o] / M[o] for o in solvent}
            xs = {o: n[o] / sum(n.values()) for o in solvent}
        elif len(tot_x) == len(solvent) - 1 and kind == "x":
            rest = [o for o in solvent if o not in tot_x][0]
            xt = {**tot_x, rest: 1 - value - sum(tot_x.values())}
            xs = {o: xt[o] / (1 - value) for o in solvent}
        elif len(tot_w) == len(solvent) - 1 and kind == "w":
            rest = [o for o in solvent if o not in tot_w][0]
            wt = {**tot_w, rest: 1 - value - sum(tot_w.values())}
            n = {o: wt[o] / M[o] for o in solvent}
            xs = {o: n[o] / sum(n.values()) for o in solvent}
        else:
            return None
    if any(v < -1e-9 for v in xs.values()):
        return None
    M_solv = sum(xs[o] * M[o] for o in solvent)
    if kind == "x":
        x = value
    elif kind == "w":
        x = (value / M[solute]) / (value / M[solute] + (1 - value) / M_solv)
    else:  # molality (mol/kg solvent)
        r = value * M_solv / 1000.0
        x = r / (1 + r)
    if not 0 < x < 1:
        return None
    out = {names[solute][0]: x}
    for o in solvent:
        out[names[o][0]] = out.get(names[o][0], 0.0) + (1 - x) * xs[o]
    return out


def harvest(inchikeys: dict[str, str], mapping, log=print) -> tuple[pd.DataFrame, dict]:
    """All solubility points for the given solutes ({name: InChIKey})."""
    dois: set[str] = set()
    for name, key in inchikeys.items():
        try:
            found = search_dois(key)
            flat = f"{first_block(key)}-UHFFFAOYSA-N"  # stereo-unspecified entry of the same compound
            if flat != key:
                found = list(dict.fromkeys(found + search_dois(flat)))
        except Exception as exc:  # noqa: BLE001
            log(f"  search failed for {name}: {exc}")
            continue
        log(f"  {name}: {len(found)} ThermoML records")
        dois |= set(found)
    rows, skipped = [], {}
    for i, doi in enumerate(sorted(dois)):
        try:
            p = fetch(doi)
            r, s = parse(p, mapping)
        except Exception as exc:  # noqa: BLE001
            log(f"  {doi}: {exc}")
            continue
        rows += r
        for k, v in s.items():
            skipped[k] = skipped.get(k, 0) + v
    df = pd.DataFrame(rows)
    if len(df):
        df = df[df["solute"].isin(inchikeys)].drop_duplicates(["doi", "solute", "T", "x_exp", "composition"])
    return df.reset_index(drop=True), {"records": len(dois), "skipped": skipped}


def first_block(key: str) -> str:
    return (key or "").split("-")[0]


def norm_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())
