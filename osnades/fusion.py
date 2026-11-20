"""Look up experimental fusion data (ΔH_fus, T_m) for a compound.

Source: Cheméo (https://www.chemeo.com), which collects measured values from the NIST WebBook and
the literature (each value carries its reference). Rows that Cheméo marks as *Calculated Property*
(Joback etc.) are ignored. The compound is found by its InChIKey (from the SMILES) and, failing
that, by name. When only a melting point is known, ΔH_fus can be estimated with Walden's rule
(ΔS_fus ≈ 56.5 J/mol/K); such values are flagged as estimates.
"""
from __future__ import annotations

import html
import re
import statistics
import urllib.parse
import urllib.request
from dataclasses import dataclass

BASE = "https://www.chemeo.com"
WALDEN_DS = 56.5  # J/(mol K)
_UA = {"User-Agent": "osnades-simulator/1.0 (fusion data lookup)"}


@dataclass
class FusionRecord:
    dh_fus: float | None  # J/mol
    t_fus: float | None  # K
    source: str  # reference title (or "NIST")
    compound: str  # Cheméo compound name the value belongs to
    url: str
    doi: str | None = None
    exact: bool = True  # compound matched the full InChIKey (stereo included)
    estimated: bool = False  # ΔH_fus from Walden's rule

    @property
    def label(self) -> str:
        parts = []
        if self.dh_fus:
            parts.append(f"ΔH_fus {self.dh_fus / 1000:.2f} kJ/mol" + (" (Walden estimate)" if self.estimated else ""))
        if self.t_fus:
            parts.append(f"T_m {self.t_fus:.2f} K")
        ref = f"doi:{self.doi}" if self.doi else self.source
        return f"{' · '.join(parts)} — {ref}" + ("" if self.exact else f" [entry: {self.compound}]")

    @property
    def note(self) -> str:
        """Text stored with the species so the origin of the numbers stays visible."""
        ref = f"{self.source} (doi:{self.doi})" if self.doi else self.source
        s = f"Fusion data via Cheméo: {ref}"
        if self.estimated:
            s += f"; ΔH_fus estimated from T_m with Walden's rule (ΔS_fus = {WALDEN_DS} J/mol/K)"
        if not self.exact:
            s += f"; value belongs to the Cheméo entry '{self.compound}' (stereo not matched) – verify"
        return s


def _ssl_context():
    import ssl

    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:  # pragma: no cover
        return ssl.create_default_context()


def _get(url: str, timeout: float) -> str:
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout, context=_ssl_context()) as r:  # noqa: S310 (fixed https host)
        return r.read().decode("utf-8", "replace")


def _search(query: str, timeout: float) -> list[str]:
    page = _get(f"{BASE}/search?q={urllib.parse.quote(query)}", timeout)
    out: list[str] = []
    for m in re.finditer(r'href="(/cid/[0-9-]+/[^"./]+)"', page):
        if m.group(1) not in out:
            out.append(m.group(1))
    return out


def _text(s: str) -> str:
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", s))).strip()


def parse_page(page: str, url: str = "", exact: bool = True) -> list[FusionRecord]:
    """Experimental ΔH_fus / T_fus rows of one Cheméo compound page, paired by reference."""
    m = re.search(r"<h1[^>]*>(.*?)</h1>", page, re.S)
    compound = _text(m.group(1)) if m else url.rsplit("/", 1)[-1]
    compound = re.sub(r"^Chemical Properties of\s+", "", compound)
    rows: dict[str, dict] = {}
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", page, re.S):
        tds = re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)
        if len(tds) < 4:
            continue
        prop = _text(tds[0]).replace(" ", "")
        key = {"ΔfusH°": "dh", "Tfus": "tm"}.get(prop)
        if not key or "Calculated Property" in tds[3]:
            continue
        val = _text(tds[1])
        if not re.fullmatch(r"-?\d+(\.\d+)?", val.split("±")[0].strip()):
            continue  # ranges like "[27.42; 34.30]" are summaries of the rows below
        v = float(val.split("±")[0])
        unit = _text(tds[2])
        if key == "dh":
            v *= 1000 if unit == "kJ/mol" else 1  # → J/mol
        src_m = re.search(r'class="src-name"[^>]*href="#ref-([^"]*)"[^>]*title="([^"]*)"', tds[3])
        ref_id, title = (src_m.group(1), html.unescape(src_m.group(2))) if src_m else ("", _text(tds[3]))
        r = rows.setdefault(ref_id or title, {"title": title, "doi": None})
        if ref_id.startswith("doi_"):
            r["doi"] = ref_id[4:].replace("_", "/", 1)
        r.setdefault(key, v)
    return [FusionRecord(dh_fus=r.get("dh"), t_fus=r.get("tm"), source=r["title"], compound=compound,
                         url=url, doi=r["doi"], exact=exact)
            for r in rows.values() if r.get("dh") or r.get("tm")]


def _inchikey(smiles: str | None = None, molblock: str | None = None) -> str | None:
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles) if smiles else Chem.MolFromMolBlock(molblock or "")
    return Chem.MolToInchiKey(mol) if mol is not None else None


def lookup(smiles: str | None = None, name: str | None = None, timeout: float = 15.0,
           max_pages: int = 12) -> list[FusionRecord]:
    """All experimental fusion records for a compound (empty list if none). Raises on network errors.

    Cheméo's search also returns similar compounds, so each hit's structure is checked: same InChIKey
    = exact match; same skeleton (first InChIKey block) = stereoisomer / unspecified-stereo entry."""
    key = _inchikey(smiles) if smiles else None
    hits = _search(key, timeout) if key else []
    if not hits and name:
        hits = _search(name, timeout)
    recs: list[FusionRecord] = []
    for path in hits[:max_pages]:
        url = BASE + path
        exact = None
        if key:
            try:
                k = _inchikey(molblock=_get(url + ".2dmol", timeout))
            except Exception:  # noqa: BLE001
                k = None
            if k == key:
                exact = True
            elif k and k.split("-")[0] == key.split("-")[0]:
                exact = False
            else:
                continue
        recs += parse_page(_get(url, timeout), url, exact if exact is not None else False)
    return recs


def _clusters(recs: list[FusionRecord], tol: float = 3.0) -> list[list[FusionRecord]]:
    out: list[list[FusionRecord]] = []
    for r in sorted((r for r in recs if r.t_fus), key=lambda r: r.t_fus):
        if out and r.t_fus - out[-1][-1].t_fus <= tol:
            out[-1].append(r)
        else:
            out.append([r])
    return out


def best(records: list[FusionRecord], estimate: bool = True) -> FusionRecord | None:
    """Pick one record. Melting points are grouped (±3 K) and the group backed by the most independent
    sources wins (data collections sometimes hold values of a glass transition or another polymorph).
    Within it, a reference reporting both ΔH_fus and T_m is preferred (DSC measurement); if there is
    none, any paired reference; otherwise the group's median T_m with a Walden's-rule estimate."""
    for exact in (True, False):
        pool = [r for r in records if r.exact == exact]
        groups = _clusters(pool)
        if not groups:
            continue
        def score(g):
            srcs = {r.source for r in g}
            return len(srcs) + 0.5 * any(s == "NIST" for s in srcs) + 0.25 * any(r.dh_fus for r in g)
        g = max(groups, key=score)
        med = statistics.median(r.t_fus for r in g)
        paired = [r for r in g if r.dh_fus] or [r for r in pool if r.dh_fus and r.t_fus]
        if paired:
            return min(paired, key=lambda r: (r.doi is None, abs(r.t_fus - med)))
        if estimate:
            r = min(g, key=lambda r: abs(r.t_fus - med))
            return FusionRecord(dh_fus=WALDEN_DS * r.t_fus, t_fus=r.t_fus, source=r.source, compound=r.compound,
                                url=r.url, doi=r.doi, exact=r.exact, estimated=True)
    return None


def find(smiles: str | None = None, name: str | None = None, timeout: float = 15.0,
         estimate: bool = True) -> FusionRecord | None:
    return best(lookup(smiles, name, timeout), estimate=estimate)
