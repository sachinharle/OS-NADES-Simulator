"""OS-NADES Simulator – Streamlit web application.

Run with ``osnades app`` or ``streamlit run osnades/app.py``.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

if __package__ in (None, ""):  # executed by `streamlit run path/to/app.py`
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

from osnades import __version__, orca
from osnades.config import load_settings
from osnades.jobs import JobManager
from osnades.library import Library
from osnades.sigma import descriptors, sigma_profile
from osnades.thermo import Engine, SolventSystem, ideal_rhs, presets, preset_systems

st.set_page_config(page_title="OS-NADES Simulator", page_icon="⚗️", layout="wide")


# --------------------------------------------------------------------------- state
@st.cache_resource
def job_manager() -> JobManager:
    return JobManager(workers=1)


def settings():
    if "settings" not in st.session_state:
        st.session_state.settings = load_settings()
    return st.session_state.settings


def library() -> Library:
    return Library(settings())


@st.cache_resource
def _engine(home: str, par: str) -> Engine:
    return Engine(Library(load_settings()), par)


def engine(par: str = "openCOSMORS24a") -> Engine:
    return _engine(str(settings().home), par)


def mol_view(xyz: str, height: int = 360) -> None:
    import py3Dmol

    v = py3Dmol.view(width=480, height=height)
    v.addModel(xyz, "xyz")
    v.setStyle({}, {"stick": {"radius": 0.12}, "sphere": {"scale": 0.28}})
    v.setBackgroundColor("0xffffff", 0.0)
    v.zoomTo()
    components.html(v._make_html(), height=height + 10)


def queue_compute(name: str, **kw) -> None:
    lib = library()
    job_manager().submit(f"Compute {name}", lib.compute, name, **kw)


def missing_species(systems: list[SolventSystem]) -> list[str]:
    lib = library()
    return sorted({c for s in systems for c in s.species if not lib.has(c)})


# --------------------------------------------------------------------------- sidebar
st.sidebar.title("⚗️ OS-NADES Simulator")
st.sidebar.caption(f"v{__version__} · ORCA + OPI + openCOSMO-RS")
page = st.sidebar.radio(
    "Navigate",
    ["NADES workflow", "Overview", "Compound library", "Solubility", "Temperature curves", "NADES designer", "Jobs",
     "Settings"],
    label_visibility="collapsed",
)
active = job_manager().active()
if active:
    st.sidebar.info(f"⏳ {len(active)} ORCA job(s) running/queued")
st.sidebar.markdown("---")
st.sidebar.caption(
    "Workflow: SMILES → RDKit 3D → GFN2-xTB opt → ORCA `!COSMORS` (BP86/def2-TZVPD + CPCM) "
    "→ σ-surface → openCOSMO-RS 24a → γ, solubility."
)


# =========================================================================== pages
def page_overview():
    st.title("OS-NADES Simulator")
    st.markdown(
        """
Predict **solid solubilities** of compounds in conventional solvents and **(natural) deep eutectic
solvents** from first principles, following the FACCTs OPI *openCOSMO-RS* tutorial and extending it to
multi-component solvents.

**Start with the NADES workflow page**: it guides you from a solute and a set of candidate components to
the best mole fractions over a temperature range (σ-surfaces and fusion data are computed or looked up on
the way). The other pages are tools for individual tasks:

1. **Compound library** – add molecules/ions (SMILES, XYZ or existing `.orcacosmo` files). ORCA computes
   their COSMO σ-surfaces.
2. **Solubility** – choose a solute, its fusion enthalpy and melting point, and screen solvents/NADES over
   temperature.
3. **NADES designer** – build custom HBA:HBD mixtures (± water), check species activity coefficients and
   binary excess Gibbs energies, and compute eutectic phase diagrams of 2 or 3 components.
"""
    )
    lib = library()
    sp = lib.all()
    systems = preset_systems(lib)
    ready = [s for s in systems if not missing_species([s])]
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Species in library", len(sp))
    c2.metric("Solvent presets ready", f"{len(ready)}/{len(systems)}")
    c3.metric("Solutes with fusion data", sum(1 for s in sp if s.dh_fus and s.t_fus))
    st_ = settings()
    c4.metric("ORCA backend", st_.resolved_mode)

    st.subheader("Quick start")
    miss = missing_species(systems)
    pre = presets()["species"]
    if miss:
        st.write(f"{len(miss)} preset species are not computed yet: " + ", ".join(miss))
        cols = st.columns([1, 1, 2])
        if cols[0].button("Compute all preset species", type="primary"):
            for m in miss:
                p = pre[m]
                queue_compute(m, smiles=p["smiles"], **{k: p[k] for k in ("dh_fus", "t_fus", "notes") if k in p})
            st.success(f"Queued {len(miss)} ORCA jobs – see the Jobs page.")
        cols[1].caption("Each species takes ~1–3 min under Docker emulation (more for large molecules).")
    else:
        st.success("All preset solvent species are available.")


# --------------------------------------------------------------------------- library
def fusion_caption(sp) -> None:
    """Show where a species' ΔH_fus / T_m come from (and warn about estimates)."""
    src = sp.fusion_source or (sp.notes if sp.dh_fus and "usion" in (sp.notes or "") else "")
    if not (sp.dh_fus or sp.t_fus):
        return
    if "Walden" in src:
        st.warning(f"ΔH_fus of {sp.name} is an estimate (Walden's rule) – replace it with a measured value if you can. "
                   f"{src}")
    elif "verify" in src:
        st.caption(f"⚠️ Fusion data: {src}")
    elif src:
        st.caption(f"Fusion data: {src}")


def fusion_editor(lib: Library, sp) -> None:
    """Fusion-data form of one species with an online lookup that lists every published value."""
    name = sp.name
    with st.form(f"fus_{name}"):
        st.markdown("**Solid–liquid data** (needed when used as solute)")
        a, b, c = st.columns(3)
        dh = a.number_input("ΔH_fus (kJ/mol)", value=float(sp.dh_fus / 1000) if sp.dh_fus else 0.0, min_value=0.0)
        tm = b.number_input("T_m (K)", value=float(sp.t_fus or 0.0), min_value=0.0)
        dcp = c.number_input("ΔCp_fus (J/mol/K)", value=float(sp.dcp_fus or 0.0))
        notes = st.text_input("Notes", value=sp.notes)
        if st.form_submit_button("Save"):
            changed = (dh * 1000 or None, tm or None) != (sp.dh_fus, sp.t_fus)
            lib.update(name, dh_fus=dh * 1000 or None, t_fus=tm or None, dcp_fus=dcp or None, notes=notes,
                       fusion_source="entered by user" if changed else sp.fusion_source)
            st.success("Saved")
            st.rerun()
    fusion_caption(sp)
    if sp.charge != 0:
        return
    key = f"fus_cands_{name}"
    if st.button("Look up published values", key=f"fus_lookup_{name}",
                 help="Searches Cheméo (NIST WebBook + literature DSC data) by InChIKey."):
        from osnades import fusion

        with st.spinner(f"Searching fusion data for {name}…"):
            try:
                st.session_state[key] = fusion.lookup(sp.smiles, name)
            except Exception as exc:  # noqa: BLE001
                st.error(f"Lookup failed: {exc}")
    cands = st.session_state.get(key)
    if cands is None:
        return
    from osnades import fusion

    pick = fusion.choose(cands, name)
    if pick is None:
        st.info(f"No published melting data found for {name}. Enter literature or DSC values above.")
        return
    opts = [pick] + [r for r in cands if r is not pick and r.t_fus and r.dh_fus] + \
           [r for r in cands if r is not pick and not (r.t_fus and r.dh_fus)]
    i = st.selectbox(f"Published values for {name} ({len(cands)} found; first = suggested)", range(len(opts)),
                     format_func=lambda i: opts[i].label, key=f"fus_pick_{name}")
    r = opts[i]
    if not r.dh_fus:
        r = fusion.FusionRecord(dh_fus=fusion.WALDEN_DS * r.t_fus, t_fus=r.t_fus, source=r.source, compound=r.compound,
                                url=r.url, doi=r.doi, match=r.match, estimated=True)
        st.caption(f"This reference gives only T_m; ΔH_fus would be estimated as {r.dh_fus / 1000:.1f} kJ/mol (Walden).")
    if st.button("Use these values", key=f"fus_use_{name}", type="primary"):
        lib.update(name, dh_fus=r.dh_fus or sp.dh_fus, t_fus=r.t_fus or sp.t_fus, fusion_source=r.note)
        st.session_state.pop(key, None)
        st.toast(f"Fusion data of {name} updated")
        st.rerun()


def fusion_status(sp) -> str:
    src = sp.fusion_source or ""
    if not (sp.dh_fus and sp.t_fus):
        return "missing"
    if "Walden" in src or "estimate" in src:
        return "estimate"
    if src.startswith("entered by user"):
        return "entered by you"
    return "measured"


def fusion_data_panel(lib: Library) -> None:
    """Fusion data of every neutral species: status overview, bulk lookup, review and apply."""
    from osnades import fusion

    neutral = [s for s in lib.all() if s.charge == 0]
    if not neutral:
        st.info("No neutral species in the library yet.")
        return
    status = {s.name: fusion_status(s) for s in neutral}
    counts = pd.Series(status).value_counts()
    cols = st.columns(4)
    for c, k in zip(cols, ["measured", "entered by you", "estimate", "missing"]):
        c.metric(k.capitalize(), int(counts.get(k, 0)))
    st.dataframe(pd.DataFrame([{
        "species": s.name, "status": status[s.name],
        "ΔH_fus (kJ/mol)": round(s.dh_fus / 1000, 2) if s.dh_fus else None, "T_m (K)": s.t_fus,
        "source": s.fusion_source or (s.notes if s.dh_fus else "")} for s in neutral]),
        hide_index=True, use_container_width=True,
        column_config={"source": st.column_config.TextColumn(width="large")})

    st.markdown("**Search published data** – Cheméo and the NIST WebBook, matched by structure (InChIKey). "
                "Ions are skipped; salts are handled as their ions on the eutectic page.")
    scope = st.radio("Species to search", ["Missing", "Missing + estimates", "All (re-check everything)"],
                     horizontal=True, key="fus_scope",
                     help="Values you entered yourself are only replaced if you tick them below.")
    want = {"Missing": {"missing"}, "Missing + estimates": {"missing", "estimate"}}.get(
        scope, {"missing", "estimate", "measured", "entered by you"})
    targets = [s.name for s in neutral if status[s.name] in want]
    if st.button(f"Search {len(targets)} species", type="primary", disabled=not targets, key="fus_scan_btn"):
        bar, found = st.progress(0.0), {}
        for i, n in enumerate(targets):
            bar.progress(i / len(targets), text=f"{n} ({i + 1}/{len(targets)})")
            sp = lib.get(n)
            try:
                found[n] = (fusion.choose(fusion.lookup(sp.smiles, n), n), None)
            except Exception as exc:  # noqa: BLE001
                found[n] = (None, str(exc))
        bar.empty()
        st.session_state["fus_scan"] = found
    if not targets:
        st.success("Nothing to search in this category.")

    found = st.session_state.get("fus_scan")
    if not found:
        return
    rows = []
    for n, (rec, err) in found.items():
        if not lib.has(n):
            continue
        sp = lib.get(n)
        cur = fusion_status(sp)
        if rec is None:
            rows.append({"apply": False, "species": n, "current": cur,
                         "ΔH now": sp.dh_fus and round(sp.dh_fus / 1000, 2), "T_m now": sp.t_fus,
                         "ΔH new": None, "T_m new": None, "kind": "error" if err else "nothing found",
                         "source": err or ""})
            continue
        same = sp.dh_fus and sp.t_fus and abs(sp.dh_fus - rec.dh_fus) < 1 and abs(sp.t_fus - rec.t_fus) < 0.01
        kind = "estimate" if rec.estimated else ("check isomer" if rec.needs_check else "measured")
        default = not same and (cur == "missing" or (cur == "estimate" and not rec.estimated))
        rows.append({"apply": bool(default), "species": n, "current": "unchanged" if same else cur,
                     "ΔH now": sp.dh_fus and round(sp.dh_fus / 1000, 2), "T_m now": sp.t_fus,
                     "ΔH new": round(rec.dh_fus / 1000, 2), "T_m new": rec.t_fus, "kind": kind, "source": rec.note})
    if not rows:
        return
    st.markdown("**Results** – ticked rows will be saved. Defaults: fill missing values and replace estimates "
                "by measured ones; other rows are left for you to decide.")
    edited = st.data_editor(pd.DataFrame(rows), hide_index=True, use_container_width=True, key="fus_scan_editor",
                            disabled=[c for c in rows[0] if c != "apply"],
                            column_config={"apply": st.column_config.CheckboxColumn("save"),
                                           "source": st.column_config.TextColumn(width="large")})
    picked = [r for r in edited.to_dict("records") if r["apply"] and found[r["species"]][0] is not None]
    a, b = st.columns([1, 3])
    if a.button(f"Save {len(picked)} selected", type="primary", disabled=not picked, key="fus_apply"):
        for r in picked:
            rec = found[r["species"]][0]
            lib.update(r["species"], dh_fus=rec.dh_fus, t_fus=rec.t_fus, fusion_source=rec.note)
        st.session_state.pop("fus_scan", None)
        st.toast(f"Saved fusion data for {len(picked)} species")
        st.rerun()
    if b.button("Discard results", key="fus_discard"):
        st.session_state.pop("fus_scan", None)
        st.rerun()


def page_library():
    st.title("Compound library")
    lib = library()
    all_sp = lib.all()

    tab_list, tab_fus, tab_add, tab_import, tab_sigma = st.tabs(
        ["Species", "Fusion data", "Add with ORCA", "Import .orcacosmo", "σ-profiles"])

    with tab_fus:
        fusion_data_panel(lib)

    with tab_list:
        if not all_sp:
            st.info("Library is empty. Add species in the next tab or from the Overview page.")
        else:
            df = pd.DataFrame(
                [{"name": s.name, "SMILES": s.smiles, "formula": s.formula, "charge": s.charge,
                  "M (g/mol)": round(s.molar_mass, 2), "conformers": len(s.conformers),
                  "ΔH_fus (kJ/mol)": s.dh_fus / 1000 if s.dh_fus else None, "T_m (K)": s.t_fus,
                  "source": s.source, "optimiser": s.optimizer, "created": s.created} for s in all_sp]
            )
            st.dataframe(df, hide_index=True, use_container_width=True)
            my_solvents_panel(lib)
            n_miss = sum(1 for s in all_sp if s.charge == 0 and not (s.dh_fus and s.t_fus))
            if n_miss:
                st.info(f"{n_miss} neutral species have no fusion data – use the **Fusion data** tab to look it up.")
            name = st.selectbox("Inspect species", [s.name for s in all_sp])
            sp = lib.get(name)
            c1, c2 = st.columns([1, 1])
            with c1:
                try:
                    mol_view(lib.xyz(name))
                except Exception as exc:  # noqa: BLE001
                    st.warning(f"Cannot show structure: {exc}")
            with c2:
                d = descriptors(lib.paths(name)[0])
                st.markdown(f"**{name}** · {sp.formula} · charge {sp.charge} · {d['method']}")
                st.write(pd.DataFrame({
                    "property": ["Area (Å²)", "Volume (Å³)", "E_dielectric (kJ/mol)", "σ-moments M0..M2",
                                 "HB donor moments", "HB acceptor moments"],
                    "value": [f"{d['area_A2']:.2f}", f"{d['volume_A3']:.2f}", f"{d['energy_dielectric_kJmol']:.3f}",
                              ", ".join(f"{v:.4g}" for v in d["sigma_moments"][:3]),
                              ", ".join(f"{v:.4g}" for v in d["hb_donor_moments"]),
                              ", ".join(f"{v:.4g}" for v in d["hb_acceptor_moments"])],
                }).set_index("property"))
                fusion_editor(lib, sp)
                preset_pure = next((k for k, v in presets()["solvent_systems"].items() if v == {name: 1}), None)
                if preset_pure:
                    st.caption(f"Already available as the preset solvent “{preset_pure}”.")
                elif sp.charge == 0:
                    on = st.toggle("Use as a solvent", value=lib.is_solvent(name), key=f"solv_{name}",
                                   help="Adds this species as a pure solvent to the Solubility, Temperature curves "
                                        "and NADES pages (group 'My solvents / custom').")
                    if on != lib.is_solvent(name):
                        try:
                            if on:
                                lib.add_solvent(name)
                            else:
                                lib.remove_solvent(name)
                            st.toast(f"{name} {'added to' if on else 'removed from'} solvents")
                        except ValueError as exc:
                            st.error(str(exc))
                else:
                    st.caption("Ions cannot be a solvent on their own – combine them in a mixture "
                               "(NADES designer → Custom → Save as solvent).")
                for p in lib.paths(name):
                    st.download_button(f"Download {p.name}", p.read_bytes(), file_name=f"{name}_{p.name}", key=str(p))
                if st.button(f"Delete {name}", type="secondary"):
                    lib.remove(name)
                    st.rerun()

    with tab_add:
        pre = presets()["species"]
        mode = st.radio("Input", ["Preset", "SMILES", "XYZ"], horizontal=True)
        smiles = xyz = None
        meta = {}
        if mode == "Preset":
            choices = [k for k in pre if not lib.has(k)]
            if not choices:
                st.info("All presets are in the library.")
                return
            name = st.selectbox("Preset", choices)
            smiles = pre[name]["smiles"]
            meta = {k: pre[name][k] for k in ("dh_fus", "t_fus", "notes") if k in pre[name]}
            st.code(smiles)
        elif mode == "SMILES":
            name = st.text_input("Name")
            smiles = st.text_input("SMILES", placeholder="CC(=O)Nc1ccc(O)cc1")
        else:
            name = st.text_input("Name")
            up = st.file_uploader("XYZ file", type=["xyz"])
            xyz = up.read().decode() if up else st.text_area("…or paste XYZ", height=150) or None
        c1, c2, c3, c4 = st.columns(4)
        optimizer = c1.selectbox("Pre-optimisation", ["xtb", "mmff", "none"],
                                 help="GFN2-xTB(ALPB water) via ORCA, MMFF from RDKit, or use geometry as is")
        nconf = c2.number_input("Conformers", 1, 10, 1, help="Lowest-energy RDKit conformers are computed; openCOSMO-RS uses the lowest-energy one after DFT")
        charge = c3.number_input("Charge (XYZ only)", value=0, step=1) if mode == "XYZ" else None
        mult = c4.number_input("Multiplicity (XYZ only)", value=1, min_value=1) if mode == "XYZ" else None
        if smiles:
            try:
                from osnades.structures import smiles_info, smiles_to_conformers

                info = smiles_info(smiles)
                st.caption(f"{info['formula']} · M = {info['molar_mass']:.2f} g/mol · charge {info['charge']}")
                if info["n_fragments"] == 1:
                    with st.expander("Preview 3D (RDKit)"):
                        mol_view(smiles_to_conformers(smiles)[0].to_xyz(), 300)
                else:
                    st.error("Multi-fragment SMILES: add each ion separately (e.g. choline and chloride).")
            except Exception as exc:  # noqa: BLE001
                st.error(str(exc))
        st.markdown("**Fusion data** (needed when the compound is used as a solute)")
        f1, f2, f3 = st.columns([1, 1, 2])
        dh_in = f1.number_input("ΔH_fus (kJ/mol)", value=float(meta["dh_fus"] / 1000) if meta.get("dh_fus") else None,
                                min_value=0.1, placeholder="look up", key=f"add_dh_{name}")
        tm_in = f2.number_input("T_m (K)", value=float(meta["t_fus"]) if meta.get("t_fus") else None, min_value=1.0,
                                placeholder="look up", key=f"add_tm_{name}")
        auto_fus = f3.checkbox("Look up published values for empty fields", value=True, key="add_lookup",
                               help="After the ORCA job, measured ΔH_fus and T_m are fetched from Cheméo (NIST WebBook "
                                    "and literature DSC data, matched by InChIKey). If only T_m is published, ΔH_fus "
                                    "is estimated with Walden's rule and flagged as an estimate.")
        if dh_in:
            meta["dh_fus"] = dh_in * 1000
        if tm_in:
            meta["t_fus"] = tm_in
        if (dh_in or tm_in) and not meta.get("fusion_source"):
            meta["fusion_source"] = "entered by user" if mode != "Preset" else meta.get("notes", "preset")
        overwrite = st.checkbox("Recompute if this species already exists (replaces its σ-surfaces)", key="add_ow")
        as_solvent = st.checkbox("Also use it as a solvent", key="add_solv",
                                 help="Listed with the solvent systems once the ORCA job has finished.")
        if st.button("Queue ORCA calculation", type="primary", disabled=not name or not (smiles or xyz)):
            clash = lib.clash(name)
            if clash:
                st.error(f"“{name}” would overwrite the files of existing species “{clash}” (names differ only "
                         f"in case/spaces). Use the name “{clash}” with the recompute option instead.")
            elif lib.has(name) and not overwrite:
                st.error(f"{name} already exists – tick the recompute option to replace it.")
            else:
                from osnades.structures import smiles_info as _info

                q = charge if charge is not None else (_info(smiles)["charge"] if smiles else 0)
                if as_solvent and q:
                    st.warning("Charged species cannot be a solvent on their own – not added as solvent.")
                    as_solvent = False
                queue_compute(name, smiles=smiles, xyz=xyz, charge=charge, multiplicity=mult,
                              n_conformers=int(nconf), optimizer=optimizer, overwrite=overwrite,
                              lookup_fusion=auto_fus, **meta)
                if as_solvent:
                    try:
                        lib.add_solvent(name, overwrite=True)
                    except ValueError as exc:
                        st.warning(f"Not added as solvent: {exc}")
                st.success(f"Queued {name}" + (" (also registered as a solvent)" if as_solvent else "")
                           + ". Track progress on the Jobs page.")

    with tab_import:
        st.write("Import σ-surfaces computed elsewhere (e.g. the tutorial's `*.solute.orcacosmo`). "
                 "Several files = several conformers of the same species.")
        name = st.text_input("Species name", key="imp_name")
        smi = st.text_input("SMILES (optional, for molar mass)", key="imp_smi")
        files = st.file_uploader("orcacosmo file(s)", accept_multiple_files=True)
        imp_fus = st.checkbox("Look up published fusion data (needs the SMILES or a well-known name)", value=True,
                              key="imp_lookup")
        if st.button("Import", disabled=not (name and files)):
            with tempfile.TemporaryDirectory() as td:
                paths = []
                for f in files:
                    p = Path(td) / f.name
                    p.write_bytes(f.read())
                    paths.append(p)
                try:
                    sp = lib.import_orcacosmo(name, paths, smiles=smi or None)
                    st.success(f"Imported {sp.name} ({len(sp.conformers)} conformer(s), charge {sp.charge})")
                    if imp_fus and sp.charge == 0:
                        with st.spinner("Looking up fusion data…"):
                            msgs = []
                            sp = lib.lookup_fusion(sp.name, log=msgs.append)
                        st.info(" ".join(msgs) or "No fusion data found.")
                except Exception as exc:  # noqa: BLE001
                    st.error(str(exc))

    with tab_sigma:
        names = st.multiselect("Species", [s.name for s in all_sp], default=[s.name for s in all_sp][:3])
        fig = go.Figure()
        for n in names:
            sig, area = sigma_profile(lib.paths(n)[0])
            fig.add_trace(go.Scatter(x=sig, y=area, mode="lines", name=n))
        for xv in (-0.0085, 0.0085):
            fig.add_vline(x=xv, line_dash="dot", line_color="gray")
        fig.update_layout(xaxis_title="σ (e/Å²)", yaxis_title="area p(σ) (Å²)", xaxis_range=[-0.03, 0.03],
                          height=450, legend_title="", margin=dict(t=30))
        st.plotly_chart(fig, use_container_width=True)
        st.caption("Dotted lines: ±0.0085 e/Å² hydrogen-bond thresholds (σ<0 donor, σ>0 acceptor regions).")


# --------------------------------------------------------------------------- solubility
def system_builder(key: str, lib: Library) -> SolventSystem | None:
    names = lib.names()
    st.markdown("**Custom mixture** (molar ratio)")
    rows = []
    for i in range(4):
        a, b = st.columns([3, 1])
        s = a.selectbox(f"Species {i + 1}", ["—"] + names, key=f"{key}_s{i}")
        r = b.number_input("ratio", min_value=0.0, value=1.0, step=0.5, key=f"{key}_r{i}")
        if s != "—" and r > 0:
            rows.append((s, r))
    label = st.text_input("Label", value=":".join(s for s, _ in rows), key=f"{key}_lbl")
    if not rows:
        return None
    system = SolventSystem(label or "custom", {s: r for s, r in rows})
    if st.button("Save as solvent", key=f"{key}_save", help="Keep this mixture in the solvent lists of every page"):
        try:
            lib.add_solvent(system.name, dict(system.components))
            st.success(f"Saved “{system.name}” to your solvents.")
        except ValueError as exc:
            st.error(str(exc))
    return system


def my_solvents_panel(lib: Library) -> None:
    user = lib.user_solvents()
    with st.expander(f"My solvents ({len(user)})"):
        if not user:
            st.caption("None yet. Toggle “Use as a solvent” on a species below, or save a custom mixture "
                       "with “Save as solvent”.")
            return
        for k, v in user.items():
            a, b = st.columns([5, 1])
            missing = [c for c in v if not lib.has(c)]
            a.markdown(f"**{k}** · " + " : ".join(f"{c} {r:g}" for c, r in v.items())
                       + (f" · ⚠️ not computed: {', '.join(missing)}" if missing else ""))
            if b.button("Remove", key=f"rm_solv_{k}"):
                lib.remove_solvent(k)
                st.rerun()


def page_solubility():
    st.title("Solid solubility")
    lib = library()
    all_sp = lib.all()
    if not all_sp:
        st.info("Add compounds first.")
        return
    c1, c2 = st.columns([1, 2])
    with c1:
        solute = st.selectbox("Solute", [s.name for s in all_sp],
                              index=next((i for i, s in enumerate(all_sp) if s.dh_fus), 0))
        sp = lib.get(solute)
        dh = st.number_input("ΔH_fus (kJ/mol)", value=float(sp.dh_fus / 1000) if sp.dh_fus else 25.0, min_value=0.1,
                             key=f"sol_dh_{solute}")
        tm = st.number_input("T_m (K)", value=float(sp.t_fus or 400.0), min_value=1.0, key=f"sol_tm_{solute}")
        dcp = st.number_input("ΔCp_fus (J/mol/K)", value=float(sp.dcp_fus or 0.0), key=f"sol_dcp_{solute}")
        if sp.dh_fus and sp.t_fus:
            fusion_caption(sp)
        else:
            st.caption(f"No fusion data stored for {solute}; the values above are placeholders – use "
                       "“Look up published values” in the Compound library or enter literature values.")
        method = st.radio("Method", ["iterative", "infinite_dilution"], horizontal=True,
                          help="iterative: γ at the saturated composition; infinite_dilution: γ∞ (tutorial's non-iterative)")
        water = st.slider("Add water to every solvent (wt%)", 0.0, 50.0, 0.0, 1.0)
        par = st.selectbox("Parameterisation", ["openCOSMORS24a"])
    with c2:
        presets_all = preset_systems(lib)
        avail = [s for s in presets_all if not missing_species([s]) and solute not in s.components]
        unavailable = [s for s in presets_all if missing_species([s])]
        chosen = st.multiselect("Solvent systems", [s.name for s in avail], default=[s.name for s in avail][:12])
        with st.expander("Add a custom mixture"):
            custom = system_builder("sol", lib)
        if unavailable:
            with st.expander(f"{len(unavailable)} preset systems need species that are not computed yet"):
                miss = missing_species(unavailable)
                st.write(", ".join(miss))
                if st.button("Queue missing species"):
                    pre = presets()["species"]
                    for m in miss:
                        if m in pre:
                            queue_compute(m, smiles=pre[m]["smiles"])
                    st.success("Queued – see Jobs.")
        tmode = st.radio("Temperature", ["single", "range"], horizontal=True)
        if tmode == "single":
            temps = [st.number_input("T (K)", value=298.15)]
        else:
            a, b, n = st.columns(3)
            t0 = a.number_input("from (K)", value=283.15)
            t1 = b.number_input("to (K)", value=333.15)
            nt = n.number_input("points", 2, 50, 6)
            temps = list(np.linspace(t0, t1, int(nt)))

    systems = [s for s in presets_all if s.name in chosen]
    if custom is not None:
        systems.append(custom)
    if not systems:
        st.warning("Select at least one solvent system.")
        return
    if water:
        systems = [s.with_water(water, lib) for s in systems]

    if st.button("Calculate", type="primary"):
        eng = engine(par)
        prog = st.progress(0.0)
        frames = []
        for i, T in enumerate(temps):
            frames.append(eng.screen(solute, systems, float(T), method=method, dh_fus=dh * 1000, t_fus=tm, dcp_fus=dcp))
            prog.progress((i + 1) / len(temps))
        st.session_state.sol_result = pd.concat(frames, ignore_index=True)
        st.session_state.sol_meta = (solute, method)

    df = st.session_state.get("sol_result")
    if df is not None and len(df):
        solute_, method_ = st.session_state.sol_meta
        st.subheader(f"{solute_} · {method_}")
        if "x" not in df:
            st.error(df.to_string())
            return
        show = df[["solvent", "T", "x", "w", "g_per_kg", "ln_gamma", "ln_gamma_inf", "x_ideal", "note"]].copy()
        if df["T"].nunique() == 1:
            fig = px.bar(df.sort_values("x"), x="x", y="solvent", orientation="h", log_x=True,
                         labels={"x": "mole fraction solubility x", "solvent": ""}, height=40 * len(df) + 120)
            fig.add_vline(x=df["x_ideal"].iloc[0], line_dash="dash", annotation_text="ideal")
            st.plotly_chart(fig, use_container_width=True)
        else:
            fig = px.line(df, x="T", y="x", color="solvent", log_y=True, markers=True,
                          labels={"T": "T (K)", "x": "mole fraction solubility x"})
            st.plotly_chart(fig, use_container_width=True)
        st.dataframe(show.style.format({"x": "{:.4e}", "w": "{:.4e}", "g_per_kg": "{:.4g}", "ln_gamma": "{:.3f}",
                                        "ln_gamma_inf": "{:.3f}", "x_ideal": "{:.4e}", "T": "{:.2f}"}),
                     hide_index=True, use_container_width=True)
        st.download_button("Download CSV", df.to_csv(index=False).encode(), file_name=f"solubility_{solute_}.csv")
        st.caption("x counts every solvent species (ions individually). g_per_kg = g solute per kg solvent. "
                   "Ideal-solid SLE model; polymorphs/solvates and solute dissociation are not considered.")


# --------------------------------------------------------------------------- temperature curves
def page_curves():
    from osnades.plots import compute_curves, group_systems, matplotlib_png, plotly_panels

    st.title("Solubility vs temperature")
    st.caption("Curves for every selected solvent system, grouped by solvent class. Each panel has its own "
               "y-scale: linear when values are close together, log otherwise.")
    lib = library()
    solutes = [s for s in lib.all() if s.charge == 0]
    if not solutes:
        st.info("Add a neutral compound to the library first.")
        return
    c1, c2, c3 = st.columns([1.2, 1, 1])
    with c1:
        names = [s.name for s in solutes]
        solute = st.selectbox("Solute", names, index=next((i for i, s in enumerate(solutes) if s.dh_fus and s.t_fus), 0),
                              key="cv_solute", format_func=lambda n: n if lib.get(n).dh_fus and lib.get(n).t_fus
                              else f"{n}  (no fusion data)")
        sp = lib.get(solute)
        a, b, c = st.columns(3)
        dh = a.number_input("ΔH_fus (kJ/mol)", value=float(sp.dh_fus / 1000) if sp.dh_fus else None, min_value=0.1,
                            key=f"cv_dh_{solute}", placeholder="required")
        tm = b.number_input("T_m (K)", value=float(sp.t_fus) if sp.t_fus else None, min_value=1.0,
                            key=f"cv_tm_{solute}", placeholder="required")
        dcp = c.number_input("ΔCp_fus (J/mol/K)", value=float(sp.dcp_fus or 0.0), key=f"cv_dcp_{solute}")
        fusion_caption(sp)
        if dh is None or tm is None:
            st.warning(f"Enter the melting enthalpy and melting point of {solute} (literature or DSC) to "
                       "compute its solubility curves.")
            return
        if (sp.dh_fus, sp.t_fus, sp.dcp_fus or 0.0) != (dh * 1000, tm, dcp) and st.button(
                "Save fusion data to library", key="cv_save"):
            lib.update(solute, dh_fus=dh * 1000, t_fus=tm, dcp_fus=dcp or None, fusion_source="entered by user")
            st.toast(f"Saved fusion data for {solute}")
            st.rerun()
    with c2:
        t0, t1 = st.slider("Temperature range (°C)", -20, 150, (0, 80), key="cv_range")
        npts = st.slider("Points per curve", 5, 41, 17, key="cv_n")
        if t1 + 273.15 >= tm:
            st.warning("Range reaches T_m; above it the solute is liquid and SLE does not apply.")
    with c3:
        method = st.radio("Method", ["iterative", "infinite_dilution"], key="cv_method")
        water = st.slider("Add water (wt%)", 0.0, 50.0, 0.0, 1.0, key="cv_water")

    ready = [s for s in preset_systems(lib) if not missing_species([s]) and solute not in s.components]
    groups = group_systems([s.name for s in ready])
    with st.expander(f"Solvent systems ({len(ready)} ready)", expanded=False):
        chosen = []
        cols = st.columns(2)
        for i, (g, names) in enumerate(groups.items()):
            with cols[i % 2]:
                chosen += st.multiselect(g, names, default=names, key=f"cv_g_{i}")
        custom = system_builder("cv", lib)
    systems = [s for s in ready if s.name in chosen]
    if custom is not None:
        systems.append(custom)
    if water:
        systems = [s.with_water(water, lib) for s in systems]
    if not systems:
        st.warning("Select at least one solvent system.")
        return

    if st.button(f"Compute {len(systems)} curves", type="primary", key="cv_go"):
        temps = np.linspace(t0, t1, int(npts)) + 273.15
        bar = st.progress(0.0)
        df = compute_curves(engine(), solute, systems, temps, method=method,
                            progress=lambda f, n: bar.progress(f, text=n), dh_fus=dh * 1000, t_fus=tm, dcp_fus=dcp)
        bar.empty()
        st.session_state.cv = dict(df=df, solute=solute, dh=dh * 1000, tm=tm, dcp=dcp, method=method)

    res = st.session_state.get("cv")
    if not res:
        return
    df = res["df"]
    st.subheader(f"{res['solute']} · {res['method']}")
    if df["note"].fillna("").str.len().sum():
        st.warning("Some points carry solver notes – see the table.")
    figs = plotly_panels(df, res["dh"], res["tm"], res["dcp"])
    cols = st.columns(2)
    for i, (g, fig, log, ideal) in enumerate(figs):
        with cols[i % 2]:
            st.plotly_chart(fig, use_container_width=True, key=f"cv_fig_{i}")
            if not ideal:
                st.caption("Ideal solubility lies outside this panel's range.")

    lo, hi = df["T"].min(), df["T"].max()
    gain = (df[df["T"] == hi].set_index("solvent")["x"] / df[df["T"] == lo].set_index("solvent")["x"])
    summary = pd.DataFrame({
        f"x at {lo - 273.15:.0f} °C": df[df["T"] == lo].set_index("solvent")["x"],
        f"x at {hi - 273.15:.0f} °C": df[df["T"] == hi].set_index("solvent")["x"],
        "gain (×)": gain,
    }).sort_values(f"x at {hi - 273.15:.0f} °C", ascending=False)
    st.markdown("**Temperature gain** – a high gain marks good candidates for cooling crystallisation.")
    st.dataframe(summary.style.format("{:.4g}"), use_container_width=True)
    d1, d2 = st.columns(2)
    d1.download_button("Download data (CSV)", df.to_csv(index=False).encode(),
                       file_name=f"{res['solute']}_solubility_vs_T.csv")
    d2.download_button("Download figure (PNG)", matplotlib_png(df, res["solute"], res["dh"], res["tm"], res["dcp"]),
                       file_name=f"{res['solute']}_solubility_vs_T.png", mime="image/png")


# --------------------------------------------------------------------------- NADES designer
def nades_solute_panel(sysm: SolventSystem, T: float, lib: Library) -> None:
    """Ideal vs openCOSMO-RS solubility of one solute in the designed mixture."""
    st.subheader("Solute in this mixture")
    cands = [s for s in lib.all() if s.charge == 0 and s.name not in sysm.components]
    if not cands:
        st.info("No neutral species outside this mixture are in the library.")
        return
    c1, c2, c3, c4 = st.columns(4)
    names = [s.name for s in cands]
    solute = c1.selectbox("Solute", names, index=next((i for i, s in enumerate(cands) if s.dh_fus), 0), key="nd_sol")
    sp = lib.get(solute)
    dh = c2.number_input("ΔH_fus (kJ/mol)", value=float(sp.dh_fus / 1000) if sp.dh_fus else 25.0,
                         min_value=0.1, key=f"nd_dh_{solute}")
    tm = c3.number_input("T_m (K)", value=float(sp.t_fus or 400.0), min_value=1.0, key=f"nd_tm_{solute}")
    dcp = c4.number_input("ΔCp_fus (J/mol/K)", value=float(sp.dcp_fus or 0.0), key=f"nd_dcp_{solute}")
    if not sp.dh_fus:
        st.caption(f"No fusion data stored for {solute}; the defaults above are placeholders – enter literature values.")
    else:
        fusion_caption(sp)

    x_id = float(np.exp(ideal_rhs(T, dh * 1000, tm, dcp)))
    mw = sp.molar_mass or 0.0
    mw_solv = sum(f * (lib.get(s).molar_mass or 0.0) for s, f in zip(sysm.species, sysm.fractions))
    g_kg_id = 1000 * x_id * mw / ((1 - x_id) * mw_solv) if mw and mw_solv and x_id < 1 else float("nan")
    fus = dict(dh_fus=dh * 1000, t_fus=tm, dcp_fus=dcp)
    try:
        res = engine().solubility(solute, sysm, T, method="iterative", **fus)
        inf = engine().solubility(solute, sysm, T, method="infinite_dilution", **fus)
    except Exception as e:  # noqa: BLE001
        st.error(f"Solubility calculation failed: {e}")
        return

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("x ideal", f"{x_id:.4g}", help="exp[−ΔH/R(1/T−1/T_m) + ΔCp/R(T_m/T−1−ln T_m/T)]; "
              "solvent-independent (γ = 1)")
    m2.metric("x predicted", f"{res.x:.4g}", f"{res.x / x_id:.2f}× ideal", delta_color="off")
    m3.metric("ln γ at saturation", f"{res.ln_gamma:.3f}", help="Pure-liquid solute reference. ln γ < 0 → "
              "solubility above ideal")
    m4.metric("g / kg solvent", f"{res.g_per_kg:.4g}", f"ideal {g_kg_id:.4g}", delta_color="off")

    tbl = pd.DataFrame([
        {"case": "ideal (γ = 1)", "x": x_id, "ln_gamma": 0.0, "g_per_kg": g_kg_id, "x / x_ideal": 1.0},
        {"case": "openCOSMO-RS iterative", "x": res.x, "ln_gamma": res.ln_gamma, "g_per_kg": res.g_per_kg,
         "x / x_ideal": res.x / x_id},
        {"case": "openCOSMO-RS infinite dilution", "x": inf.x, "ln_gamma": inf.ln_gamma_inf, "g_per_kg": inf.g_per_kg,
         "x / x_ideal": inf.x / x_id},
    ])
    st.dataframe(tbl.style.format({"x": "{:.4e}", "ln_gamma": "{:.3f}", "g_per_kg": "{:.4g}", "x / x_ideal": "{:.3f}"}),
                 hide_index=True, use_container_width=True)
    if res.note:
        st.warning(res.note)

    with st.expander("Ideal vs predicted over a temperature range"):
        a, b, n = st.columns(3)
        t0 = a.number_input("from (K)", value=283.15, key="nd_t0")
        t1 = b.number_input("to (K)", value=min(353.15, tm - 1), key="nd_t1")
        nt = int(n.number_input("points", 2, 30, 8, key="nd_nt"))
        if st.button("Compute curve", key="nd_curve"):
            temps = np.linspace(t0, t1, nt)
            df = engine().solubility_curve(solute, sysm, temps, method="iterative", **fus)
            st.session_state.nd_curve_df = (solute, sysm.name, df)
        cur = st.session_state.get("nd_curve_df")
        if cur and cur[0] == solute and cur[1] == sysm.name:
            df = cur[2]
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=df["T"] - 273.15, y=df["x_ideal"], name="ideal (γ = 1)",
                                     line=dict(color="#5f5e58", dash="dash")))
            fig.add_trace(go.Scatter(x=df["T"] - 273.15, y=df["x"], name=f"openCOSMO-RS in {sysm.name}",
                                     mode="lines+markers", line=dict(color="#2a78d6", width=2)))
            fig.update_layout(xaxis_title="Temperature (°C)", yaxis_title="mole fraction x", yaxis_type="log",
                              height=380, legend=dict(orientation="h", yanchor="top", y=-0.25, x=0))
            st.plotly_chart(fig, use_container_width=True)
            st.download_button("Download CSV", df.to_csv(index=False).encode(),
                               file_name=f"{solute}_in_{sysm.name}.csv".replace(" ", "_"), key="nd_dl")
    st.caption("x counts every mixture species individually (ions separately). The ideal mole fraction depends only "
               "on the solute's fusion properties and T; x/x_ideal = 1/γ shows how much the NADES helps (>1) or hurts (<1).")


def page_nades():
    st.title("NADES designer")
    lib = library()
    if len(lib.names()) < 2:
        st.info("Add at least two species.")
        return
    t1, t2, t3 = st.tabs(["Mixture activity", "Binary G^E", "Eutectic phase diagram"])
    with t1:
        presets_ok = [s for s in preset_systems(lib) if not missing_species([s]) and len(s.components) > 1]
        src = st.radio("System", ["Preset", "Custom"], horizontal=True)
        if src == "Preset" and presets_ok:
            name = st.selectbox("Preset", [s.name for s in presets_ok])
            sysm = next(s for s in presets_ok if s.name == name)
        else:
            sysm = system_builder("nd", lib)
        water = st.slider("Water (wt%)", 0.0, 60.0, 0.0, 1.0, key="nd_w")
        T = st.number_input("T (K)", value=298.15, key="nd_T")
        if sysm is not None:
            if water:
                sysm = sysm.with_water(water, lib)
            q = sysm.charge_balance(lib)
            if abs(q) > 1e-6:
                st.error(f"Mixture is not electroneutral (net charge per mole = {q:+.3f}). Balance cations and anions.")
            else:
                df = engine().mixture_activity(sysm, T)
                st.dataframe(df.style.format({c: "{:.4f}" for c in df.columns if c != "species"}),
                             hide_index=True, use_container_width=True)
                fig = px.bar(df, x="species", y=[c for c in df.columns if c.startswith("ln_gamma_") and c != "ln_gamma_tot"],
                             barmode="group", labels={"value": "ln γ (COSMO reference)", "variable": ""})
                st.plotly_chart(fig, use_container_width=True)
                st.caption("Reference state: ideal conductor (COSMO). Strongly negative residual ln γ = strong "
                           "stabilisation in the mixture (e.g. H-bonding between HBA and HBD).")
                fig2 = go.Figure()
                tot = None
                for s, f in zip(sysm.species, sysm.fractions):
                    sig, area = sigma_profile(engine().conformer(s))
                    tot = area * f if tot is None else tot + area * f
                fig2.add_trace(go.Scatter(x=sig, y=tot, mode="lines", name="mixture σ-profile", fill="tozeroy"))
                fig2.update_layout(xaxis_title="σ (e/Å²)", yaxis_title="x-weighted p(σ) (Å²)", xaxis_range=[-0.03, 0.03], height=350)
                st.plotly_chart(fig2, use_container_width=True)
                nades_solute_panel(sysm, T, lib)
    with t2:
        neutral = [s.name for s in lib.all() if s.charge == 0]
        a, b, c = st.columns(3)
        s1 = a.selectbox("Component 1", neutral, key="b1")
        s2 = b.selectbox("Component 2", [n for n in neutral if n != s1], key="b2")
        T = c.number_input("T (K)", value=298.15, key="bT")
        if s1 and s2:
            df = engine().binary(s1, s2, T, 41)
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=df.x1, y=df.iloc[:, 1], name=f"ln γ {s1}"))
            fig.add_trace(go.Scatter(x=df.x1, y=df.iloc[:, 2], name=f"ln γ {s2}"))
            fig.add_trace(go.Scatter(x=df.x1, y=df.GE_RT, name="Gᴱ/RT", line=dict(dash="dash")))
            fig.update_layout(xaxis_title=f"x({s1})", yaxis_title="", height=420)
            st.plotly_chart(fig, use_container_width=True)
            st.caption("Negative Gᴱ: favourable mixing (typical for good HBA/HBD eutectic partners).")
    with t3:
        eutectic_panel(lib)


def eutectic_panel(lib: Library) -> None:
    from osnades.eutectic import component_names, fusion_defaults, phase_diagram, plotly_figure

    names = component_names(lib)
    if len(names) < 2:
        st.info("Add at least two components.")
        return
    kind = st.radio("Components", ["2 (binary)", "3 (ternary)"], horizontal=True, key="eu_kind",
                    disabled=len(names) < 3)
    if kind.startswith("3"):
        ternary_panel(lib, names)
        return
    st.markdown("Solid–liquid equilibrium of an **HBA + HBD** pair: liquidus of each component, the eutectic "
                "point and the composition window that is liquid at your operating temperature.")
    cols = st.columns(2)
    comps = []
    defaults = ("choline chloride", "urea") if "choline chloride" in names else (names[0], names[1])
    for i, (col, lab) in enumerate(zip(cols, ("A (HBA)", "B (HBD)"))):
        with col:
            opts = names if i == 0 else [n for n in names if n != comps[0][0]]
            d = defaults[i] if defaults[i] in opts else opts[0]
            name = st.selectbox(f"Component {lab}", opts, index=opts.index(d), key=f"eu_c{i}")
            dh0, tm0, dcp0, src = fusion_defaults(name, lib)
            a, b, c = st.columns(3)
            dh = a.number_input("ΔH_fus (kJ/mol)", value=float(dh0 / 1000) if dh0 else 20.0, min_value=0.1,
                                key=f"eu_dh{i}_{name}")
            tm = b.number_input("T_m (K)", value=float(tm0 or 350.0), min_value=1.0, key=f"eu_tm{i}_{name}")
            dcp = c.number_input("ΔCp (J/mol/K)", value=float(dcp0 or 0.0), key=f"eu_dcp{i}_{name}")
            st.caption(f"Fusion data: {src}" if dh0 and tm0 else "⚠️ No stored fusion data – enter literature values.")
            comps.append((name, dh, tm, dcp))
    a, b, c = st.columns(3)
    T_op = a.number_input("Operating T (K)", value=298.15, key="eu_Top")
    n = int(b.number_input("Composition points", 11, 101, 41, step=10, key="eu_n"))
    c.write("")
    run = c.button("Compute phase diagram", type="primary", key="eu_run")
    if run:
        from osnades.eutectic import resolve_component

        try:
            A, B = (resolve_component(nm, lib, dh * 1000, tm, dcp) for nm, dh, tm, dcp in comps)
            prog = st.progress(0.0, text="Solving liquidus curves…")
            res = phase_diagram(engine(), A, B, n=n, progress=lambda f: prog.progress(f))
            prog.empty()
            st.session_state.eu_res = res
        except Exception as e:  # noqa: BLE001
            st.error(str(e))
    res = st.session_state.get("eu_res")
    if res is None:
        return
    m1, m2, m3, m4 = st.columns(4)
    if np.isfinite(res.T_e):
        m1.metric("Eutectic x_A", f"{res.x_e:.3f}")
        m2.metric("Eutectic T", f"{res.T_e - 273.15:.1f} °C", f"{res.T_e:.1f} K", delta_color="off")
    else:
        m1.metric("Eutectic x_A", f"≈ {res.x_e:.2f}" if np.isfinite(res.x_e) else "—")
        m2.metric("Eutectic T", f"< {res.T_e_upper:.0f} K" if np.isfinite(res.T_e_upper) else "—")
    m3.metric("Ideal eutectic", f"{res.T_e_ideal - 273.15:.1f} °C", f"x = {res.x_e_ideal:.3f}", delta_color="off")
    dep = res.depression if np.isfinite(res.depression) else (res.T_e_ideal - res.T_e_upper)
    m4.metric("Depression", (f"{dep:.1f} K" if np.isfinite(res.depression) else f"> {dep:.0f} K")
              if np.isfinite(dep) else "—", help="Positive = deeper than ideal (a 'deep' eutectic)")
    st.plotly_chart(plotly_figure(res, T_op), use_container_width=True)
    w = res.liquid_window(T_op)
    if w:
        st.success(f"Liquid at {T_op - 273.15:.1f} °C for x_{res.a.name} = {w[0]:.3f} – {w[1]:.3f} "
                   f"(A:B molar ratio {w[0] / (1 - w[0]):.2f} – {w[1] / max(1e-9, 1 - w[1]):.2f}).")
    else:
        st.warning(f"No composition is fully liquid at {T_op - 273.15:.1f} °C.")
    for note in res.notes:
        st.info(note)
    with st.expander("Data"):
        st.dataframe(res.data.style.format(precision=4), hide_index=True, use_container_width=True)
        st.download_button("Download CSV", res.data.to_csv(index=False).encode(),
                           file_name=f"eutectic_{res.a.name}_{res.b.name}.csv".replace(" ", "_"), key="eu_dl")
    st.caption("Ideal-solid SLE with immiscible solids (no co-crystals or solid solutions). Salts are split into ions "
               "and referenced to their pure molten salt; COSMO-RS tends to over-stabilise ions, so ChCl eutectics "
               "come out too deep – compare candidates relative to each other. Fusion data drive the result; "
               "preset values are approximate literature numbers.")


def fusion_inputs(lib: Library, name: str, key: str) -> tuple[float, float, float]:
    """ΔH_fus (kJ/mol), T_m (K), ΔCp inputs pre-filled from the library / presets."""
    from osnades.eutectic import fusion_defaults

    dh0, tm0, dcp0, src = fusion_defaults(name, lib)
    a, b, c = st.columns(3)
    dh = a.number_input("ΔH_fus (kJ/mol)", value=float(dh0 / 1000) if dh0 else 20.0, min_value=0.1,
                        key=f"{key}_dh_{name}")
    tm = b.number_input("T_m (K)", value=float(tm0 or 350.0), min_value=1.0, key=f"{key}_tm_{name}")
    dcp = c.number_input("ΔCp (J/mol/K)", value=float(dcp0 or 0.0), key=f"{key}_dcp_{name}")
    st.caption(f"Fusion data: {src}" if dh0 and tm0 else "⚠️ No stored fusion data – enter literature values.")
    return dh, tm, dcp


def ternary_panel(lib: Library, names: list[str]) -> None:
    from osnades.eutectic import (ratio_label, resolve_component, ternary_cost, ternary_diagram,
                                  ternary_figure, ternary_png)

    st.markdown("Liquidus surface of a **three-component** mixture (e.g. HBA + two HBDs): the melting temperature "
                "of every composition, the ternary eutectic and the region that is liquid at your operating "
                "temperature.")
    pick = ["choline chloride", "urea", "glycerol"] if "choline chloride" in names else names[:3]
    pick = [p for p in pick if p in names] + [n for n in names if n not in pick]
    cols = st.columns(3)
    comps, chosen = [], []
    for i, (col, lab) in enumerate(zip(cols, "ABC")):
        with col:
            opts = [n for n in names if n not in chosen]
            d = next(p for p in pick if p in opts)
            name = st.selectbox(f"Component {lab}", opts, index=opts.index(d), key=f"te_c{i}")
            chosen.append(name)
            comps.append((name, *fusion_inputs(lib, name, f"te{i}")))
    a, b, c = st.columns(3)
    T_op = a.number_input("Operating T (K)", value=298.15, key="te_Top")
    step = b.select_slider("Grid step", [0.2, 0.1, 0.05], value=0.1, key="te_step",
                           help="Composition spacing (mole fraction). The eutectic is then refined locally.")
    T_lo = c.number_input("Lowest T checked (K)", value=float(max(173.15, 0.5 * min(x[2] for x in comps))),
                          min_value=50.0, key=f"te_Tlo_{'_'.join(chosen)}",
                          help="Liquidus temperatures below this are reported as 'below'.")
    from osnades.designer import species_of

    ions = any(len(species_of(n)) > 1 for n in chosen)
    est = ternary_cost(step, ions)
    run = st.button("Compute ternary diagram", type="primary", key="te_run",
                    help=f"Estimated time ≈ {est:.0f} s" + (" (salts are slower)" if ions else ""))
    st.caption(f"Estimated time ≈ {est:.0f} s" + (" (salts are slower)" if ions else ""))
    if run:
        try:
            cs = [resolve_component(nm, lib, dh * 1000, tm, dcp) for nm, dh, tm, dcp in comps]
            prog = st.progress(0.0, text="Solving liquidus surfaces…")
            res = ternary_diagram(engine(), cs, step=step, T_lo=T_lo,
                                  progress=lambda f: prog.progress(min(1.0, f), text="Solving liquidus surfaces…"))
            prog.empty()
            st.session_state.te_res = res
        except Exception as e:  # noqa: BLE001
            st.error(str(e))
    res = st.session_state.get("te_res")
    if res is None:
        return
    m2, m3, m4 = st.columns(3)
    ok_e = np.all(np.isfinite(res.z_e))
    if np.isfinite(res.T_e):
        m2.metric("Eutectic T", f"{res.T_e - 273.15:.1f} °C", f"{res.T_e:.1f} K", delta_color="off")
    else:
        m2.metric("Eutectic T", f"< {res.T_e_upper - 273.15:.0f} °C" if np.isfinite(res.T_e_upper) else "—")
    m3.metric("Ideal eutectic", f"{res.T_e_ideal - 273.15:.1f} °C", ratio_label(res.z_e_ideal), delta_color="off")
    T_ref = res.T_e if np.isfinite(res.T_e) else res.T_e_upper
    dep = res.T_e_ideal - T_ref
    m4.metric("Depression", (f"{dep:.1f} K" if np.isfinite(res.T_e) else f"> {dep:.0f} K") if np.isfinite(dep)
              else "—", help="Positive = deeper than ideal")
    if ok_e:
        st.markdown(f"**Eutectic composition** {':'.join(res.names)} = **{ratio_label(res.z_e)}** (mole fractions "
                    + ", ".join(f"{z:.3f}" for z in res.z_e) + ")")
    st.plotly_chart(ternary_figure(res, T_op), use_container_width=True)
    liq = res.liquid(T_op)
    if len(liq):
        st.success(f"{res.liquid_fraction(T_op):.0%} of the composition triangle is liquid at "
                   f"{T_op - 273.15:.1f} °C (inside the dotted line); the lowest-melting is listed first below.")
        with st.expander(f"Liquid compositions at {T_op - 273.15:.1f} °C ({len(liq)})"):
            show = liq.rename(columns={f"z_{k}": f"x {n}" for k, n in zip("ABC", res.names)})
            show = show.assign(**{"liquidus (°C)": np.where(show.below_range, np.nan, show.T_liquidus - 273.15)})
            st.dataframe(show[[f"x {n}" for n in res.names] + ["liquidus (°C)", "limiting"]]
                         .sort_values("liquidus (°C)").style.format(precision=3, na_rep="below range"),
                         hide_index=True, use_container_width=True)
    else:
        st.warning(f"No computed composition is fully liquid at {T_op - 273.15:.1f} °C.")
    for note in res.notes:
        st.info(note)
    with st.expander("Data and downloads"):
        st.dataframe(res.data.style.format(precision=4), hide_index=True, use_container_width=True)
        stem = "ternary_" + "_".join(res.names).replace(" ", "_")
        c1, c2 = st.columns(2)
        c1.download_button("Download CSV", res.data.to_csv(index=False).encode(), file_name=f"{stem}.csv",
                           key="te_dl")
        c2.download_button("Download PNG", ternary_png(res, T_op), file_name=f"{stem}.png", key="te_png")
    if ok_e:
        with st.form("te_save"):
            sp: dict[str, float] = {}
            for n, z in zip(res.names, res.z_e):
                for s_, nu in species_of(n).items():
                    sp[s_] = sp.get(s_, 0.0) + round(float(z), 4) * nu
            nm = st.text_input("Save the eutectic composition as a solvent",
                               value=f"{':'.join(res.names)} ({ratio_label(res.z_e)})")
            if st.form_submit_button("Save as solvent"):
                try:
                    lib.add_solvent(nm, {k: v for k, v in sp.items() if v > 0})
                    st.success(f"Saved “{nm}” to your solvents.")
                except ValueError as exc:
                    st.error(str(exc))
    st.caption("Ideal-solid SLE with immiscible solids (no co-crystals or solid solutions). ln γ is computed on a "
               "temperature grid and interpolated in 1/T; along each edge the result matches the binary diagram "
               "within ~0.5 K.")


# --------------------------------------------------------------------------- jobs
@st.fragment(run_every=3)
def jobs_panel():
    jobs = job_manager().jobs()
    if not jobs:
        st.info("No jobs yet.")
        return
    for j in jobs:
        icon = {"queued": "🕒", "running": "⏳", "done": "✅", "failed": "❌"}[j.status]
        with st.expander(f"{icon} #{j.id} {j.label} – {j.status} ({j.elapsed:.0f} s)", expanded=j.status == "running"):
            if j.error:
                st.error(j.error[:2000])
            st.code("\n".join(j.log[-40:]) or "(waiting)", language="text")


def page_jobs():
    st.title("Jobs")
    if st.button("Clear finished"):
        job_manager().clear_finished()
    jobs_panel()


# --------------------------------------------------------------------------- settings
def page_settings():
    st.title("Settings")
    s = settings()
    with st.form("settings"):
        home = st.text_input("ORCA folder", value=str(s.orca_home or ""))
        mode = st.selectbox("Backend", ["auto", "native", "docker"], index=["auto", "native", "docker"].index(s.orca_mode))
        ncores = st.number_input("Cores per ORCA job", 1, 64, int(s.ncores))
        image = st.text_input("Docker image", value=s.docker_image)
        solv = st.text_input("COSMORS helper solvent", value=s.cosmors_solvent,
                             help="ORCA's !COSMORS keyword needs a solvent; only the solute σ-surface is used.")
        if st.form_submit_button("Save"):
            s.orca_home = Path(home).expanduser() if home else None
            s.orca_mode, s.ncores, s.docker_image, s.cosmors_solvent = mode, int(ncores), image, solv
            st.success(f"Saved to {s.save()}")
    st.write(f"Resolved backend: **{s.resolved_mode}** · data folder: `{s.home}`")
    c1, c2, c3 = st.columns(3)
    if c1.button("Check backend"):
        st.session_state.status = orca.check_backend(s)
    if c2.button("Run test calculation"):
        with st.spinner("Running HF/def2-SVP water …"):
            st.session_state.status = orca.check_backend(s, quick_run=True)
    if s.resolved_mode == "docker" and c3.button("Build Docker runtime image"):
        with st.spinner("docker build …"):
            try:
                orca.build_docker_image(s, log=st.write)
            except Exception as exc:  # noqa: BLE001
                st.error(str(exc))
    if stt := st.session_state.get("status"):
        (st.success if stt.ok else st.error)("\n\n".join(stt.messages))


# --------------------------------------------------------------------------- NADES workflow
WF_STEPS = ["1 Solute", "2 Components", "3 Properties", "4 Optimise", "5 Results"]
WF_PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
WF_TEMPLATES = {
    "ChCl + polyol": ["choline chloride", "glycerol", "urea"],
    "ChCl + acid": ["choline chloride", "lactic acid", "malic acid"],
    "Betaine": ["betaine", "glycerol", "lactic acid"],
    "Hydrophobic": ["thymol", "menthol", "decanoic acid"],
}


def wf_state() -> dict:
    if "wf" not in st.session_state:
        st.session_state.wf = {"step": 0, "solute": None, "smiles": {}, "components": [], "bounds": {},
                               "result": None}
    return st.session_state.wf


def wf_smiles(name: str, wf: dict) -> str | None:
    """SMILES for a species that still has to be computed (user input first, then presets)."""
    return wf["smiles"].get(name) or presets()["species"].get(name, {}).get("smiles")


def wf_ready(wf: dict, lib: Library) -> dict[int, str]:
    """Reasons why a step cannot be entered yet (empty = ready)."""
    from osnades import designer as D

    why: dict[int, str] = {}
    if not wf["solute"]:
        why[1] = why[2] = why[3] = why[4] = "Choose a solute first."
        return why
    if not wf["components"]:
        why[2] = why[3] = why[4] = "Choose at least one solvent component."
        return why
    comps = [D.make_component(c, lib) for c in wf["components"]]
    miss = D.missing_species(wf["solute"], comps, lib)
    if miss:
        why[3] = why[4] = "Compute the σ-surfaces of: " + ", ".join(miss)
    elif not (lib.get(wf["solute"]).dh_fus and lib.get(wf["solute"]).t_fus):
        why[3] = why[4] = f"Enter or look up ΔH_fus and T_m of {wf['solute']}."
    if not wf.get("result"):
        why.setdefault(4, "Run the optimisation first.")
    return why


def wf_nav(wf: dict, why: dict[int, str]) -> None:
    cols = st.columns(len(WF_STEPS))
    for i, (c, label) in enumerate(zip(cols, WF_STEPS)):
        cur = i == wf["step"]
        if c.button(label, key=f"wf_nav_{i}", type="primary" if cur else "secondary", use_container_width=True,
                    disabled=i in why and not cur, help=why.get(i)):
            wf["step"] = i
            st.rerun()


def wf_footer(wf: dict, why: dict[int, str]) -> None:
    st.markdown("---")
    a, _, b = st.columns([1, 4, 1])
    if wf["step"] > 0 and a.button("← Back", key="wf_back", use_container_width=True):
        wf["step"] -= 1
        st.rerun()
    nxt = wf["step"] + 1
    if nxt < len(WF_STEPS) and nxt != 4:
        if b.button("Next →", key="wf_next", type="primary", use_container_width=True, disabled=nxt in why,
                    help=why.get(nxt)):
            wf["step"] = nxt
            st.rerun()
        if nxt in why:
            b.caption(why[nxt])


def wf_step_solute(wf: dict, lib: Library) -> None:
    st.subheader("Which compound do you want to dissolve?")
    neutral = [s.name for s in lib.all() if s.charge == 0]
    pre = presets()["species"]
    mode = st.radio("Solute", ["From the library", "Preset compound (not computed yet)", "New compound (SMILES)"],
                    horizontal=True, key="wf_sol_mode")
    if mode == "From the library":
        idx = neutral.index(wf["solute"]) if wf["solute"] in neutral else \
            next((i for i, n in enumerate(neutral) if n == "paracetamol"), 0)
        name = st.selectbox("Compound", neutral, index=idx, key="wf_sol_lib",
                            format_func=lambda n: f"{n}  ({fusion_status(lib.get(n))} fusion data)")
    elif mode == "Preset compound (not computed yet)":
        opts = [n for n in pre if not lib.has(n) and n not in ("choline", "chloride")]
        if not opts:
            st.info("All preset compounds are already in the library.")
            return
        name = st.selectbox("Compound", opts, key="wf_sol_pre")
        st.code(pre[name]["smiles"])
    else:
        a, b = st.columns([1, 2])
        name = a.text_input("Name", key="wf_sol_name").strip()
        smi = b.text_input("SMILES", key="wf_sol_smi", placeholder="CC(=O)Nc1ccc(O)cc1").strip()
        if smi:
            try:
                from osnades.structures import smiles_info

                info = smiles_info(smi)
                st.caption(f"{info['formula']} · M = {info['molar_mass']:.2f} g/mol · charge {info['charge']}")
                if info["charge"] or info["n_fragments"] > 1:
                    st.error("The solute must be a single neutral molecule.")
                    return
            except Exception as exc:  # noqa: BLE001
                st.error(str(exc))
                return
        if name and smi:
            if lib.has(name) or lib.clash(name):
                st.error(f"“{name}” is already in the library – pick it from the library instead.")
                return
            wf["smiles"][name] = smi
        else:
            name = None
    if name and name != wf["solute"]:
        if name in wf["components"]:
            st.error(f"{name} is one of the solvent components – remove it there first.")
            return
        wf["solute"], wf["result"] = name, None
    if wf["solute"]:
        s = wf["solute"]
        if lib.has(s):
            sp = lib.get(s)
            st.success(f"Solute: **{s}** · {sp.formula} · M = {sp.molar_mass:.2f} g/mol · σ-surface ✓ · "
                       f"fusion data: {fusion_status(sp)}")
        else:
            st.info(f"Solute: **{s}** – its σ-surface will be computed with ORCA in step 3.")


def wf_step_components(wf: dict, lib: Library) -> None:
    from osnades import designer as D

    st.subheader("Which components can the solvent contain?")
    st.caption("Pick the hydrogen-bond acceptors (e.g. choline chloride, betaine), donors (polyols, amides, acids, "
               "sugars, terpenes) and optionally water. The optimiser varies their mole fractions within the bounds. "
               "Salts are weighed in as one unit and split into ions by COSMO-RS.")
    choices = [c for c in D.component_choices(lib) if c != wf["solute"]]
    for n in wf["smiles"]:
        if n not in choices and n != wf["solute"]:
            choices.append(n)
    t = st.columns(len(WF_TEMPLATES) + 1)
    for c, (label, names) in zip(t, WF_TEMPLATES.items()):
        if c.button(label, key=f"wf_tpl_{label}", use_container_width=True):
            wf["components"] = [n for n in names if n in choices]
            wf["result"] = None
            st.session_state["wf_comps"] = wf["components"]
            st.rerun()
    if "wf_comps" not in st.session_state:
        st.session_state["wf_comps"] = [c for c in wf["components"] if c in choices]
    sel = st.multiselect("Components", choices, key="wf_comps",
                         format_func=lambda n: n + ("" if all(lib.has(s) for s in D.species_of(n))
                                                    else "  (to compute)"))
    with st.expander("Add a new component (SMILES)"):
        a, b, c = st.columns([1, 2, 1])
        nm = a.text_input("Name", key="wf_new_comp_name").strip()
        smi = b.text_input("SMILES", key="wf_new_comp_smi").strip()
        if c.button("Add", key="wf_new_comp_add", disabled=not (nm and smi)):
            try:
                from osnades.structures import smiles_info

                info = smiles_info(smi)
                if info["charge"] or info["n_fragments"] > 1:
                    st.error("Add neutral molecules here; salts must exist as preset salts or as ion species.")
                elif lib.has(nm) or nm in choices:
                    st.error(f"{nm} already exists – select it in the list.")
                else:
                    wf["smiles"][nm] = smi
                    st.session_state["wf_comps"] = sel + [nm]
                    st.rerun()
            except Exception as exc:  # noqa: BLE001
                st.error(str(exc))
    if sel != wf["components"]:
        wf["components"], wf["result"] = list(sel), None
    if not sel:
        return
    st.markdown("**Mole-fraction bounds** (component basis; a salt counts as one unit)")
    rows = []
    for n in sel:
        default = (0.0, 0.3) if n == "water" else ((0.1, 0.5) if len(D.species_of(n)) > 1 else (0.0, 1.0))
        lo, hi = wf["bounds"].get(n, default)
        rows.append({"component": n, "min": lo, "max": hi})
    ed = st.data_editor(pd.DataFrame(rows), hide_index=True, key=f"wf_bounds_{'|'.join(sel)}",
                        disabled=["component"], use_container_width=False,
                        column_config={"min": st.column_config.NumberColumn(min_value=0.0, max_value=1.0, step=0.05),
                                       "max": st.column_config.NumberColumn(min_value=0.0, max_value=1.0, step=0.05)})
    bad = []
    for r in ed.to_dict("records"):
        lo, hi = float(r["min"] or 0), float(r["max"] if r["max"] is not None else 1)
        if lo > hi:
            bad.append(r["component"])
        if (lo, hi) != wf["bounds"].get(r["component"]):
            wf["bounds"][r["component"]] = (lo, hi)
            wf["result"] = None
    if bad:
        st.error("min > max for " + ", ".join(bad))
    if sum(wf["bounds"][n][0] for n in sel) > 1:
        st.error("The minimum fractions add up to more than 1.")
    if any(len(D.species_of(n)) > 1 for n in sel):
        st.caption("Salts default to 0.1–0.5: COSMO-RS over-stabilises ionic mixtures and would otherwise favour "
                   "unrealistically salt-rich (solid) compositions. Typical NADES use 1:1 to 1:4 salt:HBD.")
    if "water" not in sel:
        st.caption("Tip: add water (e.g. 0–0.3) – a little water lowers viscosity and often changes solubility a lot.")


def wf_step_properties(wf: dict, lib: Library) -> None:
    from osnades import designer as D
    from osnades import fusion

    st.subheader("Required properties")
    st.caption("Every species needs a COSMO σ-surface (ORCA). The solute needs ΔH_fus and T_m for the solubility; "
               "components with fusion data are also checked for crystallising out of the solvent.")
    solute = wf["solute"]
    comps = [D.make_component(c, lib) for c in wf["components"]]
    jobs = {j.label: j for j in job_manager().jobs()}
    rows, need = [], []
    entries = [(solute, "solute", [solute])] + [(c.name, "component", list(c.species)) for c in comps]
    for name, role, species in entries:
        for s in species:
            have = lib.has(s)
            j = jobs.get(f"Compute {s}")
            state = "✓" if have else (f"ORCA: {j.status}" if j and j.status in ("queued", "running", "failed")
                                      else "missing")
            if not have and not (j and j.status in ("queued", "running")):
                need.append(s)
        if role == "solute":
            sp = lib.get(name) if lib.has(name) else None
            fus = fusion_status(sp) if sp else "after ORCA"
            dh, tm, src = (sp.dh_fus, sp.t_fus, sp.fusion_source or sp.notes) if sp else (None, None, "")
        else:
            c = next(c for c in comps if c.name == name)
            dh, tm, src = c.dh_fus, c.t_fus, c.fusion_source
            fus = ("estimate" if "estimate" in (src or "") or "Walden" in (src or "") else "available") if c.has_fusion \
                else "none (not checked for crystallising)"
        rows.append({"compound": name, "role": role, "COSMO species": ", ".join(species),
                     "σ-surface": "✓" if all(lib.has(s) for s in species) else
                     ", ".join(f"{s}: {('✓' if lib.has(s) else (jobs.get(f'Compute {s}').status if jobs.get(f'Compute {s}') else 'missing'))}"
                               for s in species),
                     "ΔH_fus (kJ/mol)": round(dh / 1000, 2) if dh else None, "T_m (K)": tm, "fusion data": fus,
                     "source": src or ""})
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True,
                 column_config={"source": st.column_config.TextColumn(width="medium")})

    # σ-surfaces
    nosmi = [s for s in need if not wf_smiles(s, wf)]
    todo = [s for s in need if wf_smiles(s, wf)]
    if nosmi:
        st.error("No structure known for: " + ", ".join(nosmi) + " – add them in the Compound library.")
    if todo:
        st.warning(f"{len(todo)} σ-surface(s) to compute: {', '.join(todo)} (≈1–3 min each under Docker).")
        if st.button(f"Compute {len(todo)} missing σ-surface(s) with ORCA", type="primary", key="wf_compute"):
            for s in todo:
                queue_compute(s, smiles=wf_smiles(s, wf))
            st.rerun()
    running = [j for j in job_manager().active() if j.label.removeprefix("Compute ") in
               {s for _, _, sps in entries for s in sps}]
    if running:
        wf_jobs_live(set(j.id for j in running))
    failed = [jobs[f"Compute {s}"] for _, _, sps in entries for s in sps
              if not lib.has(s) and f"Compute {s}" in jobs and jobs[f"Compute {s}"].status == "failed"]
    for j in failed:
        st.error(f"{j.label} failed: {j.error[:300]} – see the Jobs page; you can retry above.")

    # fusion data
    if lib.has(solute):
        sp = lib.get(solute)
        st.markdown(f"**Fusion data of the solute ({solute})**")
        lookup_targets = [solute] if not (sp.dh_fus and sp.t_fus) else []
        lookup_targets += [c.name for c in comps if not c.has_fusion and lib.has(c.name)]
        if lookup_targets and st.button(f"Look up missing fusion data ({', '.join(lookup_targets)})",
                                        key="wf_fus_lookup"):
            msgs = []
            with st.spinner("Searching Cheméo and the NIST WebBook…"):
                for n in lookup_targets:
                    lib.lookup_fusion(n, log=msgs.append)
            st.session_state["wf_fus_msgs"] = msgs
            st.rerun()
        for m in st.session_state.pop("wf_fus_msgs", []):
            st.caption(m)
        with st.form(f"wf_fus_{solute}"):
            a, b, c, d = st.columns([1, 1, 1, 1])
            dh = a.number_input("ΔH_fus (kJ/mol)", value=float(sp.dh_fus / 1000) if sp.dh_fus else None,
                                min_value=0.1, placeholder="required")
            tm = b.number_input("T_m (K)", value=float(sp.t_fus) if sp.t_fus else None, min_value=1.0,
                                placeholder="required")
            dcp = c.number_input("ΔCp_fus (J/mol/K)", value=float(sp.dcp_fus or 0.0))
            d.markdown("&nbsp;")
            if d.form_submit_button("Save"):
                if dh and tm:
                    changed = (dh * 1000, tm) != (sp.dh_fus, sp.t_fus)
                    lib.update(solute, dh_fus=dh * 1000, t_fus=tm, dcp_fus=dcp or None,
                               fusion_source="entered by user" if changed else sp.fusion_source)
                    wf["result"] = None
                    st.rerun()
        fusion_caption(sp)


@st.fragment(run_every="5s")
def wf_jobs_live(ids: set[int]) -> None:
    jobs = [j for j in job_manager().jobs() if j.id in ids]
    for j in jobs:
        icon = {"queued": "⏳", "running": "⚙️", "done": "✅", "failed": "❌"}[j.status]
        st.caption(f"{icon} {j.label} – {j.status} ({j.elapsed:.0f} s) · {j.log[-1][11:] if j.log else ''}")
    if jobs and all(j.status in ("done", "failed") for j in jobs):
        st.rerun(scope="app")


def wf_step_optimise(wf: dict, lib: Library) -> None:
    from osnades import designer as D

    st.subheader("Conditions and optimisation")
    a, b = st.columns([2, 1])
    t_lo, t_hi = a.slider("Temperature range (°C)", -20, 120, key="wf_trange",
                          **({} if "wf_trange" in st.session_state else {"value": wf.get("trange", (15, 45))}))
    nT = b.number_input("Temperature points", 1, 9, wf.get("nT", 3), key="wf_nT",
                        help="Solubility is evaluated at these temperatures (evenly spaced).")
    temps = np.linspace(t_lo, t_hi, int(nT)) + 273.15 if nT > 1 else np.array([t_lo + 273.15])
    obj = st.radio("Goal", list(D.OBJECTIVES), format_func=D.OBJECTIVES.get, key="wf_obj",
                   index=list(D.OBJECTIVES).index(wf.get("objective", "mean_x")))
    c1, c2, c3 = st.columns(3)
    step = c1.select_slider("Grid step (mole fraction)", [0.25, 0.2, 0.1, 0.05], value=wf.get("grid", 0.1),
                            key="wf_step")
    rounds = c2.number_input("Local refinement rounds", 0, 4, wf.get("rounds", 2), key="wf_rounds",
                             help="Each round halves the step around the 3 best compositions.")
    top_n = c3.number_input("Candidates to report", 1, 20, wf.get("top_n", 6), key="wf_topn")
    liquid = st.checkbox("Solvent must be liquid over the whole range", value=wf.get("liquid", True), key="wf_liq",
                         help="Rejects compositions where a component would crystallise (needs its fusion data).")
    comps = [D.make_component(c, lib, *wf["bounds"].get(c, (0.0, 1.0))) for c in wf["components"]]
    st.caption("Bounds: " + " · ".join(f"{c.name} {c.lo:g}–{c.hi:g}" for c in comps) + " (change them in step 2)")
    n_pts = len(D.lattice(comps, step))
    salts = any(len(c.species) > 1 for c in comps)
    est = n_pts * len(temps) * 15 * (0.009 if salts else 0.005) + int(top_n) * len(temps) * 0.5
    st.caption(f"{n_pts} compositions on the grid × {len(temps)} temperatures · estimated time ≈ {est:.0f} s"
               + (" (+ refinement)" if rounds else ""))
    if n_pts > 3000:
        st.warning("That is a large grid – use a coarser step and let the local refinement zoom in.")
    if st.button("Find the best composition", type="primary", key="wf_run", disabled=n_pts == 0):
        sp = lib.get(wf["solute"])
        bar = st.progress(0.0, text="Starting…")
        try:
            des = D.Designer(engine(), wf["solute"], comps, sp.dh_fus, sp.t_fus, sp.dcp_fus or 0.0)
            res = des.run(temps, step=float(step), objective=obj, require_liquid=liquid, top_n=int(top_n),
                          local_rounds=int(rounds), progress=lambda f, m: bar.progress(min(1.0, f), text=m))
        except Exception as exc:  # noqa: BLE001
            bar.empty()
            st.error(f"Optimisation failed: {exc}")
            return
        bar.empty()
        wf.update(result=res, trange=(t_lo, t_hi), nT=int(nT), objective=obj, grid=step, rounds=int(rounds),
                  top_n=int(top_n), liquid=liquid)
        wf["step"] = 4
        st.rerun()
    if n_pts == 0:
        st.error("No composition satisfies the bounds – go back and widen them.")


def wf_step_results(wf: dict, lib: Library) -> None:
    from osnades import designer as D

    res = wf["result"]
    comps = res.components
    if res.best.empty:
        st.error("No valid composition found.")
        for n in res.notes:
            st.warning(n)
        return
    best = res.best.iloc[0]
    z = np.array([best[f"z_{c.name}"] for c in comps])
    mm = np.array([c._mm for c in comps])
    w = z * mm / (z * mm).sum()
    st.subheader(f"Best solvent for {res.solute}")
    st.markdown(f"### {best['name']}")
    st.caption(D.OBJECTIVES[res.objective] + f" · {res.temps[0] - 273.15:.0f}–{res.temps[-1] - 273.15:.0f} °C")
    cols = st.columns(min(4, len(res.temps)))
    sp = lib.get(res.solute)
    for c, T in zip(cols * 3, res.temps[:len(cols)] if len(res.temps) <= 4 else res.temps[[0, len(res.temps) // 2, -1]]):
        xT = best[f"x@{T - 273.15:.1f}°C"]
        gT = best[f"g/kg@{T - 273.15:.1f}°C"]
        xid = float(D.ideal_curve([T], sp.dh_fus, sp.t_fus, sp.dcp_fus or 0.0)[0])
        c.metric(f"{T - 273.15:.0f} °C", f"x = {xT:.3g}", f"{gT:.4g} g/kg · {xT / xid:.2g}× ideal", delta_color="off")
    tbl = pd.DataFrame({"component": [c.name for c in comps], "mole fraction": z.round(4),
                        "mass fraction (wt%)": (100 * w).round(2),
                        "molar ratio": D.ratio_text(z).split(" : ")})
    st.dataframe(tbl[tbl["mole fraction"] > 0], hide_index=True, use_container_width=False)
    ls = wf.get("liquid_search")
    if ls is not None and ls.feasible and ls.result is res:
        st.success(f"Liquid search: {len(ls.liquid_points)} compositions are liquid at {ls.T_min - 273.15:.1f} °C. "
                   f"The lowest-melting one, {ls.label}, melts at {ls.melt_text()}. The optimisation was re-run "
                   "with them, and the best liquid composition is shown.")
    if not best["liquid_all_T"]:
        st.warning("This composition is predicted to partly crystallise in the range (no fully liquid composition "
                   "was found).")
        wf_liquid_search(wf, lib, res)
    at_bound = [f"{c.name} ({'upper' if zi >= c.hi - 1e-6 else 'lower'} bound {c.hi if zi >= c.hi - 1e-6 else c.lo:g})"
                for c, zi in zip(comps, z) if (zi >= c.hi - 1e-6 and c.hi < 1) or (zi <= c.lo + 1e-6 and c.lo > 0)]
    if at_bound:
        st.warning("The optimum sits on a bound: " + ", ".join(at_bound) + ". The true optimum may lie beyond it; "
                   "widen the bound in step 2 if that composition is practical (for salts, a salt-rich optimum is "
                   "usually a COSMO-RS artefact).")
    for n in res.notes:
        if not n.startswith("Liquid search:"):  # shown above as a success message
            st.info(n)

    tab_t, tab_c, tab_p = st.tabs(["Solubility vs temperature", "Composition map", "Phase diagram"])
    with tab_p:
        wf_phase_diagram(wf, res, z)
    with tab_t:
        fig = go.Figure()
        T_line = np.linspace(res.temps[0], res.temps[-1], 30) if len(res.temps) > 1 else res.temps
        fig.add_trace(go.Scatter(x=T_line - 273.15, y=D.ideal_curve(T_line, sp.dh_fus, sp.t_fus, sp.dcp_fus or 0.0),
                                 name="ideal solution", mode="lines", line=dict(color="#5f5e58", dash="dash")))
        for i, (rank, d) in enumerate(res.curves.groupby("rank")):
            if rank > 6:
                break
            d = d.sort_values("T")
            fig.add_trace(go.Scatter(x=d["T"] - 273.15, y=d["x"], name=f"#{rank} {d['name'].iloc[0]}",
                                     mode="lines+markers", line=dict(color=WF_PALETTE[i % 8], width=3 if rank == 1 else 1.5)))
        fig.update_layout(title="Solubility of the best candidates", height=460, hovermode="x unified",
                          xaxis_title="Temperature (°C)", yaxis_title=f"mole fraction of {res.solute}",
                          legend=dict(orientation="h", y=-0.25), margin=dict(t=40, l=10, r=10))
        st.plotly_chart(fig, use_container_width=True)
    with tab_c:
        scr = res.screen
        Tm = res.temps[len(res.temps) // 2]
        col = f"x@{Tm:.2f}"
        live = scr[np.isfinite(scr[col])]
        nz = [c for c in comps if (live[f"z_{c.name}"] > 0).any()]
        if len(comps) == 3:
            fig2 = px.scatter_ternary(live, a=f"z_{comps[0].name}", b=f"z_{comps[1].name}", c=f"z_{comps[2].name}",
                                      color=col, symbol="liquid_all_T", color_continuous_scale="Viridis",
                                      symbol_map={True: "circle", False: "x"}, hover_data=["ratio"])
            fig2.update_layout(title=f"Screened compositions – colour: x at {Tm - 273.15:.0f} °C, ✕ = not liquid",
                               height=560, margin=dict(t=80, l=60, r=60), showlegend=False)
            st.plotly_chart(fig2, use_container_width=True)
        elif len(comps) == 2:
            d2 = live.sort_values(f"z_{comps[0].name}")
            fig2 = px.line(d2, x=f"z_{comps[0].name}", y=col, markers=True)
            fig2.update_traces(line_color=WF_PALETTE[0])
            fig2.update_layout(title=f"x at {Tm - 273.15:.0f} °C vs composition", height=420,
                               xaxis_title=f"mole fraction {comps[0].name}", yaxis_title="x")
            st.plotly_chart(fig2, use_container_width=True)
        else:
            dims = [f"z_{c.name}" for c in nz]
            fig2 = px.parallel_coordinates(live.nlargest(min(len(live), 300), col), dimensions=dims + [col],
                                           color=col, color_continuous_scale="Viridis")
            fig2.update_layout(title=f"Top screened compositions (x at {Tm - 273.15:.0f} °C)", height=420,
                               margin=dict(t=70, l=40, r=40))
            st.plotly_chart(fig2, use_container_width=True)
        st.caption(f"{len(scr)} compositions screened · {int(scr['liquid_all_T'].sum())} liquid over the whole range.")

    st.markdown("**All reported candidates** (exact SLE solubility)")
    st.dataframe(res.best.drop(columns=["score"]), hide_index=True, use_container_width=True)
    a, b, c = st.columns([1, 1, 2])
    a.download_button("Candidates (CSV)", res.best.to_csv(index=False).encode(), f"{res.solute}_best_nades.csv")
    b.download_button("Screening grid (CSV)", res.screen.to_csv(index=False).encode(), f"{res.solute}_screen.csv")
    with c.form("wf_save_solvent"):
        nm = st.text_input("Save the best composition as a solvent", value=f"Best for {res.solute}")
        if st.form_submit_button("Save as solvent"):
            try:
                lib.add_solvent(nm, res.best_system(0).components, overwrite=True)
                st.success(f"Saved “{nm}” – available on the Solubility, Temperature curves and NADES pages.")
            except ValueError as exc:
                st.error(str(exc))


def wf_liquid_search(wf: dict, lib: Library, res) -> None:
    """Look for a fully liquid composition; if none exists, say what to change."""
    from osnades import designer as D

    ls = wf.get("liquid_search")
    if ls is None or ls.result is not None or getattr(ls, "_for", None) is not res:
        st.markdown("The grid may have missed a narrow liquid region near the eutectic. A finer search looks for "
                    "the lowest-melting composition within your bounds and re-runs the optimisation with every "
                    f"composition that is liquid at {res.temps[0] - 273.15:.1f} °C.")
        if st.button("Find a fully liquid composition", type="primary", key="wf_find_liquid"):
            sp = lib.get(res.solute)
            bar = st.progress(0.0, text="Starting…")
            try:
                des = D.Designer(engine(), res.solute, res.components, sp.dh_fus, sp.t_fus, sp.dcp_fus or 0.0)
                ls = des.find_liquid(res.temps, step=float(wf.get("grid", 0.1)), objective=res.objective,
                                     top_n=int(wf.get("top_n", 6)), local_rounds=int(wf.get("rounds", 2)),
                                     progress=lambda f, m: bar.progress(min(1.0, f), text=m))
            except Exception as exc:  # noqa: BLE001
                bar.empty()
                st.error(f"Search failed: {exc}")
                return
            bar.empty()
            ls._for = res
            wf["liquid_search"] = ls
            if ls.feasible:
                wf["result"] = ls.result
                wf.pop("phase", None)
            st.rerun()
        return
    # searched, nothing liquid: explain and offer the fixes
    msg = ls.advice()
    st.error(msg[0])
    for m in msg[1:]:
        st.markdown("- " + m)
    t_new = ls.suggested_T_min()
    a, b = st.columns(2)
    if t_new is not None:
        t_lo, t_hi = wf.get("trange", (15, 45))
        t_new = min(t_new, 115)
        t_hi = min(120, max(t_hi, t_new + 5))
        if a.button(f"Start the range at {t_new} °C", type="primary", key="wf_raise_T",
                    help="Goes back to step 4 with the new range; run the optimisation again there."):
            wf["trange"] = (t_new, t_hi)
            st.session_state["wf_trange"] = (t_new, t_hi)
            wf["step"] = 3
            st.rerun()
    if b.button("Change components or bounds", key="wf_to_comps"):
        wf["step"] = 1
        st.rerun()


def wf_phase_diagram(wf: dict, res, z: np.ndarray) -> None:
    """Solid–liquid phase diagram of the solvent components (2 or 3) with the chosen composition marked."""
    from osnades import designer as D
    from osnades.eutectic import plotly_figure, ternary_cost, ternary_figure

    comps, zc = D.phase_components(res.components, z)
    if len(comps) not in (2, 3):
        st.info(f"The phase diagram needs 2 or 3 components; the best solvent has {len(comps)}. "
                "Use the NADES designer page for any 2 or 3 of them.")
        return
    missing = [c.name for c in comps if not c.has_fusion]
    if missing:
        st.warning("Fusion data (ΔH_fus, T_m) are missing for " + ", ".join(missing) + " – add them in step 3.")
        return
    T_cold = float(res.temps[0])
    names = [c.name for c in comps]
    st.markdown(f"Melting behaviour of the solute-free solvent **{' + '.join(names)}**: where it is liquid, its "
                f"eutectic, and where the chosen composition sits. Dotted line: {T_cold - 273.15:.0f} °C, the "
                "coldest temperature of your range.")
    key = (tuple(names), tuple((c.dh_fus, c.t_fus, c.dcp_fus) for c in comps), tuple(np.round(zc, 4)))
    ions = any(len(c.species) > 1 for c in comps)
    if len(comps) == 3:
        step = st.select_slider("Grid step", [0.2, 0.1, 0.05], value=0.1, key="wf_pd_step")
        label = f"Compute ternary phase diagram (≈ {ternary_cost(step, ions):.0f} s)"
    else:
        step, label = 0.1, "Compute phase diagram"
    if st.button(label, key="wf_pd_run", type="primary" if wf.get("phase", (None,))[0] != key else "secondary"):
        prog = st.progress(0.0, text="Solving liquidus…")
        try:
            sp_ = D.solvent_phase_diagram(engine(), comps, zc, step=step,
                                          progress=lambda f: prog.progress(min(1.0, f), text="Solving liquidus…"))
            wf["phase"] = (key, sp_)
        except Exception as exc:  # noqa: BLE001
            st.error(str(exc))
        prog.empty()
    ph = wf.get("phase")
    if not ph or ph[0] != key or len(ph) != 2:
        return
    sp_ = ph[1]
    pd_ = sp_.diagram
    m1, m2, m3 = st.columns(3)
    m1.metric("Chosen solvent melts at", sp_.melt_text(),
              f"{sp_.T_melt:.1f} K" if np.isfinite(sp_.T_melt) else None, delta_color="off")
    T_e, _, comp_e = sp_.eutectic_text().partition(" at ")
    m2.metric("Eutectic", T_e, (comp_e.split(" = ")[-1] if sp_.ternary else comp_e) or None, delta_color="off",
              help=sp_.eutectic_text())
    if sp_.ternary:
        m3.metric(f"Liquid at {T_cold - 273.15:.0f} °C", f"{pd_.liquid_fraction(T_cold):.0%}",
                  help="Share of all compositions of these three components that is liquid")
        st.plotly_chart(ternary_figure(pd_, T_cold, marks=[(sp_.z, "chosen solvent")]), use_container_width=True)
    else:
        win = pd_.liquid_window(T_cold)
        m3.metric(f"Liquid at {T_cold - 273.15:.0f} °C", f"x {win[0]:.2f}–{win[1]:.2f}" if win else "none")
        fig = plotly_figure(pd_, T_cold)
        fig.add_vline(x=float(sp_.z[0]), line=dict(color=WF_PALETTE[0], width=2),
                      annotation_text="chosen solvent", annotation_position="top right")
        st.plotly_chart(fig, use_container_width=True)
    liquid = sp_.liquid_at(T_cold)
    if liquid is True:
        st.success(f"The chosen composition is predicted to be liquid over the whole range "
                   f"(≥ {T_cold - 273.15:.0f} °C).")
    elif liquid is False:
        st.warning(f"The chosen composition is predicted to start crystallising at {sp_.melt_text()}, above "
                   f"the coldest temperature of your range ({T_cold - 273.15:.0f} °C).")
    for n in pd_.notes:
        st.caption("• " + n)
    wf_phase_downloads(res, sp_, T_cold)


def wf_phase_downloads(res, sp_, T_cold: float) -> None:
    """CSV / PNG / summary downloads of the workflow phase diagram."""
    from osnades.eutectic import matplotlib_png, ratio_label, ternary_png

    pd_ = sp_.diagram
    stem = f"{res.solute}_{'_'.join(sp_.names)}_phase".replace(" ", "_")
    if sp_.ternary:
        png = ternary_png(pd_, T_cold, marks=[(sp_.z, "chosen solvent")])
        data = pd_.data.rename(columns={f"z_{k}": f"x_{n}" for k, n in zip("ABC", sp_.names)})
    else:
        png = matplotlib_png(pd_, T_cold, marks=[(float(sp_.z[0]), "chosen solvent")])
        data = pd_.data.rename(columns={"x_A": f"x_{sp_.names[0]}"})
    rows = [("solute", res.solute), ("solvent components", " + ".join(sp_.names)),
            ("chosen composition", ratio_label(sp_.z, sp_.names)),
            *[(f"x {n}", f"{z:.4f}") for n, z in zip(sp_.names, sp_.z)],
            ("chosen solvent melts at", sp_.melt_text()), ("eutectic", sp_.eutectic_text()),
            ("coldest temperature of range", f"{T_cold - 273.15:.1f} °C"),
            ("liquid over the whole range", {True: "yes", False: "no", None: "unknown"}[sp_.liquid_at(T_cold)])]
    if sp_.ternary:
        rows.append((f"share of triangle liquid at {T_cold - 273.15:.0f} °C", f"{pd_.liquid_fraction(T_cold):.0%}"))
    rows += [("note", n) for n in pd_.notes]
    summary = pd.DataFrame(rows, columns=["item", "value"])
    st.markdown("**Export**")
    a, b, c = st.columns(3)
    a.download_button("Diagram (PNG)", png, f"{stem}.png", mime="image/png", key="wf_pd_png",
                      use_container_width=True)
    b.download_button("Liquidus data (CSV)", data.to_csv(index=False).encode(), f"{stem}.csv", mime="text/csv",
                      key="wf_pd_csv", use_container_width=True)
    c.download_button("Summary (CSV)", summary.to_csv(index=False).encode(), f"{stem}_summary.csv",
                      mime="text/csv", key="wf_pd_sum", use_container_width=True)


def page_workflow():
    st.title("NADES workflow")
    st.caption("Design a multicomponent (natural) deep eutectic solvent for a solute, step by step: solute → "
               "components → properties (ORCA σ-surfaces, fusion data) → optimise mole fractions over a "
               "temperature range → results.")
    lib = library()
    wf = wf_state()
    if wf["solute"] and not lib.has(wf["solute"]) and wf["solute"] not in wf["smiles"] \
            and wf["solute"] not in presets()["species"]:
        wf["solute"] = None
    why = wf_ready(wf, lib)
    if wf["step"] in why:
        wf["step"] = max(i for i in range(wf["step"]) if i not in why) if wf["step"] else 0
    wf_nav(wf, why)
    st.markdown("")
    [wf_step_solute, wf_step_components, wf_step_properties, wf_step_optimise, wf_step_results][wf["step"]](wf, lib)
    why_now = wf_ready(wf, lib)
    if set(why_now) != set(why):  # a choice on this step unlocked/locked other steps: refresh the step bar
        st.rerun()
    wf_footer(wf, why_now)


{"NADES workflow": page_workflow, "Overview": page_overview, "Compound library": page_library, "Solubility": page_solubility, "Temperature curves": page_curves,
 "NADES designer": page_nades, "Jobs": page_jobs, "Settings": page_settings}[page]()
