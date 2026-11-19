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
from osnades.thermo import Engine, SolventSystem, presets, preset_systems

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
    ["Overview", "Compound library", "Solubility", "Temperature curves", "NADES designer", "Jobs", "Settings"],
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

1. **Compound library** – add molecules/ions (SMILES, XYZ or existing `.orcacosmo` files). ORCA computes
   their COSMO σ-surfaces.
2. **Solubility** – choose a solute, its fusion enthalpy and melting point, and screen solvents/NADES over
   temperature.
3. **NADES designer** – build custom HBA:HBD mixtures (± water), check species activity coefficients and
   binary excess Gibbs energies.
"""
    )
    lib = library()
    sp = lib.all()
    systems = preset_systems()
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
def page_library():
    st.title("Compound library")
    lib = library()
    all_sp = lib.all()

    tab_list, tab_add, tab_import, tab_sigma = st.tabs(["Species", "Add with ORCA", "Import .orcacosmo", "σ-profiles"])

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
                with st.form(f"fus_{name}"):
                    st.markdown("**Solid–liquid data** (needed when used as solute)")
                    a, b, c = st.columns(3)
                    dh = a.number_input("ΔH_fus (kJ/mol)", value=float(sp.dh_fus / 1000) if sp.dh_fus else 0.0, min_value=0.0)
                    tm = b.number_input("T_m (K)", value=float(sp.t_fus or 0.0), min_value=0.0)
                    dcp = c.number_input("ΔCp_fus (J/mol/K)", value=float(sp.dcp_fus or 0.0))
                    notes = st.text_input("Notes / source", value=sp.notes)
                    if st.form_submit_button("Save"):
                        lib.update(name, dh_fus=dh * 1000 or None, t_fus=tm or None, dcp_fus=dcp or None, notes=notes)
                        st.success("Saved")
                        st.rerun()
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
        nconf = c2.number_input("Conformers", 1, 10, 1, help="Lowest-energy RDKit conformers; all are passed to openCOSMO-RS")
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
        if st.button("Queue ORCA calculation", type="primary", disabled=not name or not (smiles or xyz)):
            if lib.has(name):
                st.error(f"{name} already exists")
            else:
                queue_compute(name, smiles=smiles, xyz=xyz, charge=charge, multiplicity=mult,
                              n_conformers=int(nconf), optimizer=optimizer, **meta)
                st.success(f"Queued {name}. Track progress on the Jobs page.")

    with tab_import:
        st.write("Import σ-surfaces computed elsewhere (e.g. the tutorial's `*.solute.orcacosmo`). "
                 "Several files = several conformers of the same species.")
        name = st.text_input("Species name", key="imp_name")
        smi = st.text_input("SMILES (optional, for molar mass)", key="imp_smi")
        files = st.file_uploader("orcacosmo file(s)", accept_multiple_files=True)
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
    return SolventSystem(label or "custom", {s: r for s, r in rows})


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
        dh = st.number_input("ΔH_fus (kJ/mol)", value=float(sp.dh_fus / 1000) if sp.dh_fus else 25.0, min_value=0.1)
        tm = st.number_input("T_m (K)", value=float(sp.t_fus or 400.0), min_value=1.0)
        dcp = st.number_input("ΔCp_fus (J/mol/K)", value=float(sp.dcp_fus or 0.0))
        method = st.radio("Method", ["iterative", "infinite_dilution"], horizontal=True,
                          help="iterative: γ at the saturated composition; infinite_dilution: γ∞ (tutorial's non-iterative)")
        water = st.slider("Add water to every solvent (wt%)", 0.0, 50.0, 0.0, 1.0)
        par = st.selectbox("Parameterisation", ["openCOSMORS24a"])
    with c2:
        presets_all = preset_systems()
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
    solutes = [s for s in lib.all() if s.dh_fus and s.t_fus]
    if not solutes:
        st.info("No solute with fusion data yet. Set ΔH_fus and T_m for a species in the Compound library.")
        return
    c1, c2, c3 = st.columns([1.2, 1, 1])
    with c1:
        solute = st.selectbox("Solute", [s.name for s in solutes], key="cv_solute")
        sp = lib.get(solute)
        a, b, c = st.columns(3)
        dh = a.number_input("ΔH_fus (kJ/mol)", value=float(sp.dh_fus / 1000), min_value=0.1, key="cv_dh")
        tm = b.number_input("T_m (K)", value=float(sp.t_fus), min_value=1.0, key="cv_tm")
        dcp = c.number_input("ΔCp_fus (J/mol/K)", value=float(sp.dcp_fus or 0.0), key="cv_dcp")
    with c2:
        t0, t1 = st.slider("Temperature range (°C)", -20, 150, (0, 80), key="cv_range")
        npts = st.slider("Points per curve", 5, 41, 17, key="cv_n")
        if t1 + 273.15 >= tm:
            st.warning("Range reaches T_m; above it the solute is liquid and SLE does not apply.")
    with c3:
        method = st.radio("Method", ["iterative", "infinite_dilution"], key="cv_method")
        water = st.slider("Add water (wt%)", 0.0, 50.0, 0.0, 1.0, key="cv_water")

    ready = [s for s in preset_systems() if not missing_species([s]) and solute not in s.components]
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
def page_nades():
    st.title("NADES designer")
    lib = library()
    if len(lib.names()) < 2:
        st.info("Add at least two species.")
        return
    t1, t2 = st.tabs(["Mixture activity", "Binary G^E"])
    with t1:
        presets_ok = [s for s in preset_systems() if not missing_species([s]) and len(s.components) > 1]
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
                    sig, area = sigma_profile(lib.paths(s)[0])
                    tot = area * f if tot is None else tot + area * f
                fig2.add_trace(go.Scatter(x=sig, y=tot, mode="lines", name="mixture σ-profile", fill="tozeroy"))
                fig2.update_layout(xaxis_title="σ (e/Å²)", yaxis_title="x-weighted p(σ) (Å²)", xaxis_range=[-0.03, 0.03], height=350)
                st.plotly_chart(fig2, use_container_width=True)
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


{"Overview": page_overview, "Compound library": page_library, "Solubility": page_solubility, "Temperature curves": page_curves,
 "NADES designer": page_nades, "Jobs": page_jobs, "Settings": page_settings}[page]()
