"""Regression tests against the FACCTs OPI openCOSMO-RS tutorial (no ORCA needed)."""
from pathlib import Path

import numpy as np
import pytest

from osnades.config import Settings
from osnades.library import Library
from osnades.sigma import descriptors
from osnades.structures import parse_xyz, smiles_info, smiles_to_conformers
from osnades.thermo import Engine, SolventSystem, ideal_rhs

FIX = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="module")
def engine(tmp_path_factory):
    s = Settings(home=tmp_path_factory.mktemp("home"), orca_home=None)
    lib = Library(s)
    lib.import_orcacosmo("water", [FIX / "water.orcacosmo"], smiles="O")
    lib.import_orcacosmo("ethanol", [FIX / "ethanol.orcacosmo"], smiles="CCO")
    lib.import_orcacosmo("paracetamol", [FIX / "paracetamol.orcacosmo"], smiles="CC(=O)Nc1ccc(O)cc1",
                         dh_fus=27100.0, t_fus=443.6)
    return Engine(lib)


def test_descriptors_match_tutorial():
    d = descriptors(FIX / "water.orcacosmo")
    assert d["area_A2"] == pytest.approx(42.897746777517476)
    assert d["volume_A3"] == pytest.approx(25.36550524969701)
    assert d["energy_dielectric_kJmol"] == pytest.approx(-32.83919025918012)
    assert d["hb_donor_moments"] == pytest.approx([5.11873012, 2.60564439, 0.85940158], rel=1e-6)


def test_reference_state_conversion(engine):
    sp, X = ["paracetamol", "water"], np.array([[0.02, 0.98]])
    direct = engine.ln_gamma(sp, X, 298.15, refst="pure_component")[0, 0]
    assert engine.ln_gamma_pure_ref(sp, X, 298.15)[0] == pytest.approx(direct, abs=1e-9)


@pytest.mark.parametrize("solvent,x_noniter", [("Water", 7.42369e-04), ("Ethanol", 1.71201e-01)])
def test_tutorial_noniterative(engine, solvent, x_noniter):
    r = engine.solubility("paracetamol", SolventSystem.parse(solvent), 298.15, method="infinite_dilution")
    # tutorial uses R = 8.314; we use the exact CODATA value -> ~0.02 % difference
    assert r.x == pytest.approx(x_noniter, rel=1e-3)


@pytest.mark.parametrize("solvent,x_iter", [("Water", 7.670e-04), ("Ethanol", 9.983e-02)])
def test_iterative_is_selfconsistent(engine, solvent, x_iter):
    sysm = SolventSystem.parse(solvent)
    r = engine.solubility("paracetamol", sysm, 298.15)
    assert r.x == pytest.approx(x_iter, rel=2e-3)
    lng = engine.ln_gamma_pure_ref(["paracetamol"] + sysm.species, [[r.x, 1 - r.x]], 298.15)[0]
    assert np.log(r.x) + lng == pytest.approx(ideal_rhs(298.15, 27100.0, 443.6), abs=1e-6)


def test_mixture_and_water_addition(engine):
    lib = engine.library
    sysm = SolventSystem.parse("ethanol:water=1:1")
    assert np.allclose(sysm.fractions, [0.5, 0.5])
    wet = SolventSystem.parse("Ethanol").with_water(10, lib)
    m_e, m_w = lib.get("ethanol").molar_mass, lib.get("water").molar_mass
    xw = wet.fractions[1]
    assert xw * m_w / (xw * m_w + (1 - xw) * m_e) == pytest.approx(0.10, rel=1e-6)
    r = engine.solubility("paracetamol", sysm, 298.15)
    assert 1e-4 < r.x < 0.2


def test_solubility_increases_with_T(engine):
    df = engine.solubility_curve("paracetamol", SolventSystem.parse("Water"), np.array([280.0, 300.0, 320.0]))
    assert np.all(np.diff(df["x"]) > 0)


def test_binary_limits(engine):
    df = engine.binary("water", "ethanol", 298.15, n=5)
    assert df["GE_RT"].iloc[0] == pytest.approx(0, abs=1e-9) and df["GE_RT"].iloc[-1] == pytest.approx(0, abs=1e-9)
    assert df["ln_gamma_water"].iloc[-1] == pytest.approx(0, abs=1e-9)


def test_structures():
    info = smiles_info("C[N+](C)(C)CCO.[Cl-]")
    assert info["n_fragments"] == 2 and info["charge"] == 0
    g = smiles_to_conformers("CCO", n_conformers=3)
    assert len(g) >= 1 and len(g[0].symbols) == 9
    g2 = parse_xyz(g[0].to_xyz())
    assert g2.formula == "C2H6O" and g2.molar_mass == pytest.approx(46.07, abs=0.01)


# ---------------------------------------------------------------- eutectic phase diagram
def test_eutectic_phase_diagram(engine):
    from osnades.eutectic import Component, EutecticSolver, _t_from_lna
    from osnades.thermo import ideal_rhs

    para = Component("paracetamol", {"paracetamol": 1.0}, 27100.0, 443.6)
    water = Component("water", {"water": 1.0}, 6009.5, 273.15)
    sol = EutecticSolver(engine, para, water)
    pd_ = sol.solve(n=11)
    d = pd_.data
    # pure-component end points
    assert d.T_A.iloc[-1] == pytest.approx(443.6)
    assert d.T_B.iloc[0] == pytest.approx(273.15)
    # neutral component activity = x γ with the pure-liquid reference
    x, T = 0.3, 350.0
    lg = engine.ln_gamma_pure_ref(["paracetamol", "water"], np.array([[x, 1 - x]]), T)[0]
    assert sol.ln_a(x, T, 0) == pytest.approx(np.log(x) + lg, abs=1e-9)
    # the eutectic satisfies both SLE conditions
    assert np.isfinite(pd_.T_e) and 0 < pd_.x_e < 0.05 and pd_.T_e < 273.15
    for which, c in ((0, para), (1, water)):
        assert sol.ln_a(pd_.x_e, pd_.T_e, which) == pytest.approx(ideal_rhs(pd_.T_e, c.dh_fus, c.t_fus), abs=2e-3)
    # ideal curve: closed form 1/T = 1/Tm - R ln x / ΔH
    assert _t_from_lna(np.log(0.5), para) == pytest.approx(1 / (1 / 443.6 - 8.314462618 * np.log(0.5) / 27100.0))
    assert pd_.liquid_window(260.0) is None


def test_library_folder_clash(tmp_path):
    lib = Library(Settings(home=tmp_path, orca_home=None))
    lib.import_orcacosmo("ethanol", [FIX / "ethanol.orcacosmo"], smiles="CCO")
    assert lib.clash("Ethanol") == "ethanol" and lib.clash("ethanol") is None
    with pytest.raises(ValueError, match="storage folder"):
        lib.import_orcacosmo("Ethanol", [FIX / "water.orcacosmo"])
    assert lib.paths("ethanol")[0].read_bytes() == (FIX / "ethanol.orcacosmo").read_bytes()
    # removing a species never deletes files another entry still uses
    lib.put(lib.get("ethanol").__class__(**{**lib.get("ethanol").__dict__, "name": "ETHANOL"}))
    lib.remove("ETHANOL")
    assert lib.paths("ethanol")[0].exists()


def test_multi_conformer_species_uses_lowest_energy(engine):
    lib = engine.library
    lib.import_orcacosmo("water2", [FIX / "ethanol.orcacosmo", FIX / "water.orcacosmo"], smiles="O")
    e = {p.name: descriptors(p)["energy_total_kJmol"] for p in lib.paths("water2")}
    best = engine.conformer("water2")
    assert e[best.name] == min(e.values())
    assert np.isfinite(engine.ln_gamma(["water2", "ethanol"], np.array([[0.5, 0.5]]), 298.15)).all()
    lib.remove("water2")


def test_user_solvents(tmp_path):
    from osnades.thermo import preset_systems

    lib = Library(Settings(home=tmp_path, orca_home=None))
    lib.import_orcacosmo("my alcohol", [FIX / "ethanol.orcacosmo"], smiles="CCO")
    lib.import_orcacosmo("water", [FIX / "water.orcacosmo"], smiles="O")
    lib.add_solvent("my alcohol")
    lib.add_solvent("EtOH:H2O 3:1", {"my alcohol": 3, "water": 1})
    assert lib.is_solvent("my alcohol")
    names = [s.name for s in preset_systems(lib)]
    assert "my alcohol" in names and "EtOH:H2O 3:1" in names and "Water" in names
    assert SolventSystem.parse("EtOH:H2O 3:1", lib).components == {"my alcohol": 3.0, "water": 1.0}
    with pytest.raises(ValueError, match="preset"):
        lib.add_solvent("water")  # clashes with the preset "Water"
    with pytest.raises(ValueError, match="exists"):
        lib.add_solvent("My Alcohol")
    lib.remove("my alcohol")  # deleting a species drops its pure-solvent entry
    assert not lib.is_solvent("my alcohol") and "EtOH:H2O 3:1" in lib.user_solvents()


def test_fusion_page_parsing_and_choice():
    from osnades import fusion

    def row(prop, val, unit, src, ref="", calc=False):
        mark = '<span class="marker" title="Calculated Property"></span>' if calc else ""
        return (f'<tr><td><span>{prop}</span></td><td class="r">{val}</td><td class="c">{unit}</td>'
                f'<td class="r"><a class="src-name" href="#ref-{ref}" title="{src}">{src}</a> {mark}</td></tr>')

    page = "<h1>Chemical Properties of Testol (CAS 1-2-3)</h1><table>" + "".join([
        row("Δ<sub>fus</sub>H°", "29.24", "kJ/mol", "Joback", "chemeo-joback", calc=True),
        row("Δ<sub>fus</sub>H°", "[27.42; 34.30]", "kJ/mol", ""),
        row("Δ<sub>fus</sub>H°", "27.42", "kJ/mol", "Paper A", "doi_10.1016_j.jct.2013.02.013"),
        row("T<sub>fus</sub>", "492.35", "K", "Paper A", "doi_10.1016_j.jct.2013.02.013"),
        row("T<sub>fus</sub>", "489.00 ± 5.00", "K", "NIST", "nist-webbook"),
        row("T<sub>fus</sub>", "213.15", "K", "Some dataset", "ds"),
    ]) + "</table>"
    recs = fusion.parse_page(page, "https://x/cid/1", match="exact")
    assert {r.source for r in recs} == {"Paper A", "NIST", "Some dataset"}  # calculated + range rows dropped
    a = next(r for r in recs if r.source == "Paper A")
    assert a.dh_fus == 27420 and a.t_fus == 492.35 and a.doi == "10.1016/j.jct.2013.02.013"
    assert a.compound == "Testol (CAS 1-2-3)"
    pick = fusion.best(recs)
    assert pick.source == "Paper A" and not pick.estimated  # 2 sources near 490 K beat the lone 213 K value
    est = fusion.best([r for r in recs if r.source == "NIST"])
    assert est.estimated and abs(est.dh_fus - 56.5 * 489.0) < 1e-6 and "Walden" in est.note
    # an entry without stereo information is the same compound: no "verify" warning
    un = fusion.parse_page(page, "https://x/cid/1", match="unspecified")
    assert "verify" not in fusion.best(un).note and "no stereo" in fusion.best(un).note
    iso = fusion.parse_page(page, "https://x/cid/1", match="isomer")
    assert "verify" in fusion.best(iso).note
