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
