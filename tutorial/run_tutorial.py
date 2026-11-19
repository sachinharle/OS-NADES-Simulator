"""Script version of the FACCTs OPI openCOSMO-RS notebook (paracetamol solubility)."""
import os, sys, shutil, time
from pathlib import Path
import numpy as np
from scipy.optimize import fsolve
from opencosmorspy.parameterization import openCOSMORS24a
from opencosmorspy.cosmors import COSMORS
from opencosmorspy.input_parsers import SigmaProfileParser
from opi.core import Calculator
from opi.input.structures.structure import Structure

NCORES = int(os.environ.get("NCORES", "4"))
working_dir = Path(__file__).parent / "opencosmors"
shutil.rmtree(working_dir, ignore_errors=True); working_dir.mkdir()

WATER = """3

  O  -0.0007948665470900  0.4014278382603300  0.0000000000000000
  H  -0.7647999815056600  -0.2022173883142000  0.0000000000000000
  H  0.7655948480527500  -0.1992104499461300  0.0000000000000000
"""
PARA = """20
paracetamol xtb-optimised
C      2.88225    1.70760    0.08636
C      1.87583    0.57255    0.08465
H      3.88408    1.29416    0.14629
H      2.79182    2.28846   -0.82886
H      2.71377    2.36342    0.93735
O      2.19537   -0.59390    0.10126
N      0.59013    1.02699    0.05895
C     -0.60624    0.29815    0.04364
H      0.46337    2.02968    0.04757
C     -1.80437    1.01047    0.01777
C     -1.86108   -1.74887    0.03135
H     -1.78549    2.09153    0.01346
C     -3.01629    0.35366   -0.00279
C     -0.64533   -1.09159    0.05178
H      0.27418   -1.65102    0.07346
C     -3.05363   -1.03730    0.00247
H     -3.94534    0.90006   -0.02339
O     -4.26895   -1.65165   -0.02125
H     -1.88115   -2.83055    0.03727
H     -4.14859   -2.60876   -0.01690
"""

def run(basename, xyz):
    (working_dir / f"{basename}.xyz").write_text(xyz)
    calc = Calculator(basename=basename, working_dir=working_dir, version_check=False)
    calc.structure = Structure.from_xyz(working_dir / f"{basename}.xyz")
    calc.input.add_arbitrary_string("!COSMORS(ethanol)")
    calc.input.ncores = NCORES
    calc.write_input()
    t = time.time(); print(f"Running ORCA {basename} ...", end="", flush=True)
    calc.run(); print(f" {time.time()-t:.0f}s")
    out = calc.get_output()
    if not out.terminated_normally():
        raise RuntimeError(f"ORCA failed: {out.get_outfile()}")

run("opencosmors", WATER)
run("opencosmors_paracetamol", PARA)

for name, f in {"Water": "opencosmors.solute.orcacosmo", "Ethanol": "opencosmors.solvent.orcacosmo"}.items():
    spp = SigmaProfileParser(working_dir / f); spp.calculate_sigma_moments()
    print(name, "E_diel", spp["energy_dielectric"], "area", spp["area"], "vol", spp["volume"])

R, TEMP, DH, TF = 8.314, 298.15, 27.1e3, 443.6
def lng(crs, x, T):
    crs.clear_jobs(); crs.add_job(np.array(x), T, refst="pure_component")
    return crs.calculate()["tot"]["lng"][0][0]
def sol(crs, iterative):
    rhs = -DH / R * (1 / TEMP - 1 / TF)
    if not iterative:
        return np.exp(rhs - lng(crs, [0.0, 1.0], TEMP))
    def eq(x):
        x = min(1, max(1e-15, np.asarray(x).item()))
        return np.abs(lng(crs, np.array([x, 1 - x]), TEMP) + np.log(x) - rhs)
    return min(1, max(1e-15, fsolve(eq, 1e-5)[0]))

crs = COSMORS(par=openCOSMORS24a())
for label, f in [("Water", "opencosmors.solute.orcacosmo"), ("Ethanol", "opencosmors.solvent.orcacosmo")]:
    crs.clear_molecules()
    crs.add_molecule([str(working_dir / "opencosmors_paracetamol.solute.orcacosmo")])
    crs.add_molecule([str(working_dir / f)])
    print(f"{label:<10} non-iter {sol(crs, False):.5e}  iter {sol(crs, True):.5e}")
