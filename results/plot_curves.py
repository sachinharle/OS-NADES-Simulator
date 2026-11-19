"""Paracetamol solubility vs temperature in all preset solvent systems (same as `osnades curves paracetamol`)."""
from osnades.cli import main

main(["curves", "paracetamol", "--tmin", "0", "--tmax", "80", "-n", "17", "-o", "results/paracetamol_solubility_vs_T"])
