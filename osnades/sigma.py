"""σ-profiles and molecular descriptors from .orcacosmo files (openCOSMO-RS parser)."""
from __future__ import annotations

from pathlib import Path

import numpy as np


def sigma_profile(path: Path | str) -> tuple[np.ndarray, np.ndarray]:
    from opencosmorspy.input_parsers import SigmaProfileParser

    spp = SigmaProfileParser(str(path))
    sigmas, areas = spp.cluster_and_create_sigma_profile()
    return np.asarray(sigmas), np.asarray(areas)


def descriptors(path: Path | str) -> dict:
    from opencosmorspy.input_parsers import SigmaProfileParser

    spp = SigmaProfileParser(str(path))
    spp.calculate_sigma_moments()
    m = np.asarray(spp["sigma_moments"])
    return {
        "method": spp["method"],
        "charge": float(np.sum(spp["seg_charge"])),
        "area_A2": float(spp["area"]),
        "volume_A3": float(spp["volume"]),
        "energy_dielectric_kJmol": float(spp["energy_dielectric"]),
        "energy_total_kJmol": float(spp["energy_tot"]),
        "sigma_moments": m.tolist(),
        "hb_donor_moments": np.asarray(spp["sigma_hydrogen_bond_donor_moments"])[2:5].tolist(),
        "hb_acceptor_moments": np.asarray(spp["sigma_hydrogen_bond_acceptor_moments"])[2:5].tolist(),
        "n_segments": int(len(spp["seg_area"])),
    }
