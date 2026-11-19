"""Structure handling: SMILES -> 3D conformers (RDKit), XYZ parsing, fragment splitting."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")
except ImportError:  # pragma: no cover
    Chem = None

ATOMIC_MASS = {
    "H": 1.008, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "K": 39.098, "Ca": 40.078, "Zn": 65.38, "Br": 79.904, "I": 126.904,
}


@dataclass
class Geometry:
    symbols: list[str]
    coords: np.ndarray  # (n, 3) Angstrom
    comment: str = ""

    def to_xyz(self) -> str:
        lines = [str(len(self.symbols)), self.comment.replace("\n", " ")]
        for s, (x, y, z) in zip(self.symbols, self.coords):
            lines.append(f"{s:<3} {x:14.8f} {y:14.8f} {z:14.8f}")
        return "\n".join(lines) + "\n"

    @property
    def formula(self) -> str:
        from collections import Counter

        c = Counter(self.symbols)
        order = (["C", "H"] if "C" in c else []) + sorted(k for k in c if k not in ("C", "H") or "C" not in c)
        return "".join(f"{el}{c[el] if c[el] > 1 else ''}" for el in order)

    @property
    def molar_mass(self) -> float:
        return float(sum(ATOMIC_MASS.get(s, 0.0) for s in self.symbols))


def parse_xyz(text: str) -> Geometry:
    """Parse an XYZ block. Accepts files with or without the two header lines."""
    lines = [ln for ln in text.strip("\n").splitlines()]
    if not lines:
        raise ValueError("Empty XYZ input")
    comment = ""
    first = lines[0].split()
    if len(first) == 1 and first[0].isdigit():
        n = int(first[0])
        comment = lines[1] if len(lines) > 1 else ""
        body = lines[2 : 2 + n]
    else:
        body = lines
    symbols, coords = [], []
    for ln in body:
        parts = ln.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        sym = sym[0].upper() + sym[1:].lower()
        symbols.append(sym)
        coords.append([float(v) for v in parts[1:4]])
    if not symbols:
        raise ValueError("No atoms found in XYZ input")
    return Geometry(symbols, np.array(coords), comment)


def _require_rdkit():
    if Chem is None:
        raise ImportError("RDKit is required for SMILES input (pip install rdkit)")


def canonical_smiles(smiles: str) -> str:
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    return Chem.MolToSmiles(mol)


def split_fragments(smiles: str) -> list[str]:
    """Split a multi-fragment SMILES (e.g. a salt) into canonical fragment SMILES."""
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    return [Chem.MolToSmiles(f) for f in Chem.GetMolFrags(mol, asMols=True)]


def smiles_info(smiles: str) -> dict:
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    nrad = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    return {
        "smiles": Chem.MolToSmiles(mol),
        "charge": Chem.GetFormalCharge(mol),
        "multiplicity": nrad + 1,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molar_mass": Descriptors.MolWt(mol),
        "n_fragments": len(Chem.GetMolFrags(mol)),
    }


def smiles_to_conformers(
    smiles: str, n_conformers: int = 1, seed: int = 0xC0FFEE, rmsd_prune: float = 0.5,
    energy_window: float = 5.0,
) -> list[Geometry]:
    """Generate up to ``n_conformers`` low-energy MMFF conformers (sorted by energy).

    ``energy_window`` (kcal/mol) drops conformers far above the minimum.
    """
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    n_atoms = mol.GetNumAtoms()
    if n_atoms == 1:
        atom = mol.GetAtomWithIdx(0)
        return [Geometry([atom.GetSymbol()], np.zeros((1, 3)), f"{smiles} single atom")]

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = rmsd_prune
    n_embed = max(10, 5 * n_conformers) if n_conformers > 1 else 1
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_embed, params=params))
    if not cids:
        params.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_embed, params=params))
    if not cids:
        raise RuntimeError(f"RDKit could not embed {smiles!r}")

    energies: dict[int, float] = {}
    if AllChem.MMFFHasAllMoleculeParams(mol):
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    else:
        res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=2000)
    for cid, (_, e) in zip(cids, res):
        energies[cid] = e

    ranked = sorted(cids, key=lambda c: energies[c])
    emin = energies[ranked[0]]
    ranked = [c for c in ranked if energies[c] - emin <= energy_window][:n_conformers]

    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    out = []
    for cid in ranked:
        pos = mol.GetConformer(cid).GetPositions()
        out.append(Geometry(symbols, np.array(pos), f"{smiles} conf {cid} E_FF={energies[cid]:.3f} kcal/mol"))
    return out
