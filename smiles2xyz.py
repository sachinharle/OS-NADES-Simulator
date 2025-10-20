from rdkit import Chem
from rdkit.Chem import AllChem

def main():
    # Example usage
    smiles = "CCO"  # Ethanol
    filename = "ethanol.xyz"
    smiles2xyz(smiles, filename)

def smiles2xyz(smiles, filename):
    # 1. Create a 2D molecule object from SMILES
    mol = Chem.MolFromSmiles(smiles)
    
    # 2. Add hydrogens, which are essential for 3D generation
    mol = Chem.AddHs(mol)
    
    # 3. Generate a 3D conformer using the ETKDGv3 algorithm
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    
    # 4. Optimize the molecule's geometry to a low-energy state
    #    This makes the 3D structure more realistic.
    AllChem.UFFOptimizeMolecule(mol)
    
    # 5. Write the 3D coordinates to an XYZ file
    #    The file will be created in your working directory.
    Chem.MolToXYZFile(mol, filename)
    
    print(f"Molecule saved to {filename}")

if __name__ == "__main__":
    main()