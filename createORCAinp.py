import os

def main():    # Example usage
    xyz_file = 'molecule.xyz'  # Input .xyz file
    output_inp_file = 'molecule.inp'  # Output .inp file
    charge = 0  # Overall charge of the molecule
    multiplicity = 1  # Spin multiplicity of the molecule

    createORCAinp(xyz_file, output_inp_file, charge, multiplicity)

def createORCAinp(xyz_file, output_inp_file, charge=0, multiplicity=1):
    """
    Reads an .xyz file and generates a basic ORCA input file (.inp).

    Args:
        xyz_file (str): Path to the input .xyz file.
        output_inp_file (str): Path for the output .inp file.
        charge (int): The overall charge of the molecule.
        multiplicity (int): The spin multiplicity of the molecule.
    """
    try:
        with open(xyz_file, 'r') as f:
            lines = f.readlines()

        # The first two lines of an .xyz file are not coordinates
        # Line 1: Number of atoms
        # Line 2: Comment line
        coords_lines = lines[2:]
        
        # Join the coordinate lines to be inserted into the ORCA template
        coords_block = "".join(coords_lines)
        
        # Define the ORCA input file template
        # You can customize this to change the functional, basis set, etc.
        orca_template = f"""
# Basic Single-Point Energy Calculation
# Functional: B3LYP
# Basis Set: def2-TZVP
! B3LYP def2-TZVP SP

# Specify computational resources
%maxcore 2000
%pal
  nprocs 14
end

# Molecular geometry block
*xyz {charge} {multiplicity}
{coords_block}*
"""
        
        # Write the complete input file
        with open(output_inp_file, 'w') as out_f:
            out_f.write(orca_template)

        print(f"ORCA input file '{output_inp_file}' created successfully.")
    
    except FileNotFoundError:
        print(f"Error: The file '{xyz_file}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()