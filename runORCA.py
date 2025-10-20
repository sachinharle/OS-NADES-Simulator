#!/usr/bin/env python3
# #create input file for orca
import smiles2xyz
import createORCAinp
import pandas as pd
import os
import sys
  
def main(): #Example usage
    runORCA()    
   

def runORCA():
    #parse input csv file with molecule name and SMILES while calling this script
    
    if len(sys.argv) != 2:
        print("Usage: python createORCAinputs.py input.csv")
        sys.exit(1)
    input_csv = sys.argv[1] #input csv file with molecule name and SMILES
    #read input csv file

    df = pd.read_csv(input_csv)
    #iterate over rows of dataframe
    for index, row in df.iterrows():
        name = row['name']
        smiles = row['smiles']
        #convert SMILES to xyz
        xyz = smiles2xyz.smiles_to_xyz(smiles)
        #create ORCA input file
        createORCAinp.create_orca_input(name, xyz) 
        #run ORCA
        os.system(f"orca {name}.inp > {name}.out &") 
        print(f"ORCA calculation for {name} started.") #function to run orca calculations
    #for each molecule, create input file and run orca calculation
if __name__ == "__main__":
    main() 
