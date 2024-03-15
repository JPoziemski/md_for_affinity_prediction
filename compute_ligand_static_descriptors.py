import os
import argparse
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from utils import create_universe_for_crystal_structure
from rdkit.Chem import SmilesMolSupplier
from descriptors import LigandStaticDescriptor


def compute_descriptors_for_molecules(molecules):
    mol_static_descriptors = []
    for mol in molecules:
        mol_name = mol.GetProp("_Name")
        mol_desc = LigandStaticDescriptor(mol)
        mol_desc.calculate()
        mol_data = mol_desc.get_data()
        mol_data['Name'] = mol_name
        mol_static_descriptors.append(mol_data)
    mol_static_descriptors_df = pd.DataFrame(mol_static_descriptors)
    return mol_static_descriptors_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
        Script to generate static ligand descriptors (e.g. mass, logp ). 
        The script takes as input parameter smi file.
        Structure of typical smi file: 
        ----------------------------------------
        CCO	mol1_id
        CNC	mol2_id
        CC(=O)Oc1ccccc1C(=O)O mol3_id
        ....
        ----------------------------------------
        Make sure that mol_ids are unique.
        The script generates one resulting csv file in the following format:
        +---------+--------------+--------------+----+
        | Name    | descriptor 1 | descriptor 2 | ...|         
        +---------+--------------+--------------+----+
        | mol1_id | value        | value        | ...|
        | mol2_id | value        | value        | ...|
        | mol3_id | value        | value        | ...|
        | ...     | ...          | ...          | ...| 
        +-----------+--------------+--------------+----+
        ''')
    parser.add_argument("input_file", help="path to the input smi file")
    parser.add_argument("output_path", help="path to csv file with script results")

    args = parser.parse_args()

    mols = [mol for mol in SmilesMolSupplier(args.input_file,titleLine=False) if mol]
    mol_static_descriptors_df = compute_descriptors_for_molecules(mols)
    mol_static_descriptors_df = mol_static_descriptors_df.set_index('Name')
    mol_static_descriptors_df = mol_static_descriptors_df.fillna(0)
    mol_static_descriptors_df.to_csv(args.output_path)