import os
import argparse
import warnings
from rdkit.Chem import MolFromMol2File, SanitizeMol, MolFromPDBFile, AddHs
warnings.filterwarnings("ignore")
import pandas as pd
from utils import create_universe_for_crystal_structure, get_pocket
from descriptors import FullLigandDescriptor, FullPocketDescriptor, InteractionDescriptor


def compute_descriptors_for_complex(protein_path, ligand_path, name):
    universe, ligand_universe = create_universe_for_crystal_structure(protein_path, ligand_path)
    ligand_code = ligand_universe.residues.resnames[0]

    ligand_mol = ligand_universe.atoms.convert_to("RDKIT", NoImplicit=False, force=True)
    AddHs(ligand_mol)

    pocket_universe = get_pocket(universe, ligand_code, 3.5, updating=False)

    ligand_desc = FullLigandDescriptor(universe, ligand_universe.atoms, ligand_mol)
    ligand_desc.calculate()
    ligand_data = ligand_desc.get_data()


    pocket_desc = FullPocketDescriptor(universe, pocket_universe, ligand_code)
    pocket_desc.calculate()
    pocket_data = pocket_desc.get_data()

    protein_mol = MolFromPDBFile(protein_path)

    interaction_desc = InteractionDescriptor(universe, ligand_mol, protein_mol)
    interaction_desc.calculate()
    interaction_data = interaction_desc.get_data()


    complex_data = {**interaction_data, **pocket_data, **ligand_data}
    complex_data['Name'] = name

    return complex_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
        Script for generating descriptors from protein ligand complex structures. 
        The script as the input parameters take the path to the data folder with pdb files. 
        The structure of the input folder should look like:
        input_dir:
        ├── structure1_dir/.
        │ ├── protein.pdb
        │ └── ligand.mol2
        ├── structure2_dir/
        │ ├── protein.pdb
        │ └── ligand.mol2
        ....
        Such a structure is analogous to the one used in the pdbbind dataset.
        The script generates one resulting csv file in the following format:
        +-----------+--------------+--------------+----+
        | Name      | descriptor 1 | descriptor 2 | ...|         
        +-----------+--------------+--------------+----+
        | structure1| value        | value        | ...|
        | structure2| value        | value        | ...|
        | structure3| value        | value        | ...|
        | ...       | ...          | ...          | ...| 
        +-----------+--------------+--------------+----+
        ''')
    parser.add_argument("input_directory", help="path to the input folder")
    parser.add_argument("output_path", help="path to csv file with script results")

    args = parser.parse_args()

    if not os.path.exists(args.input_directory):
        raise FileNotFoundError(f"Directory {args.input_directory} is not exists")

    complexes_data = []
    for sub_directory in os.listdir(args.input_directory)[1:]:
        directory_path = os.path.join(args.input_directory, sub_directory)
        protein_path = None
        ligand_path = None
        for file in os.listdir(directory_path):

            ext = os.path.splitext(file)[1]
            if ext == '.pdb':
                protein_path = os.path.join(directory_path, file)
            elif ext == '.mol2':
                ligand_path = os.path.join(directory_path, file)
            else:
                pass

        if not (protein_path and ligand_path):
            print(f"{sub_directory} does not contain all necessary files. Skipping...")
            continue
        complex_data = compute_descriptors_for_complex(protein_path, ligand_path, sub_directory)
        complexes_data.append(complex_data)

    complexes_df = pd.DataFrame(complexes_data)
    complexes_df = complexes_df.set_index('Name')
    complexes_df.to_csv(args.output_path)
