import os
import argparse
from utils import prepare_trajectory, get_ligand_from_merged_universe, get_pocket
from descriptors import FullLigandDescriptor, FullPocketDescriptor, InteractionDescriptor
import pandas as pd


def compute_trajectory_descriptors(topology_path, trajectory_path, name):
    universe = prepare_trajectory(topology_path, trajectory_path)
    ligand_universe = get_ligand_from_merged_universe(universe)
    ligand_code = ligand_universe.residues.resnames[0]

    pocket_universe = get_pocket(universe, ligand_code, 3.5)

    ligand_desc = FullLigandDescriptor(universe, ligand_universe)
    ligand_desc.calculate()
    ligand_data = ligand_desc.get_data()

    pocket_desc = FullPocketDescriptor(universe, pocket_universe, ligand_code)
    pocket_desc.calculate()
    pocket_data = pocket_desc.get_data()

    pocket_to_interaction = universe.select_atoms('protein').atoms
    interaction_desc = InteractionDescriptor(universe, ligand_universe, pocket_to_interaction)
    interaction_desc.calculate()
    interaction_data = interaction_desc.get_data()

    complex_data = pd.DataFrame({**interaction_data, **pocket_data, **ligand_data})
    complex_data['Name'] = name
    complex_data['Frame'] = list(range(len(complex_data)))
    complex_data.insert(0, 'Name', complex_data.pop('Name'))
    complex_data.insert(1, 'Frame', complex_data.pop('Frame'))
    return complex_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
    Script for generating descriptors based on molecular dynamics trajectories. 
    The script as the input parameters take the path to the data folder with the molecular dynamics trajectories. 
    The structure of the input folder should looks:
    md_dir:
    ├── trajectory_1_dir/.
    │ ├── trajectory.xtc
    │ └── topology.gro
    ├── trajectory_2_dir/
    │ ├── trajectory.xtc
    │ └── topology.gro
    ....
    The script generates one resulting csv file in the following format:
    +------------+-------+--------------+--------------+----+
    | Name       | Frame | descriptor 1 | descriptor 2 | ...|         
    +------------+-------+--------------+--------------+----+
    | trajectory1| 1     | value        | value        | ...|
    | trajectory1| 2     | value        | value        | ...|
    | ...        | ...   | ...          | ...          | ...| 
    | trajectory2| 2     | value        | value        | ...|
    | ...        | ...   | ...          | ...          | ...|
    +------------+-------+--------------+--------------+----+
    ''')
    parser.add_argument("input_directory", help="path to the input folder")
    parser.add_argument("output_path", help="path to csv file with script results")

    args = parser.parse_args()

    if not os.path.exists(args.input_directory):
        raise FileNotFoundError(f"Directory {args.input_directory} is not exists")

    trajectory_data = []

    for sub_directory in os.listdir(args.input_directory)[:2]:
        directory_path = os.path.join(args.input_directory, sub_directory)
        topology_path = None
        trajectory_path = None
        for file in os.listdir(directory_path):
            print(file)
            ext = os.path.splitext(file)[1]
            if ext == '.xtc':
                trajectory_path = os.path.join(directory_path, file)
            elif ext == '.gro':
                topology_path = os.path.join(directory_path, file)
            else:
                pass

        if not (topology_path and trajectory_path):
            print(f"{sub_directory} does not contain all necessary files. Skipping...")
            continue
        complex_data = compute_trajectory_descriptors(topology_path, trajectory_path, sub_directory)
        trajectory_data.append(complex_data)
    traj_df = pd.concat(trajectory_data)
    traj_df.to_csv(args.output_path, index=False)