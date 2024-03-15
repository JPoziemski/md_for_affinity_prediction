import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_types, guess_bonds


def prepare_trajectory(topology_path, trajectory_path):
    universe = mda.Universe(topology_path, trajectory_path)
    universe = repair_universe(universe)
    return universe


def repair_universe(universe):
    guessed_elements = guess_types(universe.atoms.names)
    universe.add_TopologyAttr('elements', guessed_elements)

    guessed_bonds = guess_bonds(universe.atoms, universe.atoms.positions)
    universe.add_TopologyAttr('bonds', guessed_bonds)
    return universe


def get_pocket(universe, ligand_code, distance=6.0, updating=True):
    command = f"(not resname {ligand_code} HOH) and (around {str(distance)} resname {ligand_code}) and protein"
    if updating:
        pocket = universe.select_atoms(command, updating=updating)
    else:
        pocket = universe.select_atoms(command, updating=updating).residues.atoms
    print(type(pocket))
    return pocket


def get_ligand_from_merged_universe(universe, H=True, updating=True):
    if updating:
        ligand = universe.select_atoms('not protein', updating=True)
    else:
        ligand = universe.select_atoms('not protein').atoms

    if not H:
        ligand = ligand.select_atoms('not type H')
    return ligand


def get_ligand(universe, ligand_code):

    ligand = universe.atoms.select_atoms(f"resname {ligand_code}")
    return ligand


def create_universe_for_crystal_structure(pocket_path, ligand_path):
    protein_universe = mda.Universe(pocket_path)
    protein_universe = repair_universe(protein_universe)
    ligand_universe = mda.Universe(ligand_path)
    #ligand_universe = repair_universe(ligand_universe)
    universe = mda.Merge(protein_universe.atoms, ligand_universe.atoms)

    return universe, ligand_universe
