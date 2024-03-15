
from MDAnalysis.analysis.rms import rmsd

class MotionDescriptor:

    def __init__(self, universe, ligand_code, descriptors=("RMSD_pocket", "RMSD_ligand", "RMSD_ca")):
        self.universe = universe
        self.ligand_code = ligand_code
        print(f"around 3.5 resname {ligand_code}")
        self.pocket = universe.select_atoms(f"around 3.5 resname {ligand_code}").residues.atoms
        self.pocket = self.pocket.select_atoms(f"not type H")
        self.ligand = self.universe.select_atoms(f"resname {ligand_code} and not type H")
        self.descriptors = descriptors

    def _get_ligand_position_RMSD_to_previous_frame(self):

        ref_pos = self.ligand.positions
        values = []
        for _ in self.universe.trajectory:
            value = rmsd(self.ligand.positions, ref_pos)
            values.append(value)
            ref_pos = self.ligand.positions

        return values

    def _get_pocket_position_RMSD_to_previous_frame(self):

        ref_pos = self.pocket.positions
        values = []
        for _ in self.universe.trajectory:
            value = rmsd(self.pocket.positions, ref_pos)
            values.append(value)
            ref_pos = self.pocket.positions
        return values

    def _get_ca_position_RMSD_to_previous_frame(self):
        pocket_ca = self.pocket.select_atoms(f"protein and name CA")

        ref_pos = pocket_ca.positions
        values = []
        for _ in self.universe.trajectory:
            value = rmsd(pocket_ca.positions, ref_pos)
            values.append(value)
            ref_pos = pocket_ca.positions
        return values

    def get_descriptors_mapper(self):
        mapper = {
            "RMSD_pocket": self._get_pocket_position_RMSD_to_previous_frame,
            "RMSD_ligand": self._get_ligand_position_RMSD_to_previous_frame,
            "RMSD_ca": self._get_ca_position_RMSD_to_previous_frame,
        }
        return mapper

    def calculate(self):
        mapper = self.get_descriptors_mapper()
        self.data = {descriptor: mapper[descriptor]() for descriptor in self.descriptors}

    def get_data(self):
        return self.data