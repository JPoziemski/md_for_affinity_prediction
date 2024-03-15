
import prolif as plf
from collections import Counter, defaultdict
import pandas as pd



class InteractionDescriptor:
    def __init__(self, universe, ligand, pocket, interactions=('Anionic',
                                                               'CationPi',
                                                               'Cationic',
                                                               'EdgeToFace',
                                                               'FaceToFace',
                                                               'HBAcceptor',
                                                               'HBDonor',
                                                               'Hydrophobic',
                                                               'PiCation',
                                                               'PiStacking',
                                                               'VdWContact',
                                                               'XBAcceptor',
                                                               'XBDonor')):
        self.universe = universe
        self.pocket = pocket
        # print(self.pocket)
        self.ligand = ligand
        self.interactions = interactions
        self.fp = plf.Fingerprint(interactions, count=True)

    def calculate(self):
        # print(self.universe.trajectory.n_frames)
        if self.universe.trajectory.n_frames > 1:
            self.fp.run(self.universe.trajectory, self.ligand, self.pocket)
            self.data = self.fp.to_dataframe()
        else:
            ligand_mol = plf.Molecule.from_rdkit(self.ligand)
            pocket_mol = plf.Molecule.from_rdkit(self.pocket)
            ifp = self.fp.generate(ligand_mol, pocket_mol, metadata=True)
            self.data = plf.to_dataframe({0: ifp}, self.fp.interactions, count=True)

    def get_data(self):
        interactions_column_dict = self.__group_interactions(self.data)
        self.n_interactions = self.data.sum(axis=1).values
        self.data = self.__convert_to_interaction_dict(self.data, interactions_column_dict)
        if len(self.universe.trajectory) == 1:
            self.data = {key: value[0] for key, value in self.data.items()}
        return self.data

    def __group_interactions(self, data):

        interactions_column_dict = defaultdict(list)
        for col_index in data.columns:
            interactions_column_dict[col_index[2]].append(col_index)

        return interactions_column_dict

    def __convert_to_interaction_dict(self, data, interactions_column_dict):
        interactions_count_dict = {}
        # print("data", data.sum(axis=1))
        for interaction_type, col_indexes in interactions_column_dict.items():
            interactions_count_dict[interaction_type] = data[col_indexes].apply(sum, axis=1).values
        return interactions_count_dict