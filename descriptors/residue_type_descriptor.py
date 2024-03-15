from collections import defaultdict

import pandas as pd


class ResidueTypeDescriptors:

    def __init__(self, universe,  mol_obj, descriptors=("aromatic", "aliphatic", "charged", "hydrophobic", "polar")):
        self.pocket = mol_obj
        self.universe = universe
        self.check_descriptors(descriptors)
        self.descriptors = descriptors
        used_descriptors_dict = {desc: value for desc, value in self.get_aa_groups().items() if desc in descriptors}
        self._inverted_dict = self.invert_aa_group_dict(used_descriptors_dict)
        self.residues = self.get_residues()
        # print(self.residues)

    def __len__(self):
        return len(self.residues)

    def get_data(self):
        return self.data

    @classmethod
    def get_aa_groups(cls):
        aa_groups = {
            "aromatic": ["PHE", "TRP", "TYR"],
            "aliphatic": ["ALA", "ILE", "LEU", "PRO", "VAL"],
            "charged": ["ARG", "LYS", "ASP", "GLU"],
            "hydrophobic": ["ALA", "ILE", "LEU", "MET", "PHE", "VAL", "PRO", "GLY"],
            "polar": ["GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS"]

        }
        return aa_groups

    @classmethod
    def invert_aa_group_dict(cls, descriptors_dict):
        inverted_dict = defaultdict(list)
        for key, values in descriptors_dict.items():
            for aa_code in values:
                inverted_dict[aa_code].append(key)
        return inverted_dict

    @classmethod
    def list_avaliable_descriptors(cls):
        return list(cls.get_aa_groups().keys())

    def calculate(self):
        if len(self.universe.trajectory) >1:
            data = self._calculate_for_universe()
        else:
            data = self.calculate_descriptors()

        self.data = data

    def calculate_descriptors(self):
        frame_dict = defaultdict(int)
        for resname in self.residues:
            resname_descriptors = self._inverted_dict.get(resname, [])
            for desc_name in resname_descriptors:
                frame_dict[desc_name] += 1
        return frame_dict

    def _calculate_for_universe(self):
        descriptors_dict =[]
        for _ in self.universe.trajectory:
            frame_dict = defaultdict(int)
            for resname in self.residues:
                resname_descriptors = self._inverted_dict.get(resname, [])
                for desc_name in resname_descriptors:
                    frame_dict[desc_name] += 1

            descriptors_dict.append(dict(frame_dict))
        data = pd.DataFrame(descriptors_dict).to_dict('list')
        #print(data)
        return data

    def get_residues(self):
        residues = list(self.pocket.residues.resnames)
        return residues


    @classmethod
    def check_descriptors(cls, descriptors):
        unkonwn = set(descriptors) - set(cls.list_avaliable_descriptors())
        if unkonwn:
            raise NameError(f"Unknown descriptors: {', '.join(unkonwn)}")