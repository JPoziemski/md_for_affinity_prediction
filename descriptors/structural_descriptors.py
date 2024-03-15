
import pandas as pd
import MDAnalysis as mda
from rdkit.Chem import Descriptors as rdkit_desc
from collections import Counter
from functools import partial
from MDAnalysis.core.universe import Universe



class StructuralDescriptor:

    def __init__(self, mol_obj, universe, descriptors):
        self.universe = universe
        self.mol_obj = mol_obj
        print(self.mol_obj)
        self.check_descriptors(descriptors)
        self.descriptors = descriptors
        # print(self.bonds, self.rings, self.num_atoms)

    def calculate(self):
        if self.universe:
            if len(self.universe.trajectory) > 1:
                data = self._calculate_for_universe()
            else:
                self.mol_obj = self.convert_to_rdkit()
                data = {descriptor: self.calculate_descriptor(self.mol_obj, descriptor) for descriptor in
                        self.descriptors}
        else:
            data = {descriptor: self.calculate_descriptor(self.mol_obj, descriptor) for descriptor in
                    self.descriptors}

        self.data = data

    def _calculate_for_universe(self):
        all_data = []
        for _ in self.universe.trajectory:
            mol_obj = self.convert_to_rdkit()
            all_data.append({descriptor: self.calculate_descriptor(mol_obj, descriptor) for descriptor in self.descriptors})
            #print(len(self.mol_obj.atoms), frame_data)
        return all_data

    def get_data(self):
        if self.universe:
            if len(self.universe.trajectory) > 1:
                data = pd.DataFrame(self.data).to_dict('list')
            else:
                data = self.data
        else:
            data = self.data
        return data

    @classmethod
    def get_num_bonds(cls, mol):
        return mol.GetNumBonds()

    @classmethod
    def get_num_atoms(cls, mol):
        return mol.GetNumAtoms()

    def convert_to_rdkit(self):
        if isinstance(self.mol_obj, mda.core.groups.AtomGroup):
            mol_obj_rdkit = self.mol_obj.convert_to("RDKIT",  NoImplicit=False, force=True)
        elif isinstance(self.mol_obj, mda.core.groups.UpdatingAtomGroup):
            mol_obj_rdkit = self.mol_obj.residues.atoms.convert_to("RDKIT", force=True)
        else:
            pass
        return mol_obj_rdkit

    @classmethod
    def get_descriptors_mapper(cls):
        mapper = {
            "arom_rings": rdkit_desc.NumAromaticRings,
            "aliphatic_rings": rdkit_desc.NumAliphaticRings,
            "hba": rdkit_desc.NumHAcceptors,
            "hbd": rdkit_desc.NumHDonors,
            "rot_bonds": rdkit_desc.NumRotatableBonds,
            "single_bonds": partial(cls.get_bonds, bond_type="SINGLE"),
            "double_bonds": partial(cls.get_bonds, bond_type="DOUBLE"),
            "aromatic_bonds": partial(cls.get_bonds, bond_type="AROMATIC")

        }
        return mapper

    @classmethod
    def calculate_descriptor(cls, mol, descriptor):
        mapper = cls.get_descriptors_mapper()
        descriptor_value = mapper[descriptor](mol)
        # print(mapper[descriptor](m))
        # print(m, mapper)
        return descriptor_value

    @classmethod
    def _get_bond_info(cls, mol):
        bond_info = [str(bond.GetBondType()) for bond in mol.GetBonds()]
        return bond_info

    @classmethod
    def get_bonds(cls, mol, bond_type):
        bond_info = cls._get_bond_info(mol)
        bond_count = Counter(bond_info)
        bond_data = bond_count.get(bond_type, 0)
        return bond_data

    def check_descriptors(self, descriptors):
        mapper = self.get_descriptors_mapper()
        unknown_descriptors = set(descriptors) - set(mapper.keys())
        if unknown_descriptors:
            raise NameError(f"Unknown descriptors: {', '.join(unknown_descriptors)}")
