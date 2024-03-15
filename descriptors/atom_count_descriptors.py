
import MDAnalysis as mda
import pandas as pd
from rdkit.Chem import AddHs
from collections import Counter
from MDAnalysis.core.universe import Universe
from .descriptor import Descriptor

class AtomCountDescriptor(Descriptor):
    ATOMS_CONVERTER = {
        "F": "Halogen",
        "Cl": "Halogen",
        "I": "Halogen",
        "Br": "Halogen",
        "Na": "Metal",
        "K": "Metal",
        "Ca": "Metal",
        "Li": "Metal",
        "Mg": "Metal",
        "Fe": "Metal",
        "Cu": "Metal",
        "Zn": "Metal"
    }

    def __init__(self, mol_obj, universe, elements: list):
        self.universe = universe
        self.elements = set(elements)
        self.mol_obj = mol_obj

    def __len__(self):
        return len(self.mol_elements)

    def calculate(self):
        if self.universe:
            if len(self.universe.trajectory) > 1:
                data = self._calculate_for_universe()
            else:
                self.mol_obj = self.convert_to_rdkit()
                data = self.calculate_atoms(self.mol_obj)
        else:
            self.mol_obj = self.convert_to_rdkit()
            data = self.calculate_atoms(self.mol_obj)

        self.data = data

    def _calculate_for_universe(self):
        all_data = []
        for _ in self.universe.trajectory:
            mol_obj = self.convert_to_rdkit()
            all_data.append(self.calculate_atoms(mol_obj))
            #print(len(self.mol_obj.atoms), frame_data)
        return all_data

    def convert_to_rdkit(self):
        if isinstance(self.mol_obj, mda.core.groups.AtomGroup):
            mol_obj_rdkit = self.mol_obj.convert_to("RDKIT", NoImplicit=False, force=True)
        elif isinstance(self.mol_obj, mda.core.groups.UpdatingAtomGroup):
            mol_obj_rdkit = self.mol_obj.residues.atoms.convert_to("RDKIT", force=True)
        else:
            pass
        mol_obj_rdkit = AddHs(mol_obj_rdkit)
        return mol_obj_rdkit

    def calculate_atoms(self, mol_obj):
        mol_elements = list(at.GetSymbol() for at in mol_obj.GetAtoms())
        mol_elements = self._convert_atoms(mol_elements)
        mol_elements = [element for element in mol_elements if element in self.elements]
        data = dict(Counter(mol_elements))
        return data

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
    def _convert_atoms(cls, mol_elements):
        converted_mol_elements = list(
            cls.ATOMS_CONVERTER.get(at_symbol) if cls.ATOMS_CONVERTER.get(at_symbol) else at_symbol for at_symbol in
            mol_elements)
        return converted_mol_elements