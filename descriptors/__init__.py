from .atom_count_descriptors import AtomCountDescriptor
from .contact_descriptor import ContactDescriptor
from .interaction_descriptors import InteractionDescriptor
from .ligand_3d_descriptors import Ligand3DPropertyDescriptor
from .ligand_descriptors import LigandProperty
from .motion_desriptors import MotionDescriptor
from .residue_type_descriptor import ResidueTypeDescriptors
from .shape_descriptor import ShapeDescriptor
from .structural_descriptors import StructuralDescriptor
from .ecfp import ECFP


class LigandStaticDescriptor:
    def __init__(self, ligand):
        self.ligand = ligand
        self.descriptors = {
            "LigandProperty": LigandProperty(self.ligand, ["mol_weigth", "logp", "hba", "hbd", "tpsa", "mr"]),
            "StructuralDescriptor": StructuralDescriptor(self.ligand, None,
                                                         ["arom_rings", "aliphatic_rings", "rot_bonds", "single_bonds",
                                                          "double_bonds", "aromatic_bonds"]),
            "AtomCountDescriptor": AtomCountDescriptor(self.ligand, None,["H", "C", "N", "O", "S", "P", "Halogen"]),
            "ECFP4": ECFP(self.ligand)
        }

    def get_data(self):
        data = {}
        for _, descriptor_object in self.descriptors.items():
            data = {**data, **descriptor_object.get_data()}
        #print(data)
        return data

    def calculate(self):
        for _, descriptor_object in self.descriptors.items():
            descriptor_object.calculate()


class FullLigandDescriptor:
    def __init__(self, universe, ligand, ligand_for_3d=None):
        self.ligand = ligand
        self.ligand_for_3d = ligand_for_3d
        self.universe = universe
        self.descriptors = {
            "Shape": ShapeDescriptor(universe, self.ligand)
        }
        if ligand_for_3d:
             self.descriptors["Ligand3D"] = Ligand3DPropertyDescriptor(self.ligand_for_3d)
        else:
            self.descriptors["Ligand3D"] = Ligand3DPropertyDescriptor(self.ligand, universe)

    def calculate(self):
        for _, descriptor_object in self.descriptors.items():
            print(descriptor_object)
            descriptor_object.calculate()

    def get_data(self):
        data = {}
        for _, descriptor_object in self.descriptors.items():
            data = {**data, **descriptor_object.get_data()}
            #data = {f"ligand_{key}"}
        #print(data)
        return data


class FullPocketDescriptor:
    def __init__(self, universe, pocket, ligand_code, pocket_from_file=None):
        self.ligand = universe.select_atoms(f"resname {ligand_code} and not type H")
        self.universe = universe
        self.pocket = pocket
        pocket_5 = universe.select_atoms(f"around 5 resname {ligand_code}")
        #self.pocket_structural = self.pocket.atoms.select_atoms('not type H')

        self.descriptors = {
            "StructuralDescriptor": StructuralDescriptor(self.pocket, universe,
                                                         ["arom_rings", "aliphatic_rings", "rot_bonds", "single_bonds",
                                                          "double_bonds", "aromatic_bonds", "hba", "hbd"]),
            "AtomCountDescriptor": AtomCountDescriptor(self.pocket, universe, ["H", "C", "N", "O", "S"]),
            "ContactDescriptor": ContactDescriptor(self.universe, ligand_code, ["contact_pocket_all", "contact_ligand_pocket_all"]),
            "Shape": ShapeDescriptor(self.universe, self.pocket),
            "Shape_5": ShapeDescriptor(self.universe, pocket_5, suffix="_5"),
            "ShapeLigand": ShapeDescriptor(self.universe, self.ligand, suffix="_ligand"),
            "ResidueTypeDescriptors": ResidueTypeDescriptors(self.universe, self.pocket),
        }
        if pocket_from_file:
            self.descriptors["StructuralDescriptor"] = StructuralDescriptor(self.pocket,
                                                         ["arom_rings", "aliphatic_rings", "rot_bonds", "single_bonds",
                                                          "double_bonds", "aromatic_bonds"])
            self.descriptors["AtomCountDescriptor"] = AtomCountDescriptor(self.pocket, ["H", "C", "N", "O", "S"])

        if len(universe.trajectory) > 1:
            self.descriptors["MotionDescriptors"] = MotionDescriptor(self.universe, ligand_code)

    def calculate(self):
        for _, descriptor_object in self.descriptors.items():
            descriptor_object.calculate()

    def get_data(self):
        data = {}
        for _, descriptor_object in self.descriptors.items():

            data = {**data, **descriptor_object.get_data()}
        # print(data)
        return data

