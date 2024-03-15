
import MDAnalysis as mda
from collections import defaultdict
from rdkit.Chem import Descriptors3D
from MDAnalysis.core.universe import Universe

class Ligand3DPropertyDescriptor:
    def __init__(self, ligand, universe=None, descriptors=("asphericity", "eccentricity",
                                                           "pbf", "rog", "pmi1", "pmi2", "pmi3")):
        print("adsasads")
        self.ligand = ligand
        self.universe = universe
        print(type(ligand))
        self.descriptors = descriptors
        if isinstance(ligand, mda.core.groups.AtomGroup) and not self.universe:
            raise ValueError(f"if ligand is mda Atom Group you must provide universe")

    @classmethod
    def get_descriptors_mapper(cls):
        mapper = {
            "asphericity": Descriptors3D.Asphericity,
            "eccentricity": Descriptors3D.Eccentricity,
            "pbf": Descriptors3D.rdMolDescriptors.CalcPBF,
            "rog": Descriptors3D.RadiusOfGyration,
            "pmi1": Descriptors3D.PMI1,
            "pmi2": Descriptors3D.PMI2,
            "pmi3": Descriptors3D.PMI3
        }
        return mapper

    def calculate(self):
        if self.universe:
            data = self._calulate_for_universe()
        else:
            data = {descriptor: self.calculate_descriptor(self.ligand, descriptor) for descriptor in self.descriptors}

        self.data = data

    def get_data(self):
        return self.data

    def _calulate_for_universe(self):

        descriptors_dict = defaultdict(list)
        for _ in self.universe.trajectory:
            ligand_rdkit = self.ligand.atoms.convert_to('RDKIT', NoImplicit=False)
            for descriptor in self.descriptors:
                descriptor_value = self.calculate_descriptor(ligand_rdkit, descriptor)
                descriptors_dict[descriptor].append(descriptor_value)
        return dict(descriptors_dict)

    @classmethod
    def calculate_descriptor(cls, mol, descriptor):
        mapper = cls.get_descriptors_mapper()
        descriptor_value = mapper[descriptor](mol)
        # print(mapper[descriptor](m))
        # print(m, mapper)
        return descriptor_value


def compute_crystal_ligand_geometric_descriptor(ligand_data):
    desc_3d = Ligand3DPropertyDescriptor(ligand_data)
    desc_3d.calculate()
    desc_3d_data = desc_3d.get_data()
    return desc_3d_data