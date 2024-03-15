import MDAnalysis as mda
from rdkit.Chem import Descriptors as rdkit_desc
from MDAnalysis.core.universe import Universe



class LigandProperty:
    def __init__(self, mol_obj, descriptors):
        if isinstance(mol_obj, mda.core.groups.AtomGroup):
            self.ligand = mol_obj.convert_to("RDKIT")
        else:
            self.ligand = mol_obj
        self.check_descriptors(descriptors)
        self.descriptors = descriptors

    def descriptor_function_mapper(func):
        def wrapper(*args, **kwargs):
            mapper = {
                "mol_weigth": rdkit_desc.MolWt,
                "logp": rdkit_desc.MolLogP,
                "hba": rdkit_desc.NumHAcceptors,
                "hbd": rdkit_desc.NumHDonors,
                "tpsa": rdkit_desc.TPSA,
                "mr": rdkit_desc.MolMR,
                "rot_bonds": rdkit_desc.NumRotatableBonds
            }
            return func(*args, **kwargs, mapper=mapper)

        return wrapper

    def calculate(self):
        data = {desc: self.calculate_descriptor(self.ligand, desc) for desc in self.descriptors}
        setattr(self, "data", data)

    def get_data(self):
        val = getattr(self, 'data', {})
        return val

    @classmethod
    @descriptor_function_mapper
    def calculate_descriptor(cls, mol, descriptor, mapper):
        descriptor_value = mapper[descriptor](mol)
        # print(mapper[descriptor](m))
        # print(m, mapper)
        return descriptor_value

    @classmethod
    @descriptor_function_mapper
    def check_descriptors(cls, descriptors, mapper):
        unknown_descriptors = set(descriptors) - set(mapper.keys())
        if unknown_descriptors:
            raise NameError(f"Unknown descriptors: {', '.join(unknown_descriptors)}")

    @classmethod
    @descriptor_function_mapper
    def list_avaliable_descriptors(cls, mapper):
        return list(mapper.keys())


