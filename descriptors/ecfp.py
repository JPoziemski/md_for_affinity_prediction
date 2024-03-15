from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect


class ECFP:
    def __init__(self, mol_obj):
        self.mol_obj = mol_obj


    @classmethod
    def calculate_ecfp(cls, molecule, radius=2, nBits=1024):
        ecfp4 = GetMorganFingerprintAsBitVect(molecule, radius=radius, nBits=nBits)
        binary_vector = ecfp4.ToBitString()
        binary_vector_tuple = tuple(int(bit) for bit in binary_vector)
        return binary_vector_tuple

    def calculate(self):
        ecfp_vector = self.calculate_ecfp(self.mol_obj)
        ecfp_dict = {f"ECFP_{str(i)}": bit for i,bit in enumerate(ecfp_vector)}
        self.data = ecfp_dict

    def get_data(self):
        return self.data