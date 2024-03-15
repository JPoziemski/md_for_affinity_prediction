
from collections import Counter, defaultdict

from scipy.spatial import ConvexHull
class ShapeDescriptor:
    def __init__(self, universe, mol_obj, descriptors=("area", "volume"), suffix=""):
        self.mol_obj = mol_obj
        self.suffix = suffix
        #print(self.mol_obj.positions)
        self.descriptors = descriptors
        self.universe = universe

    def calculate(self):
        data = defaultdict(list)
        for _ in self.universe.trajectory:
            hull = ConvexHull(self.mol_obj.positions)
            for desc in self.descriptors:
                data[f"{desc}{self.suffix}"].append(getattr(hull, desc))
        self.data = dict(data)
        self.clean_data()

    def clean_data(self):
        if len(self.universe.trajectory) ==1:
            self.data = {key: value[0] for key, value in self.data.items()}

    def get_data(self):
        return self.data