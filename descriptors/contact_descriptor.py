
import MDAnalysis as mda
from .descriptor import Descriptor
from MDAnalysis.analysis import contacts
from functools import partial

class ContactDescriptor(Descriptor):

    def __init__(self, universe, ligand_code,
                 descriptors=("contact_pocket_new", "contact_pocket_old", "contact_pocket_all",
                              "contact_ligand_pocket_new", "contact_ligand_pocket_old", "contact_ligand_pocket_all")):
        self.universe = universe
        self.ligand_code = ligand_code
        self.pocket = universe.select_atoms(f"around 6 resname {ligand_code}").residues.atoms
        self.pocket = self.pocket.select_atoms(f"not type H")
        self.ligand = self.universe.select_atoms(f"resname {ligand_code} and not type H")
        self.descriptors = self.set_descriptors(descriptors)


    def _calculate_pocket_conatcts(self, contact_type=None):
        if contact_type == "all":
            contact_num = self.contacts_within_cutoff(self.universe, self.pocket, self.pocket)
        elif contact_type == "old" and len(self.universe.trajectory) > 1:
            self.universe.trajectory[-1]
            pocket = self.universe.select_atoms(f"around 6 resname {self.ligand_code}").residues.atoms
            pocket = pocket.select_atoms(f"not type H")
            contact_num = self._compute_new_old_contacts(self.universe, pocket, pocket)
            self.universe.trajectory[0]
        elif contact_type == "old" and len(self.universe.trajectory) > 1:
            contact_num = self._compute_new_old_contacts(self.universe, self.pocket, self.pocket)
        else:
            pass
        return contact_num

    def _calculate_pocket_ligand_contact(self, contact_type):
        if contact_type == "all":
            contact_num = self.contacts_within_cutoff(self.universe, self.ligand, self.pocket)
        elif contact_type == "old":
            self.universe.trajectory[-1]
            pocket = self.universe.select_atoms(f"around 6 resname {self.ligand_code}").residues.atoms
            pocket = pocket.select_atoms(f"not type H")
            ligand = self.universe.select_atoms(f"resname {self.ligand_code} and not type H")

            contact_num = self._compute_new_old_contacts(self.universe, ligand, pocket)
            self.universe.trajectory[0]
        else:
            contact_num = self._compute_new_old_contacts(self.universe, self.ligand, self.pocket)

        return contact_num

    @classmethod
    def contacts_within_cutoff(cls, u, group_a, group_b, radius=4.5):
        timeseries = []
        for ts in u.trajectory:
            # calculate distances between group_a and group_b
            dist = contacts.distance_array(group_a.positions, group_b.positions)
            # determine which distances <= radius
            n_contacts = contacts.contact_matrix(dist, radius).sum()
            timeseries.append(n_contacts)
        return timeseries

    def _compute_new_old_contacts(self, u, ligand, pocket, radius=4.5):
        contacts_data = contacts.Contacts(u,
                                          select=(pocket, ligand),
                                          refgroup=(pocket, ligand),
                                          radius=radius,
                                          method='radius_cut').run()
        n_ref = contacts_data.initial_contacts[0].sum()
        n_contacts = contacts_data.results.timeseries[:, 1] * n_ref
        return n_contacts

    def get_descriptors_mapper(self):
        mapper = {
            "contact_pocket_new": partial(self._calculate_pocket_conatcts, "new"),
            "contact_pocket_old": partial(self._calculate_pocket_conatcts, "old"),
            "contact_pocket_all": partial(self._calculate_pocket_conatcts, "all"),
            "contact_ligand_pocket_new": partial(self._calculate_pocket_ligand_contact, "new"),
            "contact_ligand_pocket_old": partial(self._calculate_pocket_ligand_contact, "old"),
            "contact_ligand_pocket_all": partial(self._calculate_pocket_ligand_contact, "all")
        }
        return mapper

    def calculate(self):
        mapper = self.get_descriptors_mapper()
        self.data = {descriptor: mapper[descriptor]() for descriptor in self.descriptors}
        if len(self.universe.trajectory) == 1:
            self.data = {key: value[0] for key, value in self.data.items()}

    def set_descriptors(self, descriptors):
        trajectory_desc = ["contact_ligand_pocket_old", "contact_ligand_pocket_new", 'contact_pocket_new', 'contact_pocket_old']
        cleaned_descriptors = []
        for desc in descriptors:
            if desc in trajectory_desc:
                if len(self.universe.trajectory) >1:
                    cleaned_descriptors.append(desc)
                else:
                    continue
            else:
                cleaned_descriptors.append(desc)

        return cleaned_descriptors

    def get_data(self):
        return self.data