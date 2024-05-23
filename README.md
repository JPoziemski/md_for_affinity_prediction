# MD for Affinity Prediction

\
Repository used to generate the data used in the paper [1]. 
Complexes from the pdbbind v2020 dataset were used to generate data from the crystallographic structures. 
The molecular dynamics trajectories used to generate the data are deposited at:  https://zenodo.org/records/11172815.

Usage:\
Use the following scripts to generate the data:
- compute_descriptors_for_cystal_structures.py - generates data dal crystallographic complexes
- compute_descriptors_for_trajectory.py - generates data for trajactories of molecular dynamics
- compute_ligand_static_descriptors.py - generates static ligand decryptors

Requirements:
- python 3.10, 
- Packages: prolif, rdkit, pandas, MDAnalysis, 


References:
1. Poziemski J, Yurkevych A, Siedlecki P. Assessment of molecular dynamics time series descriptors in protein-ligand affinity prediction. ChemRxiv. 2024; doi:10.26434/chemrxiv-2024-dxv36 


