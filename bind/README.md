## Technical Note on Generating PDBbind-CN data ##

PDBbind-CN (together with BindingDB and MOAD) is one of the three major databases for binding affinities in PDB. PDBbind-CN is updated annually. This document explains how to semi-automatically process PDBbind-CN data, which is different to handle fully automatically for the following reasons:

1. PDBbind-CN data can only be downloaded after login.
2. The file names are irregular.
3. The receptor information does not contain any chain ID.
4. The ligand is not always clearly defined.

#### Step 1. download raw data ####

Raw data can be downloaded from http://www.pdbbind.org.cn/download.php after login.

1. Download "Index files of PDBbind" (PDBbind_****_plain_text_index.tar.gz).
2. Download "Protein-ligand complexes: The general set minus refined set" (PDBbind_****_PL.tar.gz)
3. Download "Protein-ligand complexes: The refined set" (PDBbind_****_refined.tar.gz)

#### Step 2. decompress data ####

1. Extract the following files from "Index files of PDBbind":
```
index/INDEX_general_PL_data.*
index/INDEX_general_PL_PN.*
```
For release 2020, the file names are
```
index/INDEX_general_PL_data.2020
index/INDEX_general_PN.2020
```

2. Extract files endding with ligand.mol2 and pocket.pdb from "Protein-ligand complexes" and place them under:
```
pocket/*_pocket.pdb
ligand/*_ligand.mol2
```

#### Step 3. parse data ####
```bash
../script/curate_PDBbind.pl
```
This script write ../data/PDBbind.csv

#### Step 4. clean up ####
```bash
rm index/INDEX_*
rm pocket/*_pocket.pdb
rm ligand/*_ligand.mol2
```
