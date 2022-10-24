A local copy of ligand and binding site information created from mmCIF files of the PDB database.

### Step 0: Initialize database ###
This step is only needed the first time the database is created. It will not be necessary for subsequent database update.
```bash
script/rsyncPDB.sh
```
This will download the initial set of mmCIF files for full asymetric unit to pdb/data/structures/divided/.

### Step 1: Generate table for each type of ligand ###
```bash
script/download_ligand.pl
```
This will download CCD ligand to pdb/data/monomers/components.cif.gz and extract their summary to data/ligand.tsv.gz.

## Step 2: Download SIFTS function annotation ###
```bash
./script/download_sifts.pl
```
This will download taxonomy, uniprot accession, EC, GO and pubmed ID.
