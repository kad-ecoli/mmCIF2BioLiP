A local copy of ligand and binding site information created from mmCIF files of the PDB database.

### Step 0: Initialize database ###
This step is only needed the first time the database is created. It will not be necessary for subsequent database update.
```bash
script/rsyncPDB.sh
```
This will download the initial set of mmCIF files for full asymetric unit to pdb/data/structures/divided/.

### Step 1: Generate table for each type of ligand ###
```bash
script/rsyncLigand.sh
script/ligand2tsv.pl
```
This will download all types of ligand to pdb/refdata/chem_comp/ and extract their summary to data/ligand.tsv.gz.
