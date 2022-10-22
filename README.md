A local copy of ligand and binding site information created from mmCIF files of the PDB database.

### Step 0: Initialize database ###
This step is only needed the first time the database is created. It will not be necessary for subsequent database update.
```bash
script/rsyncPDB.sh
```
This will download the initial set of mmCIF files for full asymetric unit and ligands to pdb/data/structures/divided/ and pdb/refdata/chem_comp/, respectively.
