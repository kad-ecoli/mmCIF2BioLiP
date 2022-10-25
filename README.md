A local copy of ligand and binding site information created from mmCIF files of the PDB database.

## Installation ##
```bash
cd script/
make
```
The compilation requires GCC with C++ 11 support.

### Step 0: Initialize database ###
This step is only needed the first time the database is created. It will not be necessary for subsequent database update.
```bash
script/rsyncPDB.sh
```
This script downloads the initial set of mmCIF files for full asymetric unit to pdb/data/structures/divided/mmCIF.

## Step 1: Download SIFTS function annotation ###
```bash
./script/download_sifts.pl
```
This script downloads CCD ligand to pdb/data/monomers/components.cif.gz and extract their summary to data/ligand.tsv.gz.
This script downloads taxonomy, uniprot accession, EC, GO and pubmed ID to sifts/flatfiles/tsv and extract the data to data/\*.tsv.gz.

### Step 2: Download missing PDB entries ###
```bash
./script/download_pdb.pl
```
This script downloads the set of missing mmCIF files for full asymetric unit to pdb/data/structures/divided/mmCIF.

### Step 3: convert mmCIF to PDB format ###
```bash
./script/curate_pdb.pl
```
This script converts mmCIF files from pdb/data/structures/divided/ed pdb/data/structures/divided/mmCIF to interim/*/*.tar.gz and interim/*/*.txt

### Step 4: download pubmed abstract ###
```bash
./script/download_pubmed.pl
```
This script checks interim/*/*.txt for artifact ligand and download pubmed
abstract for corresponding entry.

### Step 5: remove artifact and unbound ligand ###
```bash
./script/curate_ligand.pl
```
This script check potential artifact ligand against pubmed abstract. It repackages receptor and biologically relevent ligand pdb files into interim/*/*.tar.bz2. Binding sites are written to interim/*/*.bsr.

### Step 6: prepare redundant dataset ###
```bash
./script/make_weekly.pl
```
