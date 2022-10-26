A local copy of ligand and binding site information created from mmCIF files of the PDB database.

## Installation ##
```bash
cd script/
make
```
The compilation requires GCC with C++ 11 support.
The binary executable of [cd-hit](https://github.com/weizhongli/cdhit) is available at ``script/cd-hit`` for 64 bit Linux.
Reinstall cd-hit if using other operating system.

### Step 0: Initialize database ###
This step is only needed the first time the database is created. It will not be necessary for subsequent database update.
```bash
script/rsyncPDB.sh
```
This script downloads the initial set of mmCIF files for full asymetric unit to ``pdb/data/structures/divided/mmCIF``.

### Step 1: Download SIFTS function annotation ###
```bash
./script/download_sifts.pl
```
This script downloads taxonomy, uniprot accession, EC, GO and pubmed ID to ``sifts/flatfiles/tsv`` and extract the data to ``data/*.tsv.gz``.
This script downloads CCD ligand to ``pdb/data/monomers/components.cif.gz`` and extract their summary to ``data/ligand.tsv.gz``.

Optionally, run the following script to download csa. The manually curated dataset of csa is updated infrequently. Therefore, it is not necessary to run it every week.
```bash
./script/download_csa.pl
```
This script downloads catalytic site to ``m-csa/api/residues.json`` and extract the summary to ``data/csa.tsv.gz``.

### Step 2: Download binding affinity ###
```bash
./script/download_bind.pl
```
This script downloads BindingDB to ``bind/BindingDB*`` and extract the summary to ``data/BindingDB.tsv.gz``

Optionally, run the following script to download MOAD. Since MOAD only updates every few years, it is not necessary to run it weekly.
```bash
./script/download_moad.pl
```
This script downloads MOAD to ``bind/every.csv`` and extract the summary to ``data/moad.tsv``.


### Step 3: Download missing PDB entries ###
```bash
./script/download_pdb.pl
```
This script downloads the set of missing mmCIF files for full asymetric unit to ``pdb/data/structures/divided/mmCIF``.

### Step 4: convert mmCIF to PDB format ###
```bash
./script/curate_pdb.pl
```
This script converts mmCIF files from ``pdb/data/structures/divided/mmCIF`` to ``interim/*/*.tar.gz`` and ``interim/*/*.txt``

### Step 5: download pubmed abstract ###
```bash
./script/download_pubmed.pl
```
This script checks ``interim/*/*.txt`` for artifact ligand and download pubmed abstract to ``pubmed/*.txt``

### Step 6: remove artifact and unbound ligand ###
```bash
./script/curate_ligand.pl
```
This script check potential artifact ligand listed by ``interim/*/*.txt`` against pubmed abstract at ``pubmed/*.txt``.
It repackages receptor and biologically relevent ligand pdb files into ``interim/*/*.tar.bz2``. 
Binding sites are written to ``interim/*/*.bsr``.

### Step 7: prepare redundant dataset ###
```bash
./script/make_weekly.pl
```
This script reads ``interim/*/*.tar.bz2`` and ``interim/*/*.bsr`` and writes ``weekly/receptor_*.tar.bz2``, ``weekly/receptor1_*.tar.bz2``, ``weekly/ligand_*.tar.bz2`` and ``weekly/BioLiP_*.bsr.gz``.
FASTA sequence of protein, peptide, rna, dna are saved to ``data/*.fasta.gz``

#### Step 8: prepare nonredundant dataset ###
```bash
./script/make_nr.pl
```
This script converts ``weekly/BioLiP_*.bsr.gz`` to ``weekly/BioLiP_*.txt`` and ``weekly/BioLiP_*_nr.txt``.
It also creates ``weekly/receptor_*_nr.tar.bz2``, ``weekly/receptor1_*_nr.tar.bz2`` and ``weekly/ligand_*_nr.tar.bz2`` from intermediate files of the previous step.
