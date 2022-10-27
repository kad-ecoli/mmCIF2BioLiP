A local copy of ligand and binding site information created from mmCIF files of the PDB database.

## Installation ##
```bash
cd script/
make
```
The compilation requires GCC with C++ 11 support.
The binary executable of [cd-hit](https://github.com/weizhongli/cdhit) is available at ``script/cd-hit`` for 64 bit Linux.
The binary executable of makeblastdb, blastn and blastp from the [NCBI BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) are available at ``script/`` for 64 bit Linux.
Reinstall cd-hit, makeblastdb, blastn and blastp if using other operating system.

## Usage ##

### Step 1: Download PDB entries ###
Run the following shell script the first time the database is created.
```bash
script/rsyncPDB.sh
```
This script downloads the initial set of mmCIF files for full asymetric unit to ``pdb/data/structures/divided/mmCIF``. It will also download the resolution of pdb to ``pdb/derived_data/index/resolu.idx``.

Alternatively, run the following script to download only the missing mmcif files and resolu.idx.
```bash
script/download_pdb.pl
```
Both ``script/rsyncPDB.sh`` and ``script/download_pdb.pl`` have the same purpose of downloading mmcif files and resolu.idx. For the first time the database is setup, ``script/rsyncPDB.sh`` must be used instead of ``script/download_pdb.pl`` because the later is much slower when a large number of files need be downloaded. For subsequent weekly update, either one of the script can be used. If you do not want to clean up the pdb/ folder, the first script ``script/rsyncPDB.sh`` is preferred. If you do want to clean up the pdb/ folder after every update, which can save ~60GB of disk space, you must use the second script ``script/download_pdb.pl`` to avoid re-downloading the whole pdb database every week.

### Step 2: Download SIFTS function annotation ###
This step can be performed in parallel to step 1.
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

### Step 3: Download binding affinity ###
This step must be performed after step 2 (regardless of step 1) because it requires the pdb to uniprot mapping file generated in the previous step.
```bash
./script/download_bind.pl
```
This script downloads BindingDB to ``bind/BindingDB*`` and extract the summary to ``data/BindingDB.tsv.gz``.
This script may have the warning: "gzip: bind/BindingDB_All_tsv.zip: invalid compressed data--length error", which is caused by decompressing zip with zcat. This warning can be ignored.

Optionally, run the following script to download MOAD. Since MOAD only updates every few years, it is not necessary to run it weekly.
```bash
./script/download_moad.pl
```
This script downloads MOAD to ``bind/every.csv`` and extract the summary to ``data/moad.tsv``.

Again optionally, follow instructions at [bind/README.md](bind/README.md) for how to semi-manually create PDBbind-CN summary at ``data/PDBbind.tsv``. This only need to be peformed annually.

### Step 4: convert mmCIF to PDB format ###
This step must be performed after step 1 (regardless of step 2 and 3).
```bash
./script/curate_pdb.pl
```
This script converts mmCIF files from ``pdb/data/structures/divided/mmCIF`` to ``interim/*/*.tar.gz`` and ``interim/*/*.txt``

### Step 5: download pubmed abstract ###
This step must be performed after step 2 and 4 (regardless of step 3).
```bash
./script/download_pubmed.pl
```
This script checks ``interim/*/*.txt`` for artifact ligand and download pubmed abstract to ``pubmed/*.txt``

### Step 6: remove artifact and unbound ligand ###
This step must be performed after step 5.
```bash
./script/curate_ligand.pl
```
This script check potential artifact ligand listed by ``interim/*/*.txt`` against pubmed abstract at ``pubmed/*.txt``.
It repackages receptor and biologically relevent ligand pdb files into ``interim/*/*.tar.bz2``. 
Binding sites are written to ``interim/*/*.bsr``.

### Step 7: prepare redundant dataset ###
This step must be performed after step 6.
```bash
./script/make_weekly.pl
```
This script reads ``interim/*/*.tar.bz2`` and ``interim/*/*.bsr`` and writes ``weekly/receptor_*.tar.bz2``, ``weekly/receptor1_*.tar.bz2``, ``weekly/ligand_*.tar.bz2`` and ``weekly/BioLiP_*.bsr.gz``.
FASTA sequence of protein, peptide, rna, dna are saved to ``data/*.fasta.gz``

### Step 8: prepare nonredundant dataset ###
This step must be performed after step 3 and 7.
```bash
./script/make_nr.pl
```
This script converts ``weekly/BioLiP_*.bsr.gz`` to ``weekly/BioLiP_*.txt`` and ``weekly/BioLiP_*_nr.txt``.
It also creates ``weekly/receptor_*_nr.tar.bz2``, ``weekly/receptor1_*_nr.tar.bz2`` and ``weekly/ligand_*_nr.tar.bz2`` from intermediate files of the previous step.
