
Each line in the annotation file BioLiP_*.txt annotates for each ligand-protein interaction site in BioLiP. 
The columns are separated by the tab key, i.e., "\t", so that you can easily import the data into MySQL database with command like this:

load data local infile 'BioLiP.dat' into TABLE biolip FIELDS TERMINATED BY '\t';



Example entry: 
966c    A       1.90    BS06    RS2     A       1       N180 L181 A182 V215 H218 E219 H222 H228 L235 Y237 P238 S239 Y240 T241   N73 L74 A75 V108 H111 E112 H115 H121 L128 Y130 P131 S132 Y133 T134      M236 E219;E219  M129 E112;E112  3.4.24.-        0004222,0006508,0008237,0008270,0031012         ki=23nM (RS2)   Ki=23nM (RS2)           P03956  10074939        RWEQTHLTYRIENYTPDLPRADVDHAIEKAFQLWSNVTPLTFTKVSEGQADIMISFVRGDHRDNSPFDGPGGNLAHAFQPGPGIGGDAHFDEDERWTNNFREYNLHRVAAHELGHSLGLSHSTDIGALMYPSYTFSGDVQLAQDDIDGIQAIYGRSQ

The columns are (from left to right):
01	PDB ID
02	Receptor chain
03	Resolution. "-1.00" stands for lack of resolution information, e.g. for NMR
04	Binding site number code
05	Ligand ID in the Chemical Component Dictionary (CCD) used by the PDB database
06	Ligand chain
07	Ligand serial number
08      Binding site residues (with PDB residue numbering)
09      Binding site residues (with residue re-numbered starting from 1)
10	Catalytic site residues (different sites are separated by ';') (with PDB residue numbering)
11      Catalytic site residues (different sites are separated by ';') (with residue re-numbered starting from 1)
12	EC number
13	GO terms
14	Binding affinity by manual survey of the original literature. The information in '()' is the PubMed ID
15	Binding affinity provided by the Binding MOAD database. The information in '()' is the ligand information in Binding MOAD
16	Binding affinity provided by the PDBbind-CN database. The information in '()' is the ligand information in PDBbind-CN
17	Binding affinity provided by the BindingDB database
18	UniProt ID
19	PubMed ID
20	Residue sequence number of the ligand (field _atom_site.auth_seq_id in PDBx/mmCIF format)
21	Receptor sequence

The ligand-protein complex structure for each line entry be obtained by using the two PDB files:
1. the receptor structure file under the "receptor" folder (name is formed with columns 01,02: i.e., 0102.pdb)
2. the corresponding ligand structure file under the "ligand" folder (name is formed with columns 01,05,06,07: i.e., 01_05_06_07.pdb)

For example, the complex structure for the above example can be obtained using the two files:
1. 966cA.pdb under the "receptor" folder
2. 966c_RS2_A_1.pdb under the "ligand" folder


Please kindly cite the following article if you use BioLiP:

Jianyi Yang, Ambrish Roy, Yang Zhang, BioLiP: a semi-manually curated database for biologically relevant ligand-protein interactions, Nucleic Acids Research, 41:D1096-D1103, 2013. 


