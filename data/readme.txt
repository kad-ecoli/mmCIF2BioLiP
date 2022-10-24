ligand.tsv.gz - ligand information in the following format
                [1] CCD ID (i.e., residue name)
                [2] formula
                [3] InCHI
                [4] InChIKey
                [5] SMILES
                [6] name and synonym

chain2ec.tsv.gz - chain to enzyme commission number mapping
                [1] pdb ID
		[2] chain ID
		[3] EC number (comma separated list)

chain2go.tsv.gz - chain to gene ontology term mapping
                [1] pdb ID
		[2] chain ID
		[3] GO term (comma separated list)

chain2taxonomy.tsv.gz - chain to NCBI taxonomy mapping
                [1] pdb ID
		[2] chain ID
		[3] taxon ID (comma separated list)

taxid2name.tsv.gz - scientific name for each taxon
                [1] taxon ID
		[2] scientific name

chain2uniprot.tsv.gz - chain to uniprot accession mapping
                [1] pdb ID
		[2] chain ID
		[3] uniprot accession (comma separated list)

pdb2pubmed.tsv.gz - primary citation pubmed ID
                [1] pdb ID
		[2] pubmed ID
