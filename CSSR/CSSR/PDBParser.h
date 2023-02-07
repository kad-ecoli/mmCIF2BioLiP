/* parse PDB file into data structure similar to Bio.PDB in biopython
 * (model - chain - residue - atom). */
#ifndef PDBParser_HPP
#define PDBParser_HPP 1

#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

struct AtomUnit    // struct for each atom entry
{
    string name;       // atom name
    vector<float> xyz; // coordinate
};

struct ResidueUnit // struct for each residue
{
    bool het;               // true - HETATM, false - ATOM
    int resi;               // residue sequence number
    char icode;             // insertion code
    string resn;            // residue name
    vector<AtomUnit> atoms; // list of atoms
};

struct ChainUnit  // struct for each chain
{
    string chainID_full;          // chain ID, might be more than 1 char
    char chainID;                 // short chain ID, must be 1 char
    string sequence;              // sequence converted from CA coordinate
    string sarst;                 // SARST (Structural similarity search
                                  // Aided by Ramachandran Sequential 
                                  // Transformation) code
    string ss;                    // secondary structure
    vector<ResidueUnit> residues; // list of residues
};

struct ModelUnit  // struct for each model in mult-model PDB
{
    vector<ChainUnit> chains; // list of chains
};

int parse_pdb_line(const string line,ModelUnit &pep, ChainUnit &chain,
    ResidueUnit &residue, AtomUnit &atom, map<char,string> &chainIDmap,
    const int atomic_detail,const int allowX);
ModelUnit read_pdb_structure(const char *filename,
    const int atomic_detail,const int allowX);
string write_pdb_structure(ChainUnit &chain,int &i);
void write_pdb_structure(const char *filename,ChainUnit &chain);
string write_pdb_structure(ModelUnit &pep);
void write_pdb_structure(const char *filename,ModelUnit &pep);
void reindex_pdb(const int startindex,ChainUnit& chain);
void reindex_pdb(const int startindex,ModelUnit& pep);
char aa3to1(const string resn,const int convertX);
string pdb2fasta(ChainUnit& chain);
string pdb2fasta(ModelUnit& pep,const string PDBid,const int ShowSeqLen);
void remove_sidechain(ResidueUnit& residue,int atomic_detail);
void remove_sidechain(ChainUnit& chain,int atomic_detail);
void remove_sidechain(ModelUnit& pep,int atomic_detail);
#endif
