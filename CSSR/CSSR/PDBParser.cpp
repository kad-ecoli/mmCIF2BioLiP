/* parse PDB file into data structure similar to Bio.PDB in biopython
 * (model - chain - residue - atom). */
#ifndef PDBParser_CPP
#define PDBParser_CPP 1

#include "PDBParser.h"
#include "pstream.h"

using namespace std;

/* parse one line in PDB file, append the data to pep. 
 * used by read_pdb_structure
 * allowX: 0 - ATOM, 1 - ATOM and MSE, converting MSE to MET
 *         2 - all, converting MSE to MET, 3 - all, no conversion
 *
 * return 0 if line not parsed, 1 if line is parsed
 */
int parse_pdb_line(const string line,ModelUnit &pep, ChainUnit &chain,
    ResidueUnit &residue, AtomUnit &atom, map<char,string> &chainIDmap,
    const int atomic_detail=2,const int allowX=1)
{
    string record_name=line.substr(0,6);
    char altLoc=line[16];

    atom.name=line.substr(12,4);
    residue.resn=line.substr(17,3);

    if ((allowX==0 && record_name!="ATOM  ")||
        (allowX==1 && record_name!="ATOM  " &&  
            !(record_name=="HETATM" && residue.resn=="MSE"))||
        (allowX>=2 && record_name!="ATOM  " && record_name!="HETATM"))
        return 0;
    
    // ignore alternatively locating residues
    if (altLoc!=' ' && altLoc!='A') return 0;

    if ((atomic_detail==0 && atom.name!=" CA " && atom.name!=" C3'")||
        (atomic_detail==1 && atom.name!=" CA " && atom.name!=" C3'" &&
         atom.name!=" N  "&& atom.name!=" C  " && atom.name!=" O  ")) return 0;
    
    if (residue.resn=="MSE" && allowX<3)
    {
        record_name="ATOM  ";
        residue.resn="MET";
    }
    residue.het=(record_name=="HETATM");
    chain.chainID=line[21];
    if (chainIDmap.find(chain.chainID)==chainIDmap.end())
        chain.chainID_full=chain.chainID;
    else
        chain.chainID_full=chainIDmap[chain.chainID];
    if (chain.chainID_full==" ") chain.chainID_full="_";
    residue.resi=atoi(line.substr(22,4).c_str());
    residue.icode=line[26];
    atom.xyz[0]=atof(line.substr(30,8).c_str());
    atom.xyz[1]=atof(line.substr(38,8).c_str());
    atom.xyz[2]=atof(line.substr(46,8).c_str());

    int chain_index=-1;
    for (size_t c=0;c<pep.chains.size();c++)
        if (pep.chains[c].chainID_full==chain.chainID_full) chain_index=c;
    if (chain_index==-1)
    {
        pep.chains.push_back(chain);
        chain_index=pep.chains.size()-1;
    }

    if (pep.chains[chain_index].residues.size()==0||
        pep.chains[chain_index].residues.back().resi !=residue.resi||
        pep.chains[chain_index].residues.back().icode!=residue.icode)
        pep.chains[chain_index].residues.push_back(residue);
    
    pep.chains[chain_index].residues.back().atoms.push_back(atom);
    return 1;
}

/* atomic_detail: 0 - CA only, 1 - backbone heavy atoms (CA C N O), 2 - all atom
 * allowX: 0 - ATOM, 1 - ATOM and MSE, converting MSE to MET, 
 *         2 - all, converting MSE to MET, 3 - all, no conversion
 * filename: full filename path, stdin if filename=="-"
 */
ModelUnit read_pdb_structure(const char *filename,
    const int atomic_detail=2,const int allowX=1)
{
    ModelUnit pep;

    string line="";
    string record_name="ATOM  ";

    AtomUnit atom;
    atom.xyz.assign(3,0);

    ResidueUnit residue;

    ChainUnit chain;

    map<char,string> chainIDmap;
    string filename_str=(string) filename;
    
#ifdef REDI_PSTREAM_H_SEEN
    if (filename_str.length()>=18 &&
        filename_str.substr(filename_str.length()-18,18)=="-pdb-bundle.tar.gz")
    {
        // best effort/minimal tarball
        redi::ipstream fp_gz("tar -xOzf "+filename_str+
            " --wildcards *-chain-id-mapping.txt");
        map <string,map<char,string> > PDBmap;
        vector<string> PDBvec;

        string PDBfile=""; // PDB format file in tarball
        char chainID;
        string chainID_full;

        while(fp_gz.good())
        {
            getline(fp_gz,line);
            if (line.length()==0 || (PDBfile.length()==0 && line[0]==' '))
                continue;
            if (line[0]!=' ') // new PDBfile
            {
                PDBfile=line.substr(0,line.length()-1);
                PDBmap[PDBfile]=chainIDmap;
                PDBvec.push_back(PDBfile);
            }
            else // new chain
            {
                istringstream ss(line);
                ss>>chainID>>chainID_full;
                PDBmap[PDBfile][chainID]=chainID_full;
            }
        }
        fp_gz.close();

        for (size_t i=0;i<PDBvec.size();i++)
        {
            PDBfile=PDBvec[i];
            redi::ipstream fp_gz2("tar -xOzf "+filename_str+' '+PDBfile);
            while(fp_gz2.good())
            {
                getline(fp_gz2,line);
                if (line.substr(0,3)=="END") break;
                if (line.length()<53) continue;
                parse_pdb_line(line,pep,chain,residue,atom,PDBmap[PDBfile],
                    atomic_detail,allowX);
            }
            fp_gz2.close();
        }

        chain.residues.clear();
        residue.atoms.clear();
        chainIDmap.clear();
        PDBmap.clear();
        return pep;
    }
#endif
    

    int use_stdin=(filename_str=="-");
    int use_pstream=0; // input is compressed

    ifstream fp;
#ifndef REDI_PSTREAM_H_SEEN
    ifstream fp_gz;
#else
    redi::ipstream fp_gz; // if file is compressed
    if (filename_str.length()>=3 && 
        filename_str.substr(filename_str.length()-3,3)==".gz")
    {
        // gzip pdb
        fp_gz.open("zcat "+filename_str);
        use_pstream=1;
    }
    else
#endif
    {
        fp.open(filename,ios::in); //ifstream fp(filename,ios::in);
    }

    while(use_stdin?cin.good():(use_pstream?fp_gz.good():fp.good()))
    {
        if (use_stdin)
            getline(cin,line);
        else if (use_pstream)
            getline(fp_gz,line);
        else
            getline(fp,line);

        if (line.substr(0,3)=="END") break;
        if (line.length()<53) continue;
        
        parse_pdb_line(line,pep,chain,residue,atom,chainIDmap,
            atomic_detail,allowX);
    }
    if (!use_stdin)
    {
        if (use_pstream==0)
            fp.close();
        else
            fp_gz.close();
    }
    
    chain.residues.clear();
    residue.atoms.clear();
    chainIDmap.clear();
    return pep;
}

/* i - first atom index */
string write_pdb_structure(ChainUnit &chain,int &i)
{
    stringstream buf;
    size_t r,a;

    for (r=0;r<chain.residues.size();r++)
        for (a=0;a<chain.residues[r].atoms.size();a++)
            buf<<setiosflags(ios::left)<<setw(6)
               <<(chain.residues[r].het?"HETATM":"ATOM")
               <<resetiosflags(ios::left)<<setw(5)<<i++<<' '
               <<chain.residues[r].atoms[a].name<<' '
               <<chain.residues[r].resn<<' '<<chain.chainID<<setw(4)
               <<chain.residues[r].resi<<chain.residues[r].icode
               <<"   " <<setiosflags(ios::fixed)<<setprecision(3)
               <<setw(8)<<chain.residues[r].atoms[a].xyz[0]
               <<setw(8)<<chain.residues[r].atoms[a].xyz[1]
               <<setw(8)<<chain.residues[r].atoms[a].xyz[2]<<endl;
    return buf.str();
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_pdb_structure(const char *filename,ChainUnit &chain)
{
    int i=1;
    if (strcmp(filename,"-")==0)
        cout<<write_pdb_structure(chain,i);
    else
    {
        ofstream fp(filename);
        fp<<write_pdb_structure(chain,i);
        fp.close();
    }
}

string write_pdb_structure(ModelUnit &pep)
{
    string txt="";
    int i=1; // atom index
    for (size_t c=0;c<pep.chains.size();c++)
        txt+=write_pdb_structure(pep.chains[c],i)+"TER\n";
    return txt;
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_pdb_structure(const char *filename,ModelUnit &pep)
{
    if (strcmp(filename,"-")==0)
        cout<<write_pdb_structure(pep);
    else
    {
        ofstream fp(filename);
        fp<<write_pdb_structure(pep);
        fp.close();
    }
}

/* renumber residue index from "startindex" */
void reindex_pdb(const int startindex,ChainUnit& chain)
{
    for (int r=0;r<chain.residues.size();r++)
        chain.residues[r].resi=r+startindex;
}

void reindex_pdb(const int startindex,ModelUnit& pep)
{
    for (size_t c=0;c<pep.chains.size();c++) 
        reindex_pdb(startindex,pep.chains[c]);
}

/* convert pdb structure to fasta sequence 
 * convertX - how to deal with non-standard amino acids
 *            0 - only 20 standard amino acids
 *            1 - 20 standard amino acids + MSE
 *            2 - non-standard amino acid with known parent, 
 *                all to legal amino acid in BLOSUM
 *            3 - non-standard amino acid with known parent
 */
char aa3to1(const string resn,const int convertX=2)
{
    // 20 standard amino acid + MSE
    if (resn[0]==' ' && (resn[1]=='D'||resn[1]==' ')) return tolower(resn[2]);
    if (resn=="ALA") return 'A';
    if (resn=="CYS") return 'C';
    if (resn=="ASP") return 'D';
    if (resn=="GLU") return 'E';
    if (resn=="PHE") return 'F';
    if (resn=="GLY") return 'G';
    if (resn=="HIS") return 'H';
    if (resn=="ILE") return 'I';
    if (resn=="LYS") return 'K';
    if (resn=="LEU") return 'L';
    if (resn=="MET") return 'M';
    if (resn=="ASN") return 'N';
    if (resn=="PRO") return 'P';
    if (resn=="GLN") return 'Q';
    if (resn=="ARG") return 'R';
    if (resn=="SER") return 'S';
    if (resn=="THR") return 'T';
    if (resn=="VAL") return 'V'; 
    if (resn=="TRP") return 'W';
    if (resn=="TYR") return 'Y';

    if (resn=="MSE" && convertX>=1) return 'M';

    if (convertX>=2)
    {
        // non-standard amino acid with known parent
        if (resn=="CHG"||resn=="HAC"||resn=="AYA"||resn=="TIH"||resn=="BNN"||
            resn=="ALM"||resn=="TPQ"||resn=="MAA"||resn=="PRR"||resn=="FLA"||
            resn=="AIB"||resn=="DAL"||resn=="CSD"||resn=="DHA"||resn=="DNP") 
            return 'A';
        if (resn=="PR3"||resn=="CCS"||resn=="C6C"||resn=="SMC"||resn=="BCS"||
            resn=="SCY"||resn=="DCY"||resn=="SCS"||resn=="CME"||resn=="CY1"||
            resn=="CYQ"||resn=="CEA"||resn=="CYG"||resn=="BUC"||resn=="PEC"||
            resn=="CYM"||resn=="CY3"||resn=="CSO"||resn=="SOC"||resn=="CSX"||
            resn=="CSW"||resn=="EFC"||resn=="CSP"||resn=="CSS"||resn=="SCH"||
            resn=="OCS"||resn=="SHC"||resn=="C5C") return 'C';
        if (resn=="DGL"||resn=="GGL"||resn=="CGU"||resn=="GMA"||resn=="5HP"||
            resn=="PCA") return 'E';
        if (resn=="ASQ"||resn=="ASB"||resn=="ASA"||resn=="ASK"||resn=="ASL"||
            resn=="2AS"||resn=="DAS"||resn=="DSP"||resn=="BHD") return 'D';
        if (resn=="PHI"||resn=="PHL"||resn=="DPN"||resn=="DAH"||resn=="HPQ")
            return 'F';
        if (resn=="GLZ"||resn=="SAR"||resn=="GSC"||resn=="GL3"||resn=="MSA"||
            resn=="MPQ"||resn=="NMC") return 'G';
        if (resn=="NEM"||resn=="NEP"||resn=="HSD"||resn=="HSP"||resn=="MHS"||
            resn=="3AH"||resn=="HIC"||resn=="HIP"||resn=="DHI"||resn=="HSE") 
            return 'H';
        if (resn=="IIL"||resn=="DIL") return 'I';
        if (resn=="DLY"||resn=="LYZ"||resn=="SHR"||resn=="ALY"||resn=="TRG"||
            resn=="LYM"||resn=="LLY"||resn=="KCX") return 'K';
        if (resn=="NLE"||resn=="CLE"||resn=="NLP"||resn=="DLE"||resn=="BUG"||
            resn=="NLN"||resn=="MLE") return 'L';
        if (resn=="FME"||resn=="CXM"||resn=="OMT") return 'M';
        if (resn=="MEN") return 'N';
        if (resn=="DPR"||resn=="HYP") return 'P';
        if (resn=="DGN") return 'Q';
        if (resn=="AGM"||resn=="ACL"||resn=="DAR"||resn=="HAR"||resn=="HMR"||
            resn=="ARM") return 'R';
        if (resn=="OAS"||resn=="MIS"||resn=="SAC"||resn=="SEL"||resn=="SVA"||
            resn=="SET"||resn=="DSN"||resn=="SEP") return 'S';
        if (resn=="DTH"||resn=="TPO"||resn=="ALO"||resn=="BMT") return 'T';
        if (resn=="DVA"||resn=="MVA"||resn=="DIV") return 'V';
        if (resn=="LTR"||resn=="DTR"||resn=="TRO"||resn=="TPL"||resn=="HTR") 
            return 'W';
        if (resn=="PAQ"||resn=="STY"||resn=="TYQ"||resn=="IYR"||resn=="TYY"||
            resn=="DTY"||resn=="TYB"||resn=="PTR"||resn=="TYS") return 'Y';
        
        // undeterminted amino acid
        if (resn=="ASX") return 'B'; // or D or N
        if (resn=="GLX") return 'Z'; // or Q or E

        // amino acid with no code in BLOSUM62
        if (convertX>=3)
        {
            if (resn=="SEC") return 'U';
            if (resn=="PYL") return 'O';
        }
        if (resn=="SEC") return 'C';
        if (resn=="PYL") return 'K';
    }
    return 'X';
}

/* only residues in 'ATOM' record with CA or C3' atoms are converted */
string pdb2fasta(ChainUnit& chain)
{
    chain.sequence="";
    size_t r,a;
    for (r=0;r<chain.residues.size();r++)
    {
        if (chain.residues[r].het==false)
        {
            //for (a=0;a<chain.residues[r].atoms.size();a++)
            //{
                //if (chain.residues[r].atoms[a].name==" CA "||
                    //chain.residues[r].atoms[a].name==" C3'")
                    chain.sequence+=aa3to1(chain.residues[r].resn);
            //}
        }
    }
    return chain.sequence;
}

/* ShowSeqLen - whether to show residue number for each chain */
string pdb2fasta(ModelUnit& pep,const string PDBid="",const int ShowSeqLen=0)
{
    stringstream buf;
    string sequence="";
    for (size_t c=0;c<pep.chains.size();c++)
    {
        sequence=pdb2fasta(pep.chains[c]);
        buf<<'>'<<PDBid<<':'<<pep.chains[c].chainID_full;
        if (ShowSeqLen) buf<<'\t'<<sequence.length();
        buf<<'\n'<<sequence<<'\n';
    }
    sequence.clear();
    return buf.str();
}

/* remove sidechain or backbone atoms 
 * atomic_detail - 1: only remove sidechain atoms
 *                 0: remove all non-CA atom*/
void remove_sidechain(ResidueUnit& residue,int atomic_detail=1)
{
    vector<AtomUnit> atoms; // list of atoms
    for (size_t a=0;a<residue.atoms.size();a++)
    {
        if ((atomic_detail==0 && residue.atoms[a].name==" CA ")||
            (atomic_detail==1 &&(residue.atoms[a].name==" CA " ||
             residue.atoms[a].name==" N  " || residue.atoms[a].name==" C  "
          || residue.atoms[a].name==" O  ")))
            atoms.push_back(residue.atoms[a]);
    }
    residue.atoms=atoms;
    atoms.clear();
}

void remove_sidechain(ChainUnit& chain,int atomic_detail=1)
{
    for (size_t r=0;r<chain.residues.size();r++)
        remove_sidechain(chain.residues[r],atomic_detail);
}

void remove_sidechain(ModelUnit& pep,int atomic_detail=1)
{
    for (size_t c=0;c<pep.chains.size();c++)
        remove_sidechain(pep.chains[c],atomic_detail);
}
#endif
