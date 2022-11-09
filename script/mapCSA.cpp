/* Compile this program by: 
 * $ g++ -O3 -std=c++11 fasta2nr.cpp -o fasta2nr
 */

const char* docstring=""
"mapCSA csa.tsv weekly/\n"
"    read m-csa catalytic sites at csa.tsv\n"
"    map the residue number to pdb files at weekly/*/{receptor,receptor_nr}\n" 
;

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
using namespace std;

/* StringTools START */

/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void Split(const string &line, vector<string> &line_vec,
    const char delimiter=' ')
{
    bool within_word = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (line[pos]==delimiter)
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

void clear_line_vec(vector<string> &line_vec)
{
    int i;
    for (i=0;i<line_vec.size();i++) line_vec[i].clear();
    line_vec.clear();
}

inline bool StartsWith(const string &longString, const string &shortString)
{
    return (longString.size()>=shortString.size() &&
            longString.substr(0,shortString.size())==shortString);
}


inline bool EndsWith(const string &longString, const string &shortString)
{
    return (longString.size()>=shortString.size() &&
            longString.substr(longString.size()-shortString.size(),
                shortString.size())==shortString);
}

string Trim(const string &inputString,const string &char_list=" ")
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(char_list);
    int idxEnd = inputString.find_last_not_of(char_list);
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    else result = "";
    return result;
}

/* StringTools END */
/* main START */

bool isfile(const string&filename)
{
    ifstream fp(filename.c_str());
    if (fp.good())
    {
        fp.close();
        return true;
    }
    return false;
}

inline char aa3to1(const string resn)
{
    if (resn.size()==1)      return tolower(resn[0]);
    else if (resn.size()==2) return tolower(resn[1]);
    else if (resn[0]==' ')   return tolower(resn[2]);
    else if (resn=="PSU")    return 'u';

    // 20 standard amino acid + MSE
    else if (resn=="ALA") return 'A';
    else if (resn=="CYS") return 'C';
    else if (resn=="ASP") return 'D';
    else if (resn=="GLU") return 'E';
    else if (resn=="PHE") return 'F';
    else if (resn=="GLY") return 'G';
    else if (resn=="HIS") return 'H';
    else if (resn=="ILE") return 'I';
    else if (resn=="LYS") return 'K';
    else if (resn=="LEU") return 'L';
    else if (resn=="MET") return 'M';
    else if (resn=="ASN") return 'N';
    else if (resn=="PRO") return 'P';
    else if (resn=="GLN") return 'Q';
    else if (resn=="ARG") return 'R';
    else if (resn=="SER") return 'S';
    else if (resn=="THR") return 'T';
    else if (resn=="VAL") return 'V'; 
    else if (resn=="TRP") return 'W';
    else if (resn=="TYR") return 'Y';

    if (resn=="MSE") return 'M';

    // non-standard amino acid with known parent
    if (resn=="CHG"||resn=="HAC"||resn=="AYA"||resn=="TIH"||resn=="BNN"||
        resn=="ALM"||resn=="TPQ"||resn=="MAA"||resn=="PRR"||resn=="FLA"||
        resn=="AIB"||resn=="DAL"||resn=="CSD"||resn=="DHA"||resn=="DNP") 
        return 'A';
    else if (resn=="PR3"||resn=="CCS"||resn=="C6C"||resn=="SMC"||resn=="BCS"||
             resn=="SCY"||resn=="DCY"||resn=="SCS"||resn=="CME"||resn=="CY1"||
             resn=="CYQ"||resn=="CEA"||resn=="CYG"||resn=="BUC"||resn=="PEC"||
             resn=="CYM"||resn=="CY3"||resn=="CSO"||resn=="SOC"||resn=="CSX"||
             resn=="CSW"||resn=="EFC"||resn=="CSP"||resn=="CSS"||resn=="SCH"||
             resn=="OCS"||resn=="SHC"||resn=="C5C") return 'C';
    else if (resn=="DGL"||resn=="GGL"||resn=="CGU"||resn=="GMA"||resn=="5HP"||
             resn=="PCA") return 'E';
    else if (resn=="ASQ"||resn=="ASB"||resn=="ASA"||resn=="ASK"||resn=="ASL"||
             resn=="2AS"||resn=="DAS"||resn=="DSP"||resn=="BHD") return 'D';
    else if (resn=="PHI"||resn=="PHL"||resn=="DPN"||resn=="DAH"||resn=="HPQ")
        return 'F';
    else if (resn=="GLZ"||resn=="SAR"||resn=="GSC"||resn=="GL3"||resn=="MSA"||
             resn=="MPQ"||resn=="NMC") return 'G';
    else if (resn=="NEM"||resn=="NEP"||resn=="HSD"||resn=="HSP"||resn=="MHS"||
             resn=="3AH"||resn=="HIC"||resn=="HIP"||resn=="DHI"||resn=="HSE") 
        return 'H';
    else if (resn=="IIL"||resn=="DIL") return 'I';
    else if (resn=="DLY"||resn=="LYZ"||resn=="SHR"||resn=="ALY"||resn=="TRG"||
             resn=="LYM"||resn=="LLY"||resn=="KCX") return 'K';
    else if (resn=="NLE"||resn=="CLE"||resn=="NLP"||resn=="DLE"||resn=="BUG"||
             resn=="NLN"||resn=="MLE") return 'L';
    else if (resn=="FME"||resn=="CXM"||resn=="OMT") return 'M';
    else if (resn=="MEN") return 'N';
    else if (resn=="DPR"||resn=="HYP") return 'P';
    else if (resn=="DGN") return 'Q';
    else if (resn=="AGM"||resn=="ACL"||resn=="DAR"||resn=="HAR"||resn=="HMR"||
             resn=="ARM") return 'R';
    else if (resn=="OAS"||resn=="MIS"||resn=="SAC"||resn=="SEL"||resn=="SVA"||
             resn=="SET"||resn=="DSN"||resn=="SEP") return 'S';
    else if (resn=="DTH"||resn=="TPO"||resn=="ALO"||resn=="BMT") return 'T';
    else if (resn=="DVA"||resn=="MVA"||resn=="DIV") return 'V';
    else if (resn=="LTR"||resn=="DTR"||resn=="TRO"||resn=="TPL"||resn=="HTR") 
        return 'W';
    else if (resn=="PAQ"||resn=="STY"||resn=="TYQ"||resn=="IYR"||resn=="TYY"||
             resn=="DTY"||resn=="TYB"||resn=="PTR"||resn=="TYS") return 'Y';
    
    // undeterminted amino acid
    else if (resn=="ASX") return 'B'; // or D or N
    else if (resn=="GLX") return 'Z'; // or Q or E
    else if (resn=="SEC") return 'U';
    else if (resn=="PYL") return 'O';
    return 'X';
}

void mapCSA(const string &infile, const string &infolder, const string &outfile)
{
    stringstream buf;
    ifstream fp;
    if (infile=="-") buf<<cin.rdbuf();
    else
    {
        fp.open(infile.c_str(),ios::in); //ifstream fp(filename,ios::in);
        buf<<fp.rdbuf();
        fp.close();
    }
    vector<string> lines;
    Split(buf.str(),lines,'\n');
    buf.str(string());

    string outtxt;
    string pdbid,asym_id;
    string line;
    string filename;
    string csaOrig,csaRenu;
    vector<string> pdblines;
    string resi;
    size_t l,a,r;
    vector<string> line_vec;
    string aa;
    for (l=0;l<lines.size();l++)
    {
        Split(lines[l],line_vec,'\t');
        if (line_vec.size()<3)
        {
            cerr<<"WARNING! Skip line "<<l+1<<"\t"<<lines[l]<<endl;
            clear_line_vec(line_vec);
            continue;
        }
        pdbid=line_vec[0];
        asym_id=line_vec[1];
        line=line_vec[2];
        clear_line_vec(line_vec);
        filename=infolder+pdbid.substr(pdbid.size()-3,2)+
            "/receptor/"+pdbid+asym_id+".pdb";
        if (!isfile(filename))
        {
            filename=infolder+pdbid.substr(pdbid.size()-3,2)+
                "/receptor_nr/"+pdbid+asym_id+".pdb";
            if (!isfile(filename))
            {
                //cerr<<"WARNING! No such file "<<infolder
                    //<<pdbid.substr(pdbid.size()-3,2)
                    //<<"/receptor*/"<<pdbid<<asym_id<<".pdb"<<endl;
                continue;
            }
        }
        Split(line,line_vec,',');
        csaOrig="";
        csaRenu="";
        
        fp.open(filename.c_str(),ios::in);
        buf<<fp.rdbuf();
        fp.close();
        Split(buf.str(),pdblines,'\n');
        buf.str(string());

        r=0; // residue number reindexed from 1
        for (a=0;a<pdblines.size();a++)
        {
            if ((!StartsWith(pdblines[a],"ATOM  ") &&
                !StartsWith(pdblines[a],"HETATM"))||
                pdblines[a].size()<54 || pdblines[a].substr(12,4)!=" CA ")
                continue;
            r++;
            resi=Trim(pdblines[a].substr(22,5));
            if (find(line_vec.begin(),line_vec.end(),resi)==line_vec.end())
                continue;
            aa=aa3to1(pdblines[a].substr(17,3));
            buf<<' '<<aa<<resi;
            csaOrig+=buf.str();
            buf.str(string());
            buf<<' '<<aa<<r;
            csaRenu+=buf.str();
            buf.str(string());
        }
        if (csaOrig.size())
            outtxt+=pdbid+'\t'+asym_id+'\t'+csaOrig.substr(1)+
                                       '\t'+csaRenu.substr(1)+'\n';
        //else cerr<<"WARNING! fail to map "<<pdbid<<':'<<asym_id
        //    <<" (L="<<r<<")\t"<<line<<endl;
        clear_line_vec(pdblines);
        clear_line_vec(line_vec);
    }

    if (outfile=="-") cout<<outtxt<<flush;
    else
    {
        ofstream fout;
        fout.open(outfile.c_str());
        fout<<outtxt<<flush;
        fout.close();
    }

    /* clean up */
    vector<string> ().swap(lines);
    string ().swap(outtxt);
    string ().swap(pdbid);
    string ().swap(asym_id);
    string ().swap(line);
    string ().swap(filename);
    string ().swap(csaOrig);
    string ().swap(csaRenu);
    vector<string> ().swap(pdblines);
    vector<string> ().swap(line_vec);
    return;
}

int main(int argc,char **argv)
{
    string infile  ="";
    string infolder="";
    string outfile ="";

    for (int a=1;a<argc;a++)
    {
        if (infile.size()==0)
            infile=argv[a];
        else if (infolder.size()==0)
            infolder=argv[a];
        else if (outfile.size()==0)
            outfile=argv[a];
        else
        {
            cerr<<"ERROR: unknown option "<<argv[a]<<endl;
            return 1;
        }
    }

    if (infolder.size()==0)
    {
        cerr<<docstring;
        return 1;
    }
    if (outfile=="") outfile="-";

    if (!EndsWith(infolder,"/") && !EndsWith(infolder,"\\"))
        infolder+="/";

    mapCSA(infile,infolder,outfile);

    /* clean up */
    string ().swap(infile);
    string ().swap(infolder);
    string ().swap(outfile);
    return 0;
}

/* main END */
