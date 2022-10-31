const char* docstring="\n"
"backpropagate input.txt is_a.csv alt_id.csv output.txt\n"
"    back-propagate parent GO terms for the list of input GO terms\n"
"\n"
"Input:\n"
"    input.txt  - list of input GO terms. GO terms separated by comma in\n"
"                 the same line will be back-propagate together. GO terms\n"
"                 at differnet lines will be be back-propagated separately.\n"
"    is_a.csv   - GO_ID Aspect is_a_direct is_a_indirect\n"
"    alt_id.csv - alt_id GO_ID\n"
"\n"
"Output:\n"
"    output.txt - list of input GO terms and their parent terms\n"
;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <malloc.h>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <iomanip>
#include <map>

using namespace std;

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

/* strip white space at the begining or end of string */
string Trim(const string &inputString)
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(" \n\r\t");
    int idxEnd = inputString.find_last_not_of(" \n\r\t");
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    return result;
}

string Join(const string &delimiter, const vector<string>&line_vec)
{
    string line="";
    size_t i;
    for (i=0;i<line_vec.size();i++)
    {
        if (i==0) line+=line_vec[0];
        else      line+=delimiter+line_vec[i];
    }
    return line;
}

void parse_alt_id_file(const string &alt_id_filename, 
    map<string,string> &alt_id_dict)
{
    bool fromStdin=(alt_id_filename=="-");
    ifstream fin;
    vector<string> line_vec;
    string line;
    if (!fromStdin) fin.open(alt_id_filename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        Split(line, line_vec, '\t');
        alt_id_dict[line_vec[0]]=line_vec[1];
        line_vec[0].clear(); line_vec[1].clear(); line_vec.clear();
    }
    if (fromStdin) fin.close();
    vector<string> ().swap(line_vec);
    string ().swap(line);
}

void parse_is_a_file(const string &is_a_filename, 
    map<string,vector<string> >&is_a_dict)
{
    bool fromStdin=(is_a_filename=="-");
    ifstream fin;
    vector<string> line_vec;
    string line;
    string GOterm,Aspect,is_a_direct,is_a_indirect;
    int i;
    if (!fromStdin) fin.open(is_a_filename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        Split(line, line_vec, '\t');
        GOterm       =line_vec[0];
        Aspect       =line_vec[1];
        is_a_direct  =line_vec[2];
        is_a_indirect=line_vec[3];
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        Split(is_a_direct,line_vec,',');
        is_a_dict[GOterm]=line_vec;
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        Split(is_a_indirect,line_vec,',');
        for (i=0;i<line_vec.size();i++)
        {
            is_a_dict[GOterm].push_back(line_vec[i]);
            line_vec[i].clear();
        }
        line_vec.clear();
        GOterm.clear();
        Aspect.clear();
        is_a_direct.clear();
        is_a_indirect.clear();
        is_a_dict[GOterm].push_back(GOterm);
    }
    if (fromStdin) fin.close();
    vector<string> ().swap(line_vec);
    string ().swap(line);
    string ().swap(GOterm);
    string ().swap(Aspect);
    string ().swap(is_a_direct);
    string ().swap(is_a_indirect);
}

void backpropagate(const string &inputfilename,
    map<string,string> &alt_id_dict, map<string,vector<string> > &is_a_dict,
    const string &outputfilename)
{
    string txt;
    vector<string> line_vec;
    vector<string> GOterm_list;
    vector<string> unknown_list;
    string line;
    string GOterm,parent;
    ifstream fin;
    ofstream fout;
    int i,j;
    bool fromStdin=(inputfilename=="-");
    if (!fromStdin) fin.open(inputfilename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        Trim(line);
        if (line.size()==0) continue;
        Split(line, line_vec, ',');
        for (i=0;i<line_vec.size();i++)
        {
            GOterm=line_vec[i];
            if (alt_id_dict.count(GOterm)==0)
            {
                if (find(unknown_list.begin(),unknown_list.end(),GOterm
                     )==unknown_list.end())
                {
                    cerr<<"Skip unknown term '"<<GOterm<<"'"<<endl;
                    unknown_list.push_back(GOterm);
                }
                continue;
            }
            GOterm=alt_id_dict[GOterm];
            if (find(GOterm_list.begin(),GOterm_list.end(),GOterm)!=
                     GOterm_list.end()) continue;
            GOterm_list.push_back(GOterm);
            for (j=0;j<is_a_dict[GOterm].size();j++)
            {
                parent=is_a_dict[GOterm][j];
                if (find(GOterm_list.begin(),GOterm_list.end(),parent
                    )!=GOterm_list.end()) continue;
                GOterm_list.push_back(parent);
            }
            line_vec[i].clear();
        }
        line_vec.clear();
        line.clear();
        txt+=Join(",",GOterm_list)+'\n';
        for (i=0;i<GOterm_list.size();i++) GOterm_list[i].clear();
        GOterm_list.clear();
    }
    if (fromStdin) fin.close();
    
    if (outputfilename=="-") cout<<txt<<flush;
    else
    {
        fout.open(outputfilename.c_str());
        fout<<txt<<flush;
        fout.close();
    }
    
    string ().swap(txt);
    vector<string> ().swap(line_vec);
    vector<string> ().swap(GOterm_list);
    vector<string> ().swap(unknown_list);
    string ().swap(line);
    string ().swap(GOterm);
    string ().swap(parent);
}

int main(int argc,char **argv)
{
    if (argc!=5)
    {
        cerr<<docstring;
        return 0;
    }

    string inputfilename  =argv[1];
    string is_a_filename  =argv[2];
    string alt_id_filename=argv[3];
    string outputfilename =argv[4];
    map<string,string> alt_id_dict;
    map<string,vector<string> > is_a_dict;

    parse_alt_id_file(alt_id_filename, alt_id_dict);
    parse_is_a_file(is_a_filename, is_a_dict);
    backpropagate(inputfilename, alt_id_dict, is_a_dict, outputfilename);
    
    string ().swap(inputfilename);
    string ().swap(is_a_filename);
    string ().swap(alt_id_filename);
    string ().swap(outputfilename);
    map<string,string> ().swap(alt_id_dict);
    map<string,vector<string> > ().swap(is_a_dict);
    return 0;
}
