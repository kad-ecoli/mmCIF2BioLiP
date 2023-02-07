#include "CSSR.h"
#include "PDBParser.h"
#include "cssr_struct.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
CSSR::CSSR() {

	// Initialize the calculation type description.
	calcType = "Calculation of probable structures by Thermodynamics parameters";

	// Initialize the nucleic acid type.
	isRNA = true;

    isSequence = true; // dummy variable
    fastOpt    = 2;

	// default output format is DSSR
    show_dot   = false;
    show_dssr  = true;
    show_conf  = false;
    show_bpseq = false;
	
    // whether to analysis inter-chain paring only. default is false
    interchain = false;

    // default is read all atoms
    atom = "auto";

	// Initialize the probable pairing threshold.
	threshold = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool CSSR::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "CSSR" );
	parser->addParameterDescription( "input file", "The name of the input PDB file." );
	parser->addParameterDescription( "output file", "The name of the output file. If '-', write to stdout" );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the molecule is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );
	
    // Add the DNA option.
	vector<string> interchainOptions;
	interchainOptions.push_back( "-c" );
	interchainOptions.push_back( "--interchain" );
	parser->addOptionFlagsNoParameters( interchainOptions, "Specify that only inter-chain base pairs are reported. Default is to include both intra- and inter-chain base pairs." );

	// Add the fast option.
	vector<string> fastOptions;
	fastOptions.push_back( "-f" );
	fastOptions.push_back( "--fast" );
	parser->addOptionFlagsWithParameters( fastOptions, "Specify fast level: 0 - (slow) thermodynamics with quadratic folding; 1 - thermodynamics with linear folding; 2 - (default) do not use thermodyanmics parameters." );

	// Add the outfmt option.
	vector<string> outfmtOptions;
	outfmtOptions.push_back( "-o" );
	outfmtOptions.push_back( "--outfmt" );
	parser->addOptionFlagsWithParameters( outfmtOptions, "The outut format, which can be added together. Valid formats are: 1 - dot bracket format. 2 - (default) DSSR format. 4 - confidence score. 8 - BPseq (single chain only)");
	
    // Add the atom option.
	vector<string> atomOptions;
	atomOptions.push_back( "-a" );
	atomOptions.push_back( "--atom" );
	parser->addOptionFlagsWithParameters( atomOptions, "atom used for assignment. default is all the following atoms: \" P  \", \" O5'\", \" O4'\", \" O3'\", \" C5'\", \" C4'\", \" C3'\", \" C2'\", \" C1'\" and N (N9 for a/g; N1 for c/t/u)\n");

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		input = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

    interchain = parser->contains( interchainOptions );

	// Get the sequence flag.
	if( !parser->isError() )
        parser->setOptionInteger( fastOptions, fastOpt );

	// Get the DNA option.
    isRNA = !parser->contains( dnaOptions );

    // Get the outfmt option.
	if( !parser->isError() ) {
        int outfmt=2;
		parser->setOptionInteger( outfmtOptions, outfmt );
        show_dot    =(outfmt%2==1); outfmt/=2;
        show_dssr   =(outfmt%2==1); outfmt/=2;
        show_conf   =(outfmt%2==1); outfmt/=2;
        show_bpseq  =(outfmt%2==1);
	}
    
    // Get the atom option.
	if( !parser->isError() ) {
		parser->setOptionString( atomOptions, atom );
        if (atom.size()==1) atom+=" ";
        if (atom.size()==2) atom+=" ";
        if (atom.size()==3) atom=" "+atom;
        if (atom.size()!=4) {
		    parser->setError( "atom should be 4-character long" ); }
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////

bool compareRR (pair<float,pair<int,int> > resi_pair1,
    pair<float,pair<int,int> > resi_pair2)
{
    return (resi_pair1.first>=resi_pair2.first);
}

void CSSR::run() {
    // parse 3D structure
    int atomic_detail=2;
    int allowX       =1; // only allow ATOM and MSE
    ios_base::sync_with_stdio(false);
    ModelUnit pdb_entry=read_pdb_structure(input.c_str(),atomic_detail,allowX);
    if (atom!="auto") removeUnspecifiedAtom(pdb_entry, atom);
    
    // parse sequence into upper case
    string sequence="";
    size_t c,r;
    const int convertX=2;
    for (size_t c=0;c<pdb_entry.chains.size();c++)
        for (size_t r=0;r<pdb_entry.chains[c].residues.size();r++)
            if (pdb_entry.chains[c].residues[r].het==false)
                sequence+=toupper(aa3to1(pdb_entry.chains[c].residues[r].resn,convertX));

    // ProbablePairRR
    vector<pair<float,pair<int,int> > >RR_list;
    
    // assign SS solely by structure
    vector<string> res_str_vec;
    vector<pair<float,vector<string> > > bp_vec;
    cssr(pdb_entry, res_str_vec, bp_vec, fastOpt, RR_list, interchain);
    size_t bp;
    vector<size_t> filtered_bp_vec;
    if (show_dot || show_dssr || show_bpseq)
        filter_bp(bp_vec, filtered_bp_vec);

    ofstream fp;
    bool usestdout=(ctFile=="-");
    if (!usestdout) fp.open(ctFile.c_str());
    if (show_dot)
    {
        vector<char> dot_bracket;
        cssr_bp2dot(res_str_vec, filtered_bp_vec, bp_vec, dot_bracket);
        if (usestdout)
        {
            for (size_t r=0;r<dot_bracket.size();r++) cout<<dot_bracket[r];
            cout<<endl;
        }
        else
        {
            for (size_t r=0;r<dot_bracket.size();r++) fp<<dot_bracket[r];
            fp<<endl;
        }
        dot_bracket.clear();
    }
    if (show_dssr)
    {
        size_t count=0;
        if (usestdout) cout<<"\n"
            <<"****************************************************************************"
            <<"\nList of "<<filtered_bp_vec.size()<<" base pairs\n"
            <<"     nt1            nt2            bp  name        Saenger   LW   DSSR"<<endl;
        else fp<<"\n"
            <<"****************************************************************************"
            <<"\nList of "<<filtered_bp_vec.size()<<" base pairs\n"
            <<"     nt1            nt2            bp  name        Saenger   LW   DSSR"<<endl;
        for (bp=0;bp<bp_vec.size();bp++)
        {
            if (find(filtered_bp_vec.begin(), filtered_bp_vec.end(),bp
                )==filtered_bp_vec.end()) continue;
            if (usestdout)
            {
                cout<<setw(4)<<right<<++count<<' '
                    <<setw(14)<<left<<bp_vec[bp].second[0]<<' '
                    <<setw(14)<<left<<bp_vec[bp].second[1]<<' '
                    <<bp_vec[bp].second[2];
                if (bp_vec[bp].second[2].substr(4)=="WC")
                    cout<<"          19-XIX    cWW  cW-W";
                else if (bp_vec[bp].second[2].substr(4)=="Wobble")
                    cout<<"      28-XXVIII cWW  cW-W";
                cout<<endl;
            }
            else
            {
                fp  <<setw(4)<<right<<++count<<' '
                    <<setw(14)<<left<<bp_vec[bp].second[0]<<' '
                    <<setw(14)<<left<<bp_vec[bp].second[1]<<' '
                    <<bp_vec[bp].second[2];
                if (bp_vec[bp].second[2].substr(4)=="WC")
                    fp  <<"          19-XIX    cWW  cW-W";
                else if (bp_vec[bp].second[2].substr(4)=="Wobble")
                    fp  <<"      28-XXVIII cWW  cW-W";
                fp  <<endl;
            }
        }
    }
    if (show_conf && bp_vec.size())
    {
        sort(bp_vec.begin(), bp_vec.end());
        for (bp=bp_vec.size()-1;bp>0;bp--)
        {
            if (usestdout) cout<<bp_vec[bp].second[0]<<'\t'
                <<bp_vec[bp].second[1]<<'\t'<<bp_vec[bp].first<<endl;
            else           fp  <<bp_vec[bp].second[0]<<'\t'
                <<bp_vec[bp].second[1]<<'\t'<<bp_vec[bp].first<<endl;
        }
    }
    if (show_bpseq)
    {
        vector<vector<int> > bpseq_vec;
        cssr_bpseq(res_str_vec, filtered_bp_vec, bp_vec, bpseq_vec);
        for (size_t r=0;r<bpseq_vec.size();r++)
        {
            if (usestdout) cout<<setw(5)<<bpseq_vec[r][0]<<' '
                <<(char)(bpseq_vec[r][1])<<' '<<setw(5)<<bpseq_vec[r][2]<<endl;
            else           fp  <<setw(5)<<bpseq_vec[r][0]<<' '
                <<(char)(bpseq_vec[r][1])<<' '<<setw(5)<<bpseq_vec[r][2]<<endl;
        }
        vector<vector<int> >().swap(bpseq_vec);
    }
    if (!usestdout) fp.close();

	/* clean up */
    vector<pair<float,pair<int,int> > >().swap(RR_list);
    vector<size_t>().swap(filtered_bp_vec);
    vector<string>().swap(res_str_vec);
    vector<pair<float,vector<string> > >().swap(bp_vec);
    vector<ChainUnit>().swap(pdb_entry.chains);
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	CSSR* runner = new CSSR();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
