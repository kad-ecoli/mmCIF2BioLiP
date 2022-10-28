#!/usr/bin/perl
# download uniprot accession, taxonomy, EC, GO and pubmed from sifts
# download ncbi taxonomy and scientific name from ncbi
# download ccd definition to ligand from pdb
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download CCD ligand\n";

system("mkdir -p $rootdir/pdb/data/monomers/");
my $infile ="$rootdir/pdb/data/monomers/components.cif.gz";
#system("wget https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz -O $infile");
if (!-s "$infile")
{
    print "ERROR! cannot download $infile\n";
    exit(1);
}
my $outfile="$rootdir/data/ligand.tsv";

my $txt="#CCD\tformula\tInChI\tInChIKey\tSMILES\tname\n";
my $_chem_comp_id="";
my $_chem_comp_name="";
my $_chem_comp_formula="";
my $InChI="";
my $InChIKey="";
my @SMILES_list=();
my $_pdbx_chem_comp_descriptor=0;

my @lines=`zcat $infile`;
for (my $l=0;$l<scalar @lines;$l++)
{
    my $line=$lines[$l];
    chomp($line);
    if ($l+1<scalar @lines && $lines[1+$l]=~/^;/)
    {
        $l++;
        $line.=" ".substr($lines[$l],1);
        chomp($line);
        while ($l+1<scalar @lines && $lines[1+$l]!~/^;/)
        {
            $l++;
            $line.=" ".$lines[$l];
            chomp($line);
        }
    }

    if ($line=~/^data_/)
    {
        if (length $_chem_comp_id)
        {
            $_chem_comp_name=~s/\?; //g;
            $_chem_comp_name=~s/; $//;
   
            my %SMILES_hash = map { $_, 1 } @SMILES_list;
            my $SMILES = join("; ",keys %SMILES_hash);
            $txt.="$_chem_comp_id\t$_chem_comp_formula\t$InChI\t";
            $txt.="$InChIKey\t$SMILES\t$_chem_comp_name\n";
        }

        $_chem_comp_id="";
        $_chem_comp_name="";
        $_chem_comp_formula="";
        $InChI="";
        $InChIKey="";
        @SMILES_list=();
        $_pdbx_chem_comp_descriptor=0;
    }
    elsif ($line=~/^#/)
    {
        $_pdbx_chem_comp_descriptor=0;
    }
    if ($line=~/^_chem_comp.id\b/)
    {
        $_chem_comp_id=&strip(substr($line,14));
    }
    elsif ($line=~/^_chem_comp.formula\b/)
    {
        $_chem_comp_formula=&strip(substr($line,19));
    }
    elsif ($line=~/^_chem_comp.name\s+/)
    {
        $_chem_comp_name.=&strip(substr($line,16))."; ";
    }
    elsif ($line=~/^_chem_comp.pdbx_synonyms\s+/)
    {
        $_chem_comp_name.=&strip(substr($line,25))."; ";
    }
    elsif ($line=~/^_pdbx_chem_comp_descriptor.descriptor/)
    {
        $_pdbx_chem_comp_descriptor=1;
    }
    elsif ($_pdbx_chem_comp_descriptor==1)
    {
        if ($line=~/^$_chem_comp_id\s+/ && $lines[$l+1]=~/^"/)
        {
            $l++;
            $line.=" $lines[$l]";
        }

        if ($line=~/^$_chem_comp_id\s+InChI\s+/ &&
               $line=~/(\S+)\s*$/)
        {
            $InChI=&strip("$1");
        }
        elsif ($line=~/^$_chem_comp_id\s+InChIKey\s+/ &&
               $line=~/(\S+)\s*$/)
        {
            $InChIKey=&strip("$1");
        }
        elsif ($line=~/^$_chem_comp_id\s+SMILES\s+/ &&
               $line=~/(\S+)\s*$/)
        {
            push(@SMILES_list,&strip("$1"));
        }
        elsif ($line=~/^$_chem_comp_id\s+SMILES_CANONICAL\s+/ &&
               $line=~/(\S+)\s*$/)
        {
            push(@SMILES_list,&strip("$1"));
        }
    }
}

if (length $_chem_comp_id)
{
    $_chem_comp_name=~s/\?; //g;
    $_chem_comp_name=~s/; $//;
   
    my %SMILES_hash = map { $_, 1 } @SMILES_list;
    my $SMILES = join("; ",keys %SMILES_hash);
    $txt.="$_chem_comp_id\t$_chem_comp_formula\t$InChI\t";
    $txt.="$InChIKey\t$SMILES\t$_chem_comp_name\n";
}

open(FP,">$outfile");
print FP $txt;
close(FP);
system("gzip -f $outfile");
exit(0);

sub strip
{
    my ($instring)=@_;
    $instring=~s/^\s*//;
    $instring=~s/\s*$//;
    $instring=~s/^"//;
    $instring=~s/"$//;
    return $instring;
}
