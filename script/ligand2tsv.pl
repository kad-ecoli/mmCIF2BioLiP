#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

my $indir  ="$rootdir/pdb/refdata/chem_comp";
my $outfile="$rootdir/data/ligand.tsv";

my $txt="#CCD\tformula\tInChI\tInChIKey\tSMILES\tname\n";

foreach my $suffix(`ls $indir`)
{
    chomp($suffix);
    print "parsing ligand *$suffix\n";
    foreach my $res(`ls $indir/$suffix`)
    {
        chomp($res);
        my $filename="$indir/$suffix/$res/$res.cif";
        print "$filename\n";
        
        my $_chem_comp_id="";
        my $_chem_comp_name="";
        my $_chem_comp_formula="";
        my $InChI="";
        my $InChIKey="";
        my @SMILES_list=();
        my $_pdbx_chem_comp_descriptor=0;

        my @lines=`cat $filename`;
        for (my $l=0;$l<scalar @lines;$l++)
        {
            my $line=$lines[$l];
            chomp($line);
            if ($line=~/^#/)
            {
                $_pdbx_chem_comp_descriptor=0;
            }
            if ($line=~/^_chem_comp.id\s+/)
            {
                $_chem_comp_id=&strip(substr($line,14));
            }
            elsif ($line=~/^_chem_comp.formula\s+/)
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
            }
        }
        $_chem_comp_name=~s/\?; //g;
        $_chem_comp_name=~s/; $//;
   
        my %SMILES_hash = map { $_, 1 } @SMILES_list;
        my @SMILES_list = keys %SMILES_hash;
        my $SMILES = join("; ",@SMILES_list);
        $txt.="$_chem_comp_id\t$_chem_comp_formula\t$InChI\t";
        $txt.="$InChIKey\t$SMILES\t$_chem_comp_name\n";
    }
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
