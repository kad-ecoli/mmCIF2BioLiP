#!/usr/bin/perl
# download enzyme definition
# download ccd definition to ligand from pdb
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download enzyme\n";
system("mkdir -p $rootdir/enzyme/") if (!-d "$rootdir/enzyme/");
system("wget -q http://ftp.expasy.org/databases/enzyme/enzyme.dat -O $rootdir/enzyme/enzyme.dat");
system("wget -q  ftp://ftp.expasy.org/databases/enzyme/enzyme.dat -O $rootdir/enzyme/enzyme.dat") if (!-s"$rootdir/enzyme/enzyme.dat");
my $ID;
my $DE;
my $txt;
foreach my $line(`cat $rootdir/enzyme/enzyme.dat |grep -P "(^//)|(^ID)|(^DE)"`)
{
    chomp($line);
    if ($line=~/^\/\//)
    {
        $txt.="$ID\t$DE\n" if (length $ID);
        $ID="";
        $DE="";
    }
    elsif ($line=~/^ID/)
    {
        $ID=substr($line,5);
    }
    elsif ($line=~/^DE/)
    {
        if (length $DE==0) { $DE=substr($line,5); }
        else          { $DE.=" ".substr($line,5); }
    }
}
open(FP,">$rootdir/data/enzyme.tsv");
print FP $txt;
close(FP);
&gzipFile("$rootdir/data/enzyme.tsv");


print "download unichem\n";
#1	chembl
#2	drugbank
#3	pdb
#7	chebi
#9	zinc
system("mkdir -p $rootdir/UniChem/");
system("wget -q http://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src3.txt.gz -O $rootdir/UniChem/src1src3.txt.gz");
system("wget -q  ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src3.txt.gz -O $rootdir/UniChem/src1src3.txt.gz") if (!-s "$rootdir/UniChem/src1src3.txt.gz");
system("wget -q http://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id2/src2src3.txt.gz -O $rootdir/UniChem/src2src3.txt.gz");
system("wget -q  ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id2/src2src3.txt.gz -O $rootdir/UniChem/src2src3.txt.gz") if (!-s "$rootdir/UniChem/src2src3.txt.gz");
system("wget -q http://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id3/src3src7.txt.gz -O $rootdir/UniChem/src3src7.txt.gz");
system("wget -q  ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id3/src3src7.txt.gz -O $rootdir/UniChem/src3src7.txt.gz") if (!-s "$rootdir/UniChem/src3src7.txt.gz");
system("wget -q http://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id3/src3src9.txt.gz -O $rootdir/UniChem/src3src9.txt.gz");
system("wget -q  ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id3/src3src9.txt.gz -O $rootdir/UniChem/src3src9.txt.gz") if (!-s "$rootdir/UniChem/src3src9.txt.gz");
my %pdb2chembl;
foreach my $line(`zcat $rootdir/UniChem/src1src3.txt.gz|tail -n +2`)
{ 
    $pdb2chembl{"$2"}="$1" if ($line=~/(\S+)\t(\S+)/);
}
my %pdb2drugbank;
foreach my $line(`zcat $rootdir/UniChem/src2src3.txt.gz|tail -n +2`)
{ 
    $pdb2drugbank{"$2"}="$1" if ($line=~/(\S+)\t(\S+)/);
}
my %pdb2zinc;
foreach my $line(`zcat $rootdir/UniChem/src3src9.txt.gz|tail -n +2`)
{ 
    $pdb2zinc{"$1"}="$2" if ($line=~/(\S+)\t(\S+)/);
}


print "download CCD ligand\n";
system("mkdir -p $rootdir/pdb/data/monomers/");
my $infile ="$rootdir/pdb/data/monomers/components.cif.gz";
system("wget -q http://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz -O $infile");
system("wget -q  ftp://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz -O $infile") if (!-s "$infile");
if (!-s "$infile")
{
    print "ERROR! cannot download $infile\n";
    exit(1);
}
my $outfile="$rootdir/data/ligand.tsv";

$txt="#CCD\tformula\tInChI\tInChIKey\tSMILES\tname\tChEMBL\tDrugBank\tZINC\n";
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
            my $SMILES  = join("; ",keys %SMILES_hash);
            my $chembl  ="";
            my $drugbank="";
            my $zinc    ="";
            $chembl  =$pdb2chembl{$_chem_comp_id}   if (exists $pdb2chembl{$_chem_comp_id});
            $drugbank=$pdb2drugbank{$_chem_comp_id} if (exists $pdb2drugbank{$_chem_comp_id});
            $zinc    =$pdb2zinc{$_chem_comp_id}     if (exists $pdb2zinc{$_chem_comp_id});
            $txt.="$_chem_comp_id\t$_chem_comp_formula\t$InChI\t$InChIKey\t";
            $txt.="$SMILES\t$_chem_comp_name\t$chembl\t$drugbank\t$zinc\n";
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

    my $chembl  ="";
    my $drugbank="";
    my $zinc    ="";
    $chembl  =$pdb2chembl{$_chem_comp_id}   if (exists $pdb2chembl{$_chem_comp_id});
    $drugbank=$pdb2drugbank{$_chem_comp_id} if (exists $pdb2drugbank{$_chem_comp_id});
    $zinc    =$pdb2zinc{$_chem_comp_id}     if (exists $pdb2zinc{$_chem_comp_id});
    $txt.="$_chem_comp_id\t$_chem_comp_formula\t$InChI\t$InChIKey\t";
    $txt.="$SMILES\t$_chem_comp_name\t$chembl\t$drugbank\t$zinc\n";
}

open(FP,">$outfile");
print FP $txt;
close(FP);
system("cat $outfile |grep -P '^\\w+\\t[A-Z][a-z]{0,1}\\tInChI=1S\\/[A-Z][a-z]{0,1}/q\\+\\d' > $rootdir/data/metal.tsv");
&gzipFile("$outfile");
&gzipFile("$rootdir/data/metal.tsv");
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

sub gzipFile
{
    my ($filename)=@_;
    my $oldNum=`zcat $filename.gz 2>/dev/null|wc -l`+0;
    my $newNum=` cat $filename   |wc -l`+0;
    if (0.8*$oldNum>$newNum)
    {
        print "WARNING! do not update $filename from $oldNum to $newNum entries\n";
        return;
    }
    print "update $filename from $oldNum to $newNum entries\n";
    system("gzip -f $filename");
    return;
}
