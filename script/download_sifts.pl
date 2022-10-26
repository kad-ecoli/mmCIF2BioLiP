#!/usr/bin/perl
# download uniprot accession, taxonomy, EC, GO and pubmed from sifts
# download ncbi taxonomy and scientific name from ncbi
# download ccd definition to ligand from pdb
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download sifts\n";
system("mkdir -p $rootdir/sifts/flatfiles/tsv");
foreach my $filename(("pdb_chain_uniprot.tsv.gz",
                      "pdb_chain_taxonomy.tsv.gz",
                      "pdb_pubmed.tsv.gz",
                      "pdb_chain_enzyme.tsv.gz",
                      "pdb_chain_go.tsv.gz"))
{
    system("wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/$filename -O $rootdir/sifts/flatfiles/tsv/$filename");
}

system("mkdir -p $rootdir/taxonomy");
system("wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O $rootdir/taxonomy/taxdump.tar.gz");
system("cd $rootdir/taxonomy; tar -xvf taxdump.tar.gz names.dmp");
if (-s "$rootdir/taxonomy/names.dmp")
{
    system("rm $rootdir/taxonomy/taxdump.tar.gz");
    system("gzip -f $rootdir/taxonomy/names.dmp");
}

system("mkdir -p $rootdir/data");
print "parsing taxonomy\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz")
{
    system("zcat $rootdir/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz |tail -n +3|cut -f1-3|uniq > $rootdir/data/chain2taxonomy.tsv");
}
my @taxid_list;
foreach my $taxid(`cat $rootdir/data/chain2taxonomy.tsv|cut -f3|sort|uniq`)
{
    chomp($taxid);
    push(@taxid_list,($taxid));
}

my %taxid_dict = map { $_ => "" } @taxid_list;
foreach my $line(`zcat $rootdir/taxonomy/names.dmp.gz`)
{
    my @items=split(/\t\|\t/,$line);
    my $taxid=$items[0];
    next if (!exists $taxid_dict{$taxid});
    my $name_type=$items[-1];
    next if (length $taxid_dict{$taxid} && $name_type!~/scientific name/);
    $taxid_dict{$taxid}=$items[1];
}
my $txt="";
foreach my $taxid(@taxid_list)
{
    $txt.="$taxid\t$taxid_dict{$taxid}\n";
}
if (scalar @taxid_list)
{
    open(FP,">$rootdir/data/taxid2name.tsv");
    print FP $txt;
    close(FP);
    system("gzip -f $rootdir/data/taxid2name.tsv");
}
&concatenate_sift("cat $rootdir/data/chain2taxonomy.tsv","$rootdir/data/chain2taxonomy.tsv");

print "parsing uniprot accession\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz")
{
    &concatenate_sift("zcat $rootdir/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz |tail -n +3|cut -f1-3|sort|uniq","$rootdir/data/chain2uniprot.tsv");
}

print "parsing EC\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_chain_enzyme.tsv.gz")
{
    &concatenate_sift("zcat $rootdir/sifts/flatfiles/tsv/pdb_chain_enzyme.tsv.gz |tail -n +3|cut -f1,2,4|sort|uniq","$rootdir/data/chain2ec.tsv");
}

print "parsing GO\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_chain_go.tsv.gz")
{
    &concatenate_sift("zcat $rootdir/sifts/flatfiles/tsv/pdb_chain_go.tsv.gz |tail -n +3|cut -f1,2,6|sort|uniq","$rootdir/data/chain2go.tsv");
}

print "parsing pubmed\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_pubmed.tsv.gz")
{
    system("zcat $rootdir/sifts/flatfiles/tsv/pdb_pubmed.tsv.gz |tail -n +3|grep -P '\\t0\\t'|cut -f1,3|sort|uniq|gzip - > $rootdir/data/pdb2pubmed.tsv.gz");
}

print "download CCD ligand\n";

system("mkdir -p $rootdir/pdb/data/monomers/");
my $infile ="$rootdir/pdb/data/monomers/components.cif.gz";
system("wget https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz -O $infile");
if (!-s "$infile")
{
    print "ERROR! cannot download $infile\n";
    exit(1);
}
my $outfile="$rootdir/data/ligand.tsv";

$txt="#CCD\tformula\tInChI\tInChIKey\tSMILES\tname\n";
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

sub concatenate_sift
{
    my ($cmd,$outfile)=@_;
    my @chain_list;
    my %chain_dict;
    foreach my $line(`$cmd`)
    {
        chomp($line);
        my @items=split("\t",$line);
        my $chain="$items[0]\t$items[1]";
        if (exists $chain_dict{$chain})
        {
            $chain_dict{$chain}.=",$items[2]";
        }
        else
        {
            push(@chain_list,($chain));
            $chain_dict{$chain}="$items[2]";
        }
    }
    my $txt="";
    foreach my $chain(@chain_list)
    {
        $txt.="$chain\t$chain_dict{$chain}\n";
    }
    open(FP,">$outfile");
    print FP $txt;
    close(FP);
    system("gzip -f $outfile");
}
