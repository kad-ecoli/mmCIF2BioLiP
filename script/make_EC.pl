#!/usr/bin/perl
# parse gene ontology
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

### Step 1: create weekly/Enzyme_*.tar.bz2 ###
my @receptor_list;
foreach my $line(`zcat $rootdir/data/pdb_all.tsv.gz |tail -n +2|cut -f1,2|sort|uniq`)
{
    if ($line=~/(\w+)\s(\w+)/)
    {
        my $chain="$1:$2";
        chomp($chain);
        push(@receptor_list,($chain));
    }
}
my %receptor_dict=map { $_, 0 } @receptor_list; 

my @chain2ec_list;
my %chain2ec_dict;
my @miss_list;
foreach my $line(`zcat $rootdir/data/chain2ec.tsv.gz|head -10000`)
{
    if ($line=~/(\w+)\s(\w+)\s(\S+)/)
    {
        my $chain="$1:$2";
        my $ec="$3";
        push(@chain2ec_list,($chain));
        $chain2ec_dict{$chain}=$ec;
        push(@miss_list,($chain)) if (!exists $receptor_dict{$chain});
    }
}
my $size=scalar @miss_list;
print "$size chains with EC but not receptor\n";

my @divided_list;
my %divided_dict;
foreach my $chain(@miss_list)
{
    if ($chain=~/(\w+):(\w+)/)
    {
        my $pdb    ="$1";
        my $asym_id="$2";
        my $divided=substr($pdb,length($pdb)-3,2);
        if (!exists $divided_dict{$divided})
        {
            push(@divided_list,($divided));
            $divided_dict{$divided}=0;
            if (-d "$rootdir/weekly/$divided/Enzyme")
            {
                system("rm -rf $rootdir/weekly/$divided/Enzyme/");
            }
            system("mkdir -p $rootdir/weekly/$divided/Enzyme/");
        }
        $divided_dict{$divided}+=1;
        my $chain="$pdb$asym_id";
        foreach my $suffix(("gz","bz2"))
        {
            next if (-s "$rootdir/weekly/$divided/Enzyme/$chain.pdb");
            next if (!-s "$rootdir/interim/$pdb.tar.$suffix");
            system("cd $rootdir/weekly/$divided/Enzyme; tar -xvf $rootdir/interim/$pdb.tar.$suffix $chain.pdb");
        }
        if (!-s "$rootdir/weekly/$divided/Enzyme/$chain.pdb")
        {
            my $cmd="cd $rootdir/weekly/$divided/Enzyme; $bindir/cif2chain $rootdir/pdb/data/structures/divided/mmCIF/$divided/$pdb.cif.gz $chain.pdb $asym_id";
            #print "$cmd\n";
            system("$cmd");
            if (!-f "$rootdir/weekly/$divided/Enzyme/$chain.pdb")
            {
                print "no such file: $rootdir/weekly/$divided/Enzyme/$chain.pdb\n";
            }
        }
    }
}

foreach my $divided(@divided_list)
{
    if (`ls $rootdir/weekly/$divided/Enzyme|wc -l`+0 >= 1)
    {
        system("cd $rootdir/weekly/$divided/; tar -cjf $rootdir/weekly/Enzyme_$divided.tar.bz2 Enzyme/");
    }
}

### Step 2: create data/ec_all.tsv.gz for all entries with EC number ####
my %resolu_dict;
foreach my $line(`cat $rootdir/pdb/derived_data/index/resolu.idx`)
{
    if ($line=~/^(\S+)\t;\t(\S+)/)
    {
        my $pdbid=lc($1);
        my $resolu="$2";
        $resolu_dict{$pdbid}=$resolu;
    }
}
$size=keys %resolu_dict;
print "$size resolution\n";

my %csa_dict;
foreach my $line(`$bindir/mapCSA $rootdir/data/csa.tsv $rootdir/weekly`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    $csa_dict{"$items[0]:$items[1]"}="$items[2]\t$items[3]";
}
$size=keys %csa_dict;
print "$size csa site\n";

my %uniprot_dict;
foreach my $line(`zcat $rootdir/data/chain2uniprot.tsv.gz`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $uniprot="$3";
        $uniprot_dict{"$pdbid:$asym_id"}=$uniprot;
    }
}
$size=keys %uniprot_dict;
print "$size uniprot\n";

my %go_dict;
foreach my $line(`zcat $rootdir/data/chain2go.tsv.gz`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $go     ="$3";
        $go=~s/GO://g;
        $go_dict{"$pdbid:$asym_id"}=$go;
    }
}
$size=keys %go_dict;
print "$size go\n";

my %pubmed_dict;
foreach my $line(`zcat $rootdir/data/pdb2pubmed.tsv.gz`)
{
    if ($line=~/^(\S+)\t(\d+)/)
    {
        my $pdbid="$1";
        my $pubmed="$2";
        $pubmed_dict{$pdbid}=$pubmed;
    }
}
$size=keys %pubmed_dict;
print "$size pubmed\n";

# [1] pdb
# [2] chain
# [3] resolution
# [4] csa
# [5] ec
# [6] go
# [7] uniprot
# [8] pubmed
my $txt="#pdb\tchain\tresolution\tcsa\tcsa(Renumbered)\tec\tgo\tuniprot\tpubmed\n";
foreach my $chain(@chain2ec_list)
{
    if ($chain=~/(\w+):(\w+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $reso   ="";
        my $csa    ="\t";
        my $ec     =$chain2ec_dict{$chain};
        my $go     ="";
        my $uniprot="";
        my $pubmed ="";
        $reso   =$resolu_dict{$pdbid}  if (exists $resolu_dict{$pdbid});
        $csa    =$csa_dict{$chain}     if (exists $csa_dict{$chain});
        $go     =$go_dict{$chain}      if (exists $go_dict{$chain});
        $uniprot=$uniprot_dict{$chain} if (exists $uniprot_dict{$chain});
        $pubmed =$pubmed_dict{$pdbid}  if (exists $pubmed_dict{$pdbid});
        $txt.="$pdbid\t$asym_id\t$reso\t$csa\t$ec\t$go\t$uniprot\t$pubmed\n";
    }
}
open(FP,">$rootdir/data/ec_all.tsv");
print FP "$txt";
close(FP);
&gzipFile("$rootdir/data/ec_all.tsv");

exit(0);

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
