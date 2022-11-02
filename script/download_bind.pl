#!/usr/bin/perl
# download binding affinity from BindingDB and PDBbind-CN
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download BindingDB\n";
system("mkdir -p $rootdir/bind");
system("wget 'https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp' -O $rootdir/bind/Download.jsp");
if (`cat $rootdir/bind/Download.jsp`=~/(BindingDB_All_\w+.tsv.zip)/)
{
    my $outfile="$1";
    system("wget 'https://www.bindingdb.org/rwd/bind/downloads/$outfile' -O $rootdir/bind/BindingDB_All_tsv.zip");
    if (-s "$rootdir/bind/BindingDB_All_tsv.zip")
    {
        system("zcat $rootdir/bind/BindingDB_All_tsv.zip| cut -f9-14,27,28,42,43,47|  head -1 > $rootdir/bind/BindingDB.tsv");
        system("zcat $rootdir/bind/BindingDB_All_tsv.zip| cut -f9-14,27,28,42,43,47|grep -P '\\S+\\t\\S+\\t\\S*\\t\\S*\\t\\S*\$' |grep -vP '\\t\\t\\t\$'  >> $rootdir/bind/BindingDB.tsv");
        &gzipFile("$rootdir/bind/BindingDB.tsv");
    }
}

my %pdb2chain;
my %chain2uniprot;
foreach my $line(`zcat $rootdir/data/chain2uniprot.tsv.gz`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)$/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $uniprot="$3";
        if (!exists $pdb2chain{$pdbid}) { $pdb2chain{$pdbid}="$asym_id"; }
        else                          { $pdb2chain{$pdbid}.=",$asym_id"; }
        $chain2uniprot{"$pdbid\t$asym_id"}=$uniprot;
    }
}

#  0 Ki (nM)
#  1 IC50 (nM)
#  2 Kd (nM)
#  3 EC50 (nM)
#  4 kon (M-1-s-1)
#  5 koff (s-1)
#  6 Ligand HET ID in PDB
#  7 PDB ID(s) for Ligand-Target Complex
#  8 UniProt (SwissProt) Primary ID of Target Chain
#  9 UniProt (SwissProt) Secondary ID(s) of Target Chain
# 10 UniProt (TrEMBL) Primary ID of Target Chain
my @header_list=split(/\t/,`zcat $rootdir/bind/BindingDB.tsv.gz|head -1`);
my @BindingDB_list;
my %BindingDB_dict;
foreach my $line(`zcat $rootdir/bind/BindingDB.tsv.gz|tail -n +2`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $uniprot_line="$items[8],$items[9],$items[10]";

    my $affinity="";
    for (my $i=0;$i<6;$i++)
    {
        if ($items[$i]=~/(\S+)/)
        {
            my $kd="$1";
            my $first=$header_list[$i];
            my $second="";
            if ($first=~/(\S+)\s\((\S+)\)/)
            {
                $first="$1";
                $second="$2";
            }
            if (length $affinity) { $affinity.=", $first=$kd$second"; }
            else                  { $affinity =  "$first=$kd$second"; }
        }
    }
    next if (length $affinity==0);

    my $ccd     =   $items[6];
    my $pdb_line=lc($items[7]);
    foreach my $pdbid(split(/,/,$pdb_line))
    {
        foreach my $asym_id(split(/,/,$pdb2chain{$pdbid}))
        {
            my $uniprot=$chain2uniprot{"$pdbid\t$asym_id"};
            if ($uniprot_line=~/\b$uniprot\b/)
            {
                my $key="$pdbid\t$asym_id\t$ccd";
                if (!exists($BindingDB_dict{$key}))
                {
                    push(@BindingDB_list,($key));
                    $BindingDB_dict{$key}=$affinity;
                }
                else
                {
                    $BindingDB_dict{$key}.=", $affinity";
                }
            }
        }
    }
}
my $txt="";
foreach my $key(@BindingDB_list)
{
    my $value=$BindingDB_dict{$key};
    my %exist_key;
    my $uniq_value="";
    foreach my $kd(split(/, /,$value))
    {
        if ($kd=~/^(\S+)=/)
        {
            my $key="$1";
            next if (exists $exist_key{$key});
            $exist_key{$key}=1;
            if (length $uniq_value==0) { $uniq_value=$kd; }
            else { $uniq_value.=", $kd"; }
        }
    }
    $txt.="$key\t$uniq_value\n";
}
open(FP,">$rootdir/data/BindingDB.tsv");
print FP "$txt";
close(FP);
&gzipFile("$rootdir/data/BindingDB.tsv");

exit();

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
