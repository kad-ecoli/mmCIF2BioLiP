#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

system("rm -rf $rootdir/weekly/ligand") if (-d "$rootdir/weekly/ligand");
foreach my $divided(`ls $rootdir/weekly/|grep -F ligand_|grep -F .tar.bz2|grep -vF _nr.tar.bz2|cut -f1 -d.|cut -f2 -d_`)
{
    chomp($divided);
    system("cd $rootdir/weekly; tar -xvf ligand_$divided.tar.bz2 --wildcards 'ligand/*_rna_*.pdb' 2>/dev/null");
}

my %cssr_dict;
foreach my $line(`zcat $rootdir/data/rna_ss.txt.gz|grep -v '^#'`)
{
    if ($line=~/^(\S+)\s(\S+)/)
    {
        my $target="$1";
        my $cssr="$2";
        $cssr_dict{$target}=$cssr;
    }
}


my @target_list;
foreach my $target(`zcat $rootdir/data/rna.fasta.gz |grep -ohP ">\\w+"|sed 's/>//g'`)
{
    chomp($target);
    print "$target\n";
    push(@target_list,("$target"));
    if (!exists $cssr_dict{$target})
    {
        my $cssr=`sed 's/HETATM/ATOM  /g' $rootdir/weekly/ligand/${target}_0.pdb|$bindir/CSSR - - -o 1`;
        chomp($cssr);
        $cssr_dict{$target}=$cssr;
    }
}

my $txt="#target\tcssr\n";
foreach my $target(@target_list)
{
    $txt.="$target\t$cssr_dict{$target}\n";
}
open(FP,">$rootdir/data/rna_ss.txt");
print FP "$txt";
close(FP);
&gzipFile("$rootdir/data/rna_ss.txt");

system("rm -rf $rootdir/weekly/ligand") if (-d "$rootdir/weekly/ligand");
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
