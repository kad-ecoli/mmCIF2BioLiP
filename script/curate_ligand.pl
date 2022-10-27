#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

my @pdb_list;
foreach my $pdb(`grep ';' $rootdir/pdb/derived_data/index/resolu.idx|cut -f1 -d';'|grep -ohP '^\\S+'`)
{
    chomp($pdb);
    $pdb=lc($pdb);
    push(@pdb_list,($pdb));
}
print "remove irrelevant ligand\n";
foreach my $pdb(@pdb_list)
{
    my $divided=substr($pdb,length($pdb)-3,2);
    my $outdir="$rootdir/interim/$divided";
    next if (!-s "$outdir/$pdb.txt" ||  -s "$outdir/$pdb.ignore" ||
             !-s "$outdir/$pdb.tar.gz");
    next if (-s "$outdir/$pdb.tar.bz2" && 
        `grep -F 'pmid:?' $outdir/$pdb.txt`=~/pmid:/);
    my $cmd="cd $outdir; $bindir/rmligand $pdb $rootdir/ligand_list $rootdir/pubmed";
    printf "$cmd\n";
    system("$cmd");
}
exit();
