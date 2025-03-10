#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download PDB\n";
system("mkdir -p $rootdir/pdb/derived_data/index/");
system("wget -q http://files.wwpdb.org/pub/pdb/derived_data/index/resolu.idx -O $rootdir/pdb/derived_data/index/resolu.idx");
system("wget -q ftp://files.wwpdb.org/pub/pdb/derived_data/index/resolu.idx -O $rootdir/pdb/derived_data/index/resolu.idx") if (!-s "$rootdir/pdb/derived_data/index/resolu.idx");
if (!-s "$rootdir/pdb/derived_data/index/resolu.idx")
{
    print "ERROR! cannot download $rootdir/pdb/derived_data/index/resolu.idx\n";
    exit(1);
}
#system("wget -q https://ftp.wwpdb.org/pub/pdb/derived_data/index/compound.idx -O $rootdir/pdb/derived_data/index/compound.idx");

foreach my $pdb(`grep ';' $rootdir/pdb/derived_data/index/resolu.idx|cut -f1 -d';'|grep -ohP '^\\S+'`)
{
    chomp($pdb);
    $pdb=lc($pdb);
    my $divided=substr($pdb,length($pdb)-3,2);
    my $inputdir="$rootdir/pdb/data/structures/divided/mmCIF/$divided";
    system("mkdir -p $inputdir") if (!-d $inputdir);
    next if (-s "$inputdir/$pdb.cif.gz");
    my $outdir="$rootdir/interim/$divided";
    next if (-s "$outdir/$pdb.txt" &&  ( -s "$outdir/$pdb.ignore" ||
             -s "$outdir/$pdb.tar.gz" || -s "$outdir/$pdb.tar.bz2"));
    my $cmd="wget -q https://files.wwpdb.org/pub/pdb/data/structures/all/mmCIF/$pdb.cif.gz -O $inputdir/$pdb.cif.gz";
    print "$cmd\n";
    system("$cmd");
    if (!-s "$inputdir/$pdb.cif.gz")
    {
        my $cmd="wget -q https://files.rcsb.org/download/$pdb.cif.gz -O $inputdir/$pdb.cif.gz";
        print "$cmd\n";
        system("$cmd");
        if (!-s "$inputdir/$pdb.cif.gz")
        {
            print "ERROR! cannot download $inputdir/$pdb.cif.gz\n";
        }
    }
}

exit();
