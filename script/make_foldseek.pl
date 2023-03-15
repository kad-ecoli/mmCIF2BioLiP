#!/usr/bin/perl
# generate database for foldseek and TMalign
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

system("mkdir -p $rootdir/foldseek/") if (! -d "$rootdir/foldseek/");
system("rm -rf $rootdir/foldseek/receptor_nr") if (-d "$rootdir/foldseek/receptor_nr");
foreach my $line(`ls $rootdir/weekly/|grep -F receptor_`)
{
    if ($line=~/receptor_(\w+)_nr.tar.bz2/)
    {
        my $divide="$1";
        my $cmd="cd $rootdir/foldseek/; tar -xf $rootdir/weekly/receptor_${divide}_nr.tar.bz2";
        print "$cmd\n";
        system("$cmd");
    }
}
system("zcat $rootdir/data/protein_nr.fasta.gz |grep -F '>'|sed 's/^>//g' > $rootdir/foldseek/list");
system("cd $rootdir/foldseek/; $bindir/pdb2xyz -dir receptor_nr/  list -suffix .pdb > receptor_nr.new");
system("cd $rootdir/foldseek/; $bindir/xyz_sfetch receptor_nr.new");
system("mv $rootdir/foldseek/receptor_nr.new $rootdir/foldseek/receptor_nr.xyz");
system("mv $rootdir/foldseek/receptor_nr.new.index $rootdir/foldseek/receptor_nr.xyz.index");
system("cd $rootdir/foldseek/; $bindir/foldseek createdb receptor_nr/ receptor_DB");
system("rm -rf $rootdir/foldseek/receptor_nr");

exit(0);
