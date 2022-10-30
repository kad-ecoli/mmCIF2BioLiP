#!/usr/bin/perl
# clean up intermediate files
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

foreach my $divided(`ls $rootdir/weekly/|grep -F _nr.txt|cut -f2 -d_`)
{
    chomp($divided);
    if (-d "$rootdir/weekly/$divided")
    {
        my $cmd="rm -rf $rootdir/weekly/$divided";
        print "$cmd\n";
        system("$cmd");
    }
    
    if (-f "$rootdir/weekly/BioLiP_$divided.bsr.gz")
    {
        my $cmd="rm -rf $rootdir/weekly/BioLiP_$divided.bsr.gz";
        print "$cmd\n";
        system("$cmd");
    }
}

# remove files older than 1 hour
foreach my $filename(`find $rootdir/output/*gz -mmin +60`)
{
    chomp($filename);
    my $cmd="rm -f $filename";
    print "$cmd\n";
    system("$cmd");
}

my $mmCIF_folder="$rootdir/pdb/data/structures/divided/mmCIF";
my $mmCIF_count =`ls $mmCIF_folder/|wc -l`+0;
if ($mmCIF_count)
{
    print <<EOF
There are $mmCIF_count folders at $mmCIF_folder/, which you can choose to delete by

    rm -rf $mmCIF_folder/*

This will save ~60GB of disk space. If you choose to do so, you will have to use 
$rootdir/script/download_pdb.pl
instead of
$rootdir/script/rsyncPDB.sh
The next time you want to perform a database update.
EOF
;
}

exit();
