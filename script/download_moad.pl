#!/usr/bin/perl
# download binding affinity from MOAD
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download MOAD\n";
system("mkdir -p $rootdir/bind");
system("wget -q 'http://www.bindingmoad.org/files/csv/every.csv' -O $rootdir/bind/every.csv");

if (!-s "$rootdir/bind/every.csv")
{
    print "ERROR! cannot download moad\n";
    exit(1);
}

my @moad_list;
my %moad_dict;
my $pdbid;
foreach my $line(`cut -f3-9 -d, $rootdir/bind/every.csv`)
{
    my @items=split(/,/,$line);
    if (length $items[0])
    {
        $pdbid=lc($items[0]);
    }
    elsif (length $items[4])
    {
        my $ligand=$items[1];
        if ($ligand=~/^(\w+):(\w+)/)
        {
            my $ccd    ="$1";
            my $asym_id="$2";
            my $kd="$items[3]$items[4]$items[5]$items[6]";
            my $key="$pdbid\t$1\t$2";
            next if (exists $moad_dict{$key});
            $moad_dict{$key}=$kd;
            push(@moad_list,($key));
        }
    }
}
my $txt="";
foreach my $key(@moad_list)
{
    $txt.="$key\t$moad_dict{$key}";
}
open(FP,">$rootdir/data/moad.tsv");
print FP "$txt";
close(FP);

exit();
