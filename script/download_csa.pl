#!/usr/bin/perl
# download catalytic site from m-csa
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download csa\n";
system("mkdir -p $rootdir/m-csa/api/");
system("wget 'http://www.ebi.ac.uk/thornton-srv/m-csa/api/homologues_residues.json' -O $rootdir/m-csa/api/homologues_residues.json");
system("wget 'http://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json' -O $rootdir/m-csa/api/residues.json");

my $chain_name="";
my $pdb_id    ="";
my $auth_resid="";

my %chain2resid;
my @chain_list;

foreach my $line(split(/{/,`cat $rootdir/m-csa/api/residues.json $rootdir/m-csa/api/homologues_residues.json`))
{
    if ($line=~/\"chain_name\":\s*\"(\w+)\"/)
    {
        $chain_name="$1"; 
    }
    else
    {
        next;
    }
    if ($line=~/\"pdb_id\":\s*\"(\w+)\"/)
    {
        $pdb_id="$1"
    }
    else
    {
        next;
    }
    if ($line=~/\"auth_resid\":\s*(\w+)/)
    {
        $auth_resid="$1";
        next if ($auth_resid eq "null");
    }
    else
    {
        print "$line\n";
        next;
    }

    my $key="$pdb_id\t$chain_name";
    if (exists($chain2resid{$key}))
    {
        if ($chain2resid{$key}!~/\b$auth_resid\b/)
        {
            $chain2resid{$key}.=",$auth_resid";
        }
    }
    else
    {
        $chain2resid{$key}="$auth_resid";
        push(@chain_list,"$key");
    }
}

my $txt="";
foreach my $key(@chain_list)
{
    my @items=sort {$a <=> $b} split(/,/,$chain2resid{$key});
    my $auth_resid=join ',', @items;
    $txt.="$key\t$auth_resid\n";
}
if (scalar @chain_list)
{
    open(FP,">$rootdir/data/csa.tsv");
    print FP "$txt";
    close(FP);
}

exit();

