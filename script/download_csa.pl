#!/usr/bin/perl
# download catalytic site from m-csa
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download csa\n";
system("mkdir -p $rootdir/m-csa/api/");
system("wget 'https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json' -O $rootdir/m-csa/api/residues.json");

my $chain_name="";
my $pdb_id    ="";
my $auth_resid="";

my %chain2resid;
my @chain_list;

foreach my $line(`grep -ohP '(\\"chain_name\\":\\s*\\"\\w+\\")|(\\"pdb_id\\":\\s*\\"\\w+\\")|(\\"auth_resid\\":\\s*\\d+)' $rootdir/m-csa/api/residues.json`)
{
    $line=~/(\S+):(\S+)/;
    my $key="$1";
    my $value="$2";
    $key=~s/^"|"$//g;
    $value=~s/^"|"$//g;
    
    if ($key eq "chain_name") { $chain_name=$value; }
    elsif ($key eq "pdb_id")  { $pdb_id    =$value; }
    elsif ($key eq "auth_resid")
    {
        $auth_resid=$value;
        $key       ="$pdb_id\t$chain_name";
        if (exists($chain2resid{$key}))
        {
            $chain2resid{$key}.=",$auth_resid";
        }
        else
        {
            $chain2resid{$key}="$auth_resid";
            push(@chain_list,"$key");
        }
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

