#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

my $txt="";
my $size;
my %pdb2receptor;
my %pdb2ligand;

foreach my $line(`cd $rootdir/bind/pocket/; grep -F ATOM *_pocket.pdb|grep -F ' CA '|sed 's/_pocket.pdb:ATOM//g'|awk '{print \$1,\$5}'|grep -ohP "^\\w+\\s\\w"|uniq`)
{
    if ($line=~/(\w+)\s(\w+)/)
    {
        my $pdbid="$1";
        my $asym_id="$2";
        if (!exists $pdb2receptor{$pdbid}){$pdb2receptor{$pdbid}=$asym_id;}
        else      { $pdb2receptor{$pdbid}.=",$asym_id";}
    }
}
$size=keys %pdb2receptor;
print "$size receptor\n";



foreach my $line(`cd $rootdir/bind/ligand/; grep -v '#' *_ligand.mol2|sed 's/_ligand.mol2:/ /g'|awk '{ print \$1, \$7, \$8, \$9; }'|grep -P '\\w+\\s[A-Z][.\\w]*\\s\\d+\\s\\w+'|cut -f1,4 -d' '|sort|uniq`)
{
    if ($line=~/(\w+)\s(\w+)/)
    {
        my $pdbid="$1";
        my $comp_id="$2";
        if (!exists $pdb2ligand{$pdbid}){$pdb2ligand{$pdbid}=$comp_id;}
        else      { $pdb2ligand{$pdbid}.=",$comp_id";}
    }
}
$size=keys %pdb2ligand;
print "$size ligand\n";



my %pdb2kd;
my @pdbid_list;
foreach my $filename(`ls $rootdir/bind/index/|grep '^INDEX'`)
{
    chomp($filename);
    my @fields;
    foreach my $line(`cat $rootdir/bind/index/$filename`)
    {
        chomp($line);
        if ($line=~/^# PDB code,/)
        {
            my @items=split(/,/,$line);
            for (my $i=0;$i<scalar @items;$i++)
            {
                if ($items[$i]=~/binding data/ || $items[$i]=~/Kd\/Ki/)
                {
                    push(@fields,("$i"));
                }
            }
        }
        elsif ($line!~/^#/)
        {
            my @items=split(/\.pdf /,$line);
            my $first=$items[0];
            my $second=$items[1];
            @items=split(/\s+/,$first);
            my $pdbid=$items[0];
            my $kd="";
            foreach (my $f=0;$f<scalar @fields;$f++)
            {
                my $value=$items[$fields[$f]];
                $value="-logKd/Ki=$value" if ($value=~/^[.\d]+/);
                if ($f==0) { $kd=$value; }
                else    { $kd.=",$value"; }
            }
            if (!exists($pdb2kd{$pdbid}))
            {
                $pdb2kd{$pdbid}="$kd\t$second";
                push(@pdbid_list,($pdbid));
            }
            else
            {
                print "WARNING! Duplicated affinity $line\n";
            }
        }
    }
}
$size=keys %pdb2kd;
print "$size binding affinity\n";



$txt="";
foreach my $pdbid(@pdbid_list)
{
    my $asym_id_list="_";
    my $comp_id_list="_";
    if (exists($pdb2receptor{$pdbid}))
    {
        $asym_id_list="$pdb2receptor{$pdbid}";
    }
    if (exists($pdb2ligand{$pdbid}))
    {
        $comp_id_list="$pdb2ligand{$pdbid}";
    }
    else
    {
        if ($pdb2kd{$pdbid}=~/RNA\b/) {$comp_id_list="rna";}
        if ($pdb2kd{$pdbid}=~/DNA\b/) 
        {
            if ($comp_id_list eq "_") { $comp_id_list="dna"; }
            else                    { $comp_id_list.=",dna"; }
        }
        if ($comp_id_list eq "_")
        {
            if    ($pdb2kd{$pdbid}=~/RNA/) {$comp_id_list="rna";}
            elsif ($pdb2kd{$pdbid}=~/DNA/) {$comp_id_list="dna";}
            else                       {$comp_id_list="rna,dna";};
        }
    }
    foreach my $asym_id(split(/,/,$asym_id_list))
    {
        foreach my $comp_id(split(/,/,$comp_id_list))
        {
            $txt.="$pdbid\t$asym_id\t$comp_id\t$pdb2kd{$pdbid}\n";
        }
    }
}
open(FP,">$rootdir/data/PDBbind.tsv");
print FP "$txt";
close(FP);
exit();
