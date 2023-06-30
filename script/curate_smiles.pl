#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

my $infile ="$rootdir/pdb/data/monomers/components.cif.gz";
if (!-s "$infile")
{
    print "ERROR! cannot find $infile\n";
    exit(1);
}
my $outfile="$rootdir/data/smiles.tsv";
my $txt="";
foreach my $line(`zcat $infile | grep -P "( SMILES )|( SMILES_CANONICAL )"| sed 's/ SMILES_CANONICAL / /g'|sed 's/ SMILES / /g'`)
{
    if ($line=~/^(\S+)\s+(\w+)\s+(\S+)\s+(\S+)/ ||
        $line=~/^(\S+)\s+\"([\w\s]+)\"\s+(\S+)\s+(\S+)/)
    {
        my $lig3   ="$1";
        my $method ="$2";
        my $version="$3";
        my $smiles ="$4";
        if ($smiles=~/\"(\S+)\"/)
        {
            $smiles="$1";
        }
        $txt.="$lig3\t$smiles\t$method $version\n";
    }
}
open(FP,">$outfile.tmp");
print FP "$txt";
close(FP);
system("cat $outfile.tmp|uniq|sort|uniq > $outfile");
system("rm  $outfile.tmp");

&gzipFile("$outfile");
exit(0);

sub strip
{
    my ($instring)=@_;
    $instring=~s/^\s*//;
    $instring=~s/\s*$//;
    $instring=~s/^"//;
    $instring=~s/"$//;
    return $instring;
}

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
