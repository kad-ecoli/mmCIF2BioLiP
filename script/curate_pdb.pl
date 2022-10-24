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
print "convert mmCIF to PDB\n";
foreach my $pdb(@pdb_list)
{
    my $divided=substr($pdb,length($pdb)-3,2);
    my $inputdir="$rootdir/pdb/data/structures/divided/mmCIF/$divided";
    my $outdir="$rootdir/interim/$divided";
    next if (-s "$outdir/$pdb.txt" && 
            (-s "$outdir/$pdb.tar.gz" || -s "$outdir/$pdb.tar.bz2"));
    system("mkdir -p $outdir") if (!-d "$outdir");
    my $cmd="cd $outdir; $bindir/cif2pdb $inputdir/$pdb.cif.gz $pdb $rootdir/ligand_list";
    printf "$cmd\n";
    system("$cmd");
}

print "update pubmed citation\n";
my %pubmed_dict;
foreach my $line(`zcat $rootdir/data/pdb2pubmed.tsv.gz`)
{
    if ($line=~/^(\w+)\t(\d+)/)
    {
        my $pdb="$1";
        my $pubmed="$2";
        $pubmed_dict{$pdb}=$pubmed;
    }
}
foreach my $pdb(@pdb_list)
{
    my $divided=substr($pdb,length($pdb)-3,2);
    my $outdir="$rootdir/interim/$divided";
    next if (!exists $pubmed_dict{$pdb} || !-s "$outdir/$pdb.txt");
    my $pubmed=`head -2 $outdir/$pdb.txt|tail -1`;
    chomp($pubmed);
    next if ($pubmed=~/^\d+$/);
    print "updating pubmed for $outdir/$pdb.txt\n";
    $pubmed=$pubmed_dict{$pdb};
    system("cat $outdir/$pdb.txt|gzip - > $outdir/$pdb.txt.backup.gz");
    my $txt=`head -1 $outdir/$pdb.txt`;
    $txt.="$pubmed\n";
    $txt.=`tail +3 $outdir/$pdb.txt`;
    open(FP,">$outdir/$pdb.txt");
    print FP $txt;
    close(FP);
}

exit();
