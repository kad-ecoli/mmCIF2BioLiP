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
    next if (-s "$outdir/$pdb.txt" &&  ( -s "$outdir/$pdb.ignore" ||
             -s "$outdir/$pdb.tar.gz" || -s "$outdir/$pdb.tar.bz2"));
    system("mkdir -p $outdir") if (!-d "$outdir");
    my $cmd="cd $outdir; $bindir/cif2pdb $inputdir/$pdb.cif.gz $pdb $rootdir/ligand_list";
    printf "$cmd\n";
    system("$cmd");

    next if (!-s "$outdir/$pdb.txt");
    my $receptorNum=0;
    my $ligandNum=0;
    my $startTable=0;
    foreach my $line(`cat $outdir/$pdb.txt|tail +3`)
    {
        if ($line=~/^#rec/)
        {
            $startTable=1;
        }
        elsif ($startTable==0 && $line=~/^>\S+\tprotein\t/)
        {
            $receptorNum++;
        }
        elsif ($startTable==1)
        {
            $ligandNum++;
        }
    }
    my $msg="";
    $msg.="no_protein\n" if ($receptorNum==0);
    $msg.="no_ligand\n" if ($ligandNum==0);
    next if (length $msg==0);
    print "ignore $outdir/$pdb.txt\n";
    open(FP,">$outdir/$pdb.ignore");
    print FP $msg;
    close(FP);
    system("rm $outdir/$pdb.tar.gz") if (-s "$outdir/$pdb.tar.gz");
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
foreach my $divided(`ls $rootdir/interim`)
{
    chomp($divided);
    foreach my $line(`grep --with-filename 'pmid:?' $rootdir/interim/$divided/*.txt`)
    {
        chomp($line);
        my $filename=substr($line,0,(length $line)-7);
        my $ignorefilename=substr($filename,0,(length $filename)-4).".ignore";
        my $targzfilename =substr($filename,0,(length $filename)-4).".tar.gz";
        $filename=~/(\w+)\.txt$/;
        my $pdb="$1";
        next if (!exists $pubmed_dict{$pdb} || -s "$ignorefilename" ||
                 !-s "$targzfilename");
        print "update citation for $filename\n";
        my $pubmed=$pubmed_dict{$pdb};
        system("cat $filename|gzip - > $filename.backup.gz");
        system("sed -i 's/^pmid:?/pmid:$pubmed/' $filename");
    }
}

exit();
