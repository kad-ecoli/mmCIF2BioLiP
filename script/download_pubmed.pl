#!/usr/bin/perl
my $docstring=<<EOF
download_pubmed.pl pmid
    download pubmed abstract to pmid.txt
EOF
;
use strict;
use File::Basename;
use Cwd 'abs_path';

if (scalar @ARGV)
{
    my $pmid=$ARGV[0];
    my $pubmeddir=".";
    $pubmeddir=$ARGV[1] if (scalar @ARGV>1);
    &download_pubmed($ARGV[$a],$pubmeddir);
    exit(0);
}


my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);
my $pubmeddir = "$rootdir/pubmed";
system("mkdir -p $pubmeddir");

my @pdb_list;
foreach my $pdb(`grep ';' $rootdir/pdb/derived_data/index/resolu.idx|cut -f1 -d';'|grep -ohP '^\\S+'`)
{
    chomp($pdb);
    $pdb=lc($pdb);
    push(@pdb_list,($pdb));
}

my @artifact_list;
my %artifact_dict;
foreach my $line(`cat $rootdir/ligand_list`)
{
    chomp($line);
    my @items = split(/\t/, $line);
    my $resn = $items[0];
    my $desc = $items[1];
    push(@artifact_list,($resn));
    $artifact_dict{$resn}=$desc;
}

foreach my $pdb(@pdb_list)
{
    my $divided=substr($pdb,length($pdb)-3,2);
    my $outdir="$rootdir/interim/$divided";
    next if (!-s "$outdir/$pdb.txt" || !-s "$outdir/$pdb.tar.gz" || -s "$outdir/$pdb.ignore");
    my $txt=`cat $outdir/$pdb.txt`;
    my $pmid="?";
    if ($txt=~/pmid:(\d+)/) { $pmid="$1"; }
    else { next; }
    next if (-s "$pubmeddir/$pmid.txt");
    my @lines=split(/\n/, $txt);
    my $startTable=0;
    my $foundArtifact="";
    for (my $l=2;$l<scalar @lines;$l++)
    {
        my $line=$lines[$l];
        if ($line=~/^#rec/)
        {
            $startTable=1;
        }
        elsif ($startTable==1 && $line=~/^\w+\t(\w+)/)
        {
            my $resn="$1";
            if (exists $artifact_dict{$resn})
            {
                $foundArtifact="$resn";
                last;
            }
        }
    }
    next if (length $foundArtifact==0);
    print "download pmid:$pmid for $pdb to check $foundArtifact\n";
    &download_pubmed($pmid,$pubmeddir);
}


exit();

sub download_pubmed
{
    my ($pmid,$pubmeddir)=@_;

    return if (-s "$pubmeddir/$pmid.txt");
    my $abstract_txt="";
    my $title_txt="";
    my $ab=0;
    my $ti=0;
    foreach my $line(`curl 'https://pubmed.ncbi.nlm.nih.gov/$pmid/?format=pubmed'`)
    {
        $line=~s/\r//g;
        chomp($line);
        if ($line=~/^AB  - /)
        {
            $ab=1;
            $line=substr($line, 6);
            $abstract_txt=$line;
        }
        elsif ($line=~/^TI  - /)
        {
            $ti=1;
            $line=substr($line, 6);
            $title_txt=$line;
        }
        elsif ($line=~/^\w+\s+- /)
        {
            $ab=2;
            $ti=2;
        }
        elsif ($ab==1)
        {
            $line=substr($line, 6);
            $abstract_txt.="$line";
        }
        elsif ($ti==1)
        {
            $line=substr($line, 6);
            $title_txt.="$line";
        }
    }
    if (length $abstract_txt || length $title_txt)
    {
        open(FP,">$pubmeddir/$pmid.txt");
        print FP "$title_txt\n$abstract_txt\n";
        close(FP);
    }
}
