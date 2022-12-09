#!/usr/bin/perl
my $docstring=<<EOF
fsearch.pl jobID
    perform foldseek and TMalign on output/jobID.cif or output/jobID.pdb
EOF
;

use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

if (scalar @ARGV!=1)
{
    print "$docstring";
    exit(1);
}
my $input =$ARGV[0];
my $output="$rootdir/output/$input";

if (!-d "$output/")
{
    system("mkdir -p $output/");
}

my $infile="input.pdb";
if (-s "$output.pdb")
{
    system("cp $output.pdb $output/input.pdb");
}
else
{
    if (!-s "$output.cif")
    {
        print "no such file $output.cif\n";
        exit(1);
    }
    system("cp $output.cif $output/input.cif");
    $infile="input.cif";
}

if (!-s "$output/aln.m8")
{
    system("rm -rf $output/tmpFolder/") if (-d "$output/tmpFolder");
    system("cd $output; $bindir/foldseek easy-search $infile $rootdir/foldseek/receptor_DB aln.m8 tmpFolder --threads 1 --format-output target,alnlen,nident,qlen,tlen,evalue -e 1");
    if (!-s "$output/aln.m8")
    {
        print "no foldseek hit\n";
        exit(2);
    }
    system("rm -rf $output/tmpFolder");
}
my @m8_list;
my %m8_dict;
my $txt;
foreach my $line(`cat $output/aln.m8`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $target=$items[0];
    $target=~s/.pdb$//;
    next if (exists $m8_dict{$target});
    my $evalue=$items[5];
    last if (scalar @m8_list>=10 && $evalue>0.001);
    push(@m8_list,("$target"));
    $m8_dict{$target}=$evalue;
    $txt.="$target\n";
}

if (!-s "$output/TMalign.tsv" && length $txt)
{
    #system("cut -f1 -d. $output/aln.m8 > $output/list");
    open(FP,">$output/list");
    print FP "$txt";
    close(FP);
    system("$bindir/xyz_sfetch $rootdir/foldseek/receptor_nr.xyz $output/list | gzip - > $output/xyz.gz");
    system("cd $output/; $bindir/USalign $infile xyz.gz -infmt2 2 -fast -outfmt 2 |grep -v '^#'|cut -f2-|sed 's/^xyz.gz://g' > TMalign.tsv");
}

#### prepare output ####
$txt=`cat $rootdir/index.html`;
my @items=split(/<!-- CONTENT START -->/,$txt);
my $html_header=$items[0];
@items   =split(/<!-- CONTENT END -->/,$txt);
my $html_end=$items[1];
@items   =split(/<!-- ==BBB== ======ending of head======== -->/,$html_header);
$html_header=~s/href="/href="..\/..\//g;
$html_header=~s/href="..\/..\/http/href="http/g;
$html_header=~s/href="..\/..\/\//href="\//g;
open(FP,">$output/index.html");
print FP <<EOF
$html_header
Search <a href=$infile>$infile</a> through BioLiP using Foldseek and US-align (-fast mode). Results are ranked in descending order of E-value.<br>

<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=5% ALIGN=center><strong> # </strong></th>
    <th width=10% ALIGN=center><strong> Hit<br>(length)</strong></th>
    <th width=10% ALIGN=center><strong> Aligned<br>length </strong></th>
    <th width=10% ALIGN=center><strong> TM-score (Identity)<br>normalized by query</strong> </th>           
    <th width=10% ALIGN=center><strong> TM-score (Identity)<br>normalized by hit</strong> </th>           
    <th width=10% ALIGN=center><strong> RMSD (Identity)<br>normalized by<br>aligned length</strong> </th>           
    <th width=10% ALIGN=center><strong> Foldseek<br>E-value </strong> </th>
    <th width=35% ALIGN=center><strong> Homologs<br>to hit</strong> </th>           
</tr><tr ALIGN=center>
EOF
;


my %hit2chain_dict;
my %hit2clust_dict;
if (scalar @m8_list)
{
    foreach my $line(`zcat $rootdir/data/protein.fasta.gz|grep '>'|cut -f1,2`)
    {
        chomp($line);
        if ($line=~/>(\w+)\t(\w+)/)
        {
            my $hit="$1";
            my $asym_id="$2";
            my $pdbid=substr($hit,0,(length $hit)-(length $asym_id));
            $hit2chain_dict{$hit}="$pdbid:$asym_id";
        }
    }
    foreach my $line(`zcat $rootdir/data/protein_nr.fasta.clust.gz`)
    {
        if ($line=~/(\S+)\t(\S+)/)
        {
            my $rep="$1";
            my $mem="$2";
            $hit2clust_dict{$rep}=$mem;
        }
    }
}

my %line_dict;
foreach my $line(`cat $output/TMalign.tsv`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    $line_dict{$items[0]}=$line;
}

my $totalNum=0;
foreach my $target(@m8_list)
{
    my @items=split(/\t/,$line_dict{$target});
    my $sacc =$items[0];
    my $TM1  =$items[1];
    my $TM2  =$items[2];
    my $RMSD =$items[3];
    my $ID1  =$items[4];
    my $ID2  =$items[5];
    my $IDali=$items[6];
    my $L1   =$items[7];
    my $L2   =$items[8];
    my $Lali =$items[9];
    my $evalue="NA";
       $evalue=$m8_dict{$sacc} if (exists $m8_dict{$sacc});
    $totalNum++;
    my $bgcolor='';
       $bgcolor='BGCOLOR="#DEDEDE"' if ($totalNum % 2==0);

    my $pdbid="$sacc";
    my $asym_id="";
    if ($hit2chain_dict{$sacc}=~/(\w+):(\w+)/)
    {
        $pdbid="$1";
        $asym_id="$2";
    }
    my $hit="<a href=../../qsearch.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$pdbid:$asym_id</a>";
    my $homolog_line;
    if (exists $hit2clust_dict{$sacc})
    {
        foreach my $mem (split(/,/,$hit2clust_dict{$sacc}))
        {
            if ($hit2chain_dict{$mem}=~/(\w+):(\w+)/)
            {
                my $pdbid="$1";
                my $asym_id="$2";
                $homolog_line.=", <a href=../../qsearch.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$pdbid:$asym_id</a>";
            }
        }
    }
    $homolog_line=substr($homolog_line,2) if (length $homolog_line);
    print FP <<EOF
<tr $bgcolor ALIGN=center>
    <td>$totalNum</td>
    <td>$hit ($L2)</td>
    <td>$Lali</td>
    <td>$TM1 ($ID1)</td>
    <td>$TM2 ($ID2)</td>
    <td>$RMSD ($IDali)</td>
    <td>$evalue</td>
    <td>$homolog_line</td>
</tr>
EOF
;
}
print FP <<EOF
</table><p></p><a href=../..>[Back]</a>
$html_end
EOF
;
close(FP);
exit(0);
