#!/usr/bin/perl
# make non redundant dataset
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

my %fasta_dict;
my @lines=`zcat $rootdir/data/protein.fasta.gz`;
for (my $l=0;$l<scalar @lines;$l++)
{
    my $header=$lines[$l];
    my $sequence=$lines[$l+1];
    chomp($sequence);
    $header=~/>(\S+)\t/;
    my $chain ="$1";
    $fasta_dict{$chain}=$sequence;
}
my $size=keys %fasta_dict;
print "$size sequence\n";

my %nr_dict;
foreach my $line(`zcat $rootdir/data/protein_nr.fasta.gz |grep '>'|cut -f1|sed 's/>//g'`)
{
    chomp($line);
    $nr_dict{$line}=1;
}
$size=keys %nr_dict;
print "$size non-redundant chain\n";

my %resolu_dict;
foreach my $line(`cat $rootdir/pdb/derived_data/index/resolu.idx`)
{
    if ($line=~/^(\S+)\t;\t(\S+)/)
    {
        my $pdbid=lc($1);
        my $resolu="$2";
        $resolu_dict{$pdbid}=$resolu;
    }
}
$size=keys %resolu_dict;
print "$size resolution\n";

my %ec_dict;
foreach my $line(`zcat $rootdir/data/chain2ec.tsv.gz`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $ec     ="$3";
        $ec_dict{$pdbid.$asym_id}=$ec;
    }
}
$size=keys %ec_dict;
print "$size ec\n";

my %go_dict;
foreach my $line(`zcat $rootdir/data/chain2go.tsv.gz`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $go     ="$3";
        $go=~s/GO://g;
        $go_dict{$pdbid.$asym_id}=$go;
    }
}
$size=keys %go_dict;
print "$size go\n";

my %uniprot_dict;
foreach my $line(`zcat $rootdir/data/chain2uniprot.tsv.gz`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $uniprot="$3";
        $uniprot_dict{$pdbid.$asym_id}=$uniprot;
    }
}
$size=keys %uniprot_dict;
print "$size uniprot\n";

my %pubmed_dict;
foreach my $line(`zcat $rootdir/data/pdb2pubmed.tsv.gz`)
{
    if ($line=~/^(\S+)\t(\d+)/)
    {
        my $pdbid="$1";
        my $pubmed="$2";
        $pubmed_dict{$pdbid}=$pubmed;
    }
}
$size=keys %pubmed_dict;
print "$size pubmed\n";

my %affman_dict;
foreach my $line(`cat $rootdir/data/affman.tsv`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+\s\S+)/)
    {
        my $pdbid="$1";
        my $asym_id="$2";
        my $affman="$3";
        $affman_dict{$pdbid.$asym_id}=$affman;
    }
}
$size=keys %affman_dict;
print "$size affinity from manual survey of literature\n";

my %moad_dict;
foreach my $line(`cat $rootdir/data/moad.tsv`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $comp_id="$2";
        my $asym_id="$3";
        my $kd     ="$4";
        $moad_dict{$pdbid.$asym_id."\t".$comp_id}=$kd;
    }
}
$size=keys %moad_dict;
print "$size affinity from moad\n";

my %bindcn_dict;
foreach my $line(`cat $rootdir/data/PDBbind.tsv`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $comp_id="$3";
        my $kd     ="$4";
        $bindcn_dict{$pdbid.$asym_id."\t".$comp_id}=$kd;
    }
}
$size=keys %bindcn_dict;
print "$size affinity from PDBbind-CN\n";

my %binddb_dict;
foreach my $line(`zcat $rootdir/data/BindingDB.tsv.gz`)
{
    chomp($line);
    $line=~s/, /,/g;
    if ($line=~/^(\S+)\t(\S+)\t(\S+)\t([\s\S]+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $comp_id="$3";
        my $kd     ="$4";
        $binddb_dict{$pdbid.$asym_id."\t".$comp_id}=$kd;
    }
}
$size=keys %binddb_dict;
print "$size affinity from BindingDB\n";

my %csa_dict;
foreach my $line(`cat $rootdir/data/csa.tsv`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my $csaOrig="$3";
        my $csaRenu="";

        my $divided=substr($pdbid,(length $pdbid)-3,2);
        my $filename="$rootdir/weekly/$divided/receptor/$pdbid$asym_id.pdb";
           $filename="$rootdir/weekly/$divided/receptor_nr/$pdbid$asym_id.pdb" if (!-s "$filename");
        next if (!-s "$filename");
        my @res_list=split(/,/,$csaOrig);
        my %res_dict=map { $_, 0 } @res_list; 
        
        my $r=0;
        foreach my $resSeq(`grep -F ' CA ' $filename|cut -c23-27`)
        {
            $r++;
            my $resi=substr($resSeq,0,4);
            $resi=~s/ //g;
            next if (!exists $res_dict{$resi});
            my $icode=substr($resSeq,4);
            $res_dict{$resi}=$r if ($icode eq " " || $res_dict{$resi}==0);
        }
        foreach my $resi(@res_list)
        {
            $csaRenu.=",$res_dict{$resi}";
        }
        $csaRenu=~s/^,//;
        $csa_dict{$pdbid.$asym_id}="$csaOrig\t$csaRenu";
    }
}
$size=keys %csa_dict;
print "$size csa site\n";


my $lig_full_all="";
my $lig_nr_all  ="";
my $pdb_full_all="";
my $pdb_nr_all  ="";
foreach my $divided(`ls $rootdir/weekly/|grep -P "BioLiP_\\w+\\.bsr\\.gz"|cut -f1 -d.|cut -f2 -d_`)
{
    chomp($divided);
    my $infile="$rootdir/weekly/BioLiP_$divided.bsr.gz";
    print "converting $infile\n";
    my $txt_full;
    my $txt_nr;
    my @receptor_nr;
    my @ligand_nr;
    foreach my $line(`zcat $infile`)
    {
        chomp($line);
        my @items  =split(/\t/,"$line");
        my $pdbid  =$items[0]; # [1] pdb id
        my $recCha =$items[1]; # [2] receptor chain
        my $bs     =$items[2]; # [4] binding site
        my $ccd    =$items[3]; # [5] ligand residue name
        my $ligCha =$items[4]; # [6] ligand chain
        my $ligIdx =$items[5]; # [7] ligand serial number
        my $resOrig=$items[6]; # [8] binding site residue, original numbering
        my $resRenu=$items[7]; # [9] binding site residue, renumbered from 1
        
        my $chain  =$pdbid.$recCha;

        my $resolu =""; # [3]  resolution
           $resolu =$resolu_dict{$pdbid} if (exists $resolu_dict{$pdbid});
        my $csaOrig=""; # [10] catalytic site atlas original numbering
        my $csaRenu=""; # [11] catalytic site atlas renumbered from 1
        my $ec     =""; # [12] EC number
           $ec     =$ec_dict{$chain} if (exists $ec_dict{$chain});
        my $go     =""; # [13] GO term
           $go     =$go_dict{$chain} if (exists $go_dict{$chain});
        my $affman =""; # [14] binding affinity from manual survey of pubmed
           $affman =$affman_dict{$chain} if (exists $affman_dict{$chain});
        my $moad   =""; # [15] binding affinity from Binding MOAD
        if (exists  $moad_dict{$chain."\t".$ccd})
        {
            $moad  =$moad_dict{$chain."\t".$ccd};
        }
        my $bindcn =""; # [16] binding affinity from PDBbind-CN
        if (exists  $bindcn_dict{$chain."\t".$ccd})
        {
            $bindcn=$bindcn_dict{$chain."\t".$ccd};
        }
        elsif (exists $bindcn_dict{$pdbid."_\t".$ccd})
        {
            $bindcn=$bindcn_dict{$pdbid."_\t".$ccd};
        }
        my $binddb =""; # [17] binding affinity from BindingDB
        if (exists  $binddb_dict{$chain."\t".$ccd})
        {
            $binddb=$binddb_dict{$chain."\t".$ccd};
        }
        my $uniprot=""; # [18] uniprot accession
           $uniprot=$uniprot_dict{$chain} if (exists $uniprot_dict{$chain});
        my $pubmed =""; # [19] pubmed id
           $pubmed =$pubmed_dict{$pdbid} if (exists $pubmed_dict{$pdbid});
        my $sequence=$fasta_dict{$chain};#[20] receptor sequence
        if (exists $csa_dict{$chain})
        {
            my $line=$csa_dict{$chain};
            if ($line=~/(\S+)\t(\S+)/)
            {
                $csaOrig="$1";
                $csaRenu="$2";
                my @csaOrig_list=split(/,/,$csaOrig);
                my @csaRenu_list=split(/,/,$csaRenu);
                $csaOrig="";
                $csaRenu="";
                for (my $i=0;$i<scalar @csaRenu_list;$i++)
                {
                    my $rOrig=$csaOrig_list[$i];
                    my $rRenu=$csaRenu_list[$i];
                    my $aa=substr($sequence,$rRenu-1,1);
                    $csaOrig.=" $aa$rOrig";
                    $csaRenu.=" $aa$rRenu";
                }
                $csaOrig=~s/^\s//g;
                $csaRenu=~s/^\s//g;
            }
        }
        $line="$pdbid\t$recCha\t$resolu\t$bs\t$ccd\t$ligCha\t$ligIdx\t".
              "$resOrig\t$resRenu\t$csaOrig\t$csaRenu\t$ec\t$go\t$affman\t".
              "$moad\t$bindcn\t$binddb\t$uniprot\t$pubmed\t$sequence\n";
        $txt_full.="$line";
        next if (!exists $nr_dict{$chain});
        $txt_nr.="$line";
        push(@receptor_nr,"$chain.pdb");
        push(@ligand_nr,"${pdbid}_${ccd}_${ligCha}_$ligIdx.pdb");
    }
    open(FP,">$rootdir/weekly/BioLiP_$divided.txt");
    print FP $txt_full;
    close(FP);
    open(FP,">$rootdir/weekly/BioLiP_${divided}_nr.txt");
    print FP $txt_nr;
    close(FP);

    $lig_full_all.=`cut -f1,2,4-9,14-17 $rootdir/weekly/BioLiP_$divided.txt`;
    $lig_nr_all.=`cut -f1,2,4-9,14-17 $rootdir/weekly/BioLiP_${divided}_nr.txt`;
    $pdb_full_all.=`cut -f1-3,10-13,18-19 $rootdir/weekly/BioLiP_$divided.txt|uniq`;
    $pdb_nr_all.=`cut -f1-3,10-13,18-19 $rootdir/weekly/BioLiP_${divided}_nr.txt|uniq`;

    system("mkdir -p $rootdir/weekly/$divided/receptor_nr");
    system("mkdir -p $rootdir/weekly/$divided/receptor_nr1");
    system("mkdir -p $rootdir/weekly/$divided/ligand_nr");
    my %filename_dict=map {$_ => 1} @receptor_nr;
    my @filename_list=keys %filename_dict;
    foreach my $filename(@filename_list)
    {
        next if (-s "$rootdir/weekly/$divided/receptor_nr/$filename"  &&
              #  -s "$rootdir/weekly/$divided/receptor_nr1/$filename" &&
              # !-s "$rootdir/weekly/$divided/receptor1/$filename"    &&
                !-s "$rootdir/weekly/$divided/receptor/$filename");
        system("mv $rootdir/weekly/$divided/receptor/$filename $rootdir/weekly/$divided/receptor_nr/$filename");
        #system("mv $rootdir/weekly/$divided/receptor1/$filename $rootdir/weekly/$divided/receptor_nr1/$filename");
    }
    %filename_dict=map {$_ => 1} @ligand_nr;
    @filename_list=keys %filename_dict;
    foreach my $filename(@filename_list)
    {
        next if (-s "$rootdir/weekly/$divided/ligand_nr/$filename" &&
                !-s "$rootdir/weekly/$divided/ligand/$filename");
        system("mv $rootdir/weekly/$divided/ligand/$filename $rootdir/weekly/$divided/ligand_nr/$filename");
    }
    #system("cd $rootdir/weekly/$divided/; tar -cjf $rootdir/weekly/receptor1_${divided}_nr.tar.bz2 receptor_nr1/");
    system("cd $rootdir/weekly/$divided/; tar -cjf $rootdir/weekly/receptor_${divided}_nr.tar.bz2 receptor_nr/");
    system("cd $rootdir/weekly/$divided/; tar -cjf $rootdir/weekly/ligand_${divided}_nr.tar.bz2 ligand_nr/");
}

print "generating $rootdir/data/{pdb,lig}_{all,nr}.tsv\n";
open(FP,">$rootdir/data/pdb_all.tsv");
print FP $pdb_full_all;
close(FP);
open(FP,">$rootdir/data/pdb_nr.tsv");
print FP $pdb_nr_all;
close(FP);
open(FP,">$rootdir/data/lig_all.tsv");
print FP $lig_full_all;
close(FP);
open(FP,">$rootdir/data/lig_nr.tsv");
print FP $lig_nr_all;
close(FP);
foreach my $prefix(("pdb_all","pdb_nr","lig_all","lig_nr"))
{
    #system("gzip -f $rootdir/data/$prefix.tsv");
    my $txt=`sort $rootdir/data/$prefix.tsv|uniq`;
    open(FP,">$rootdir/data/$prefix.tsv");
    if ($prefix=~/lig/)
    {
        print FP "#pdb\tchain\tBS\tCCD\tligandChain\tligIdx\tresidue\tresidueRenumbered\taffinity(manual)\tMOAD\tPDBbind-CN\tBindingDB\n";
    }
    else
    {
        print FP "#pdb\tchain\tresolution\tcsa\tcsa(Renumbered)\tec\tgo\tuniprot\tpubmed\n";
    }
    print FP "$txt";
    close(FP);
    system("gzip -f $rootdir/data/$prefix.tsv");
}

print "generating $rootdir/data/title.tsv\n";
my $txt="";
foreach my $pdbid(`zcat $rootdir/data/pdb_all.tsv.gz |tail -n +2|cut -f1|sort|uniq`)
{
    chomp($pdbid);
    my $divided=substr($pdbid,(length $pdbid)-3,2);
    my $title=`head -1 $rootdir/interim/$divided/$pdbid.txt`;
    chomp($title);
    next if (length $title==0);
    $txt.="$pdbid\t$title\n";
}
open(FP,">$rootdir/data/title.tsv");
print FP $txt;
close(FP);
system("gzip -f $rootdir/data/title.tsv");

print "generation $rootdir/download/lig_frequency.txt\n";
my $today=`date '+%Y-%m-%d'`;
chomp($today);
my $txt=<<EOF
#Please note that "peptide" is named III in old BioLiP; "rna" and "dna" are named "NUC" in old BioLiP
#The following statistics is based on the version $today of BioLiP.

Rank	Lig ID	Frequency
EOF
;
my $r=0;
foreach my $line(`zcat $rootdir/data/lig_all.tsv.gz |cut -f4|sort|uniq -c|sort -nr`)
{
    if ($line=~/(\d+)\s+(\S+)/)
    {
        my $freq    ="$1";
        my $comp_id ="$2";
        $r++;
        $txt.="$r\t$comp_id\t$freq\n";
    }
}
open(FP,">$rootdir/download/lig_frequency.txt");
print FP $txt;
close(FP);


print "generating $rootdir/weekly.html\n";
my @html_items=split(/<!-- CONTENT [A-Z]+ -->/,`cat $rootdir/index.html`);
my $html_head =<<EOF
<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>BioLiP</title>
</head>
<body>

<body bgcolor="#F0FFF0">
<img src=images/BioLiP1.png ></br>
<p><a href=.>[Back to database home]</a></p>
EOF
;


my $html=<<EOF
<div align="justify">
For your convenience, you can use this Perl script for automatic download of all the datasets: <a href=download/download_all_sets.pl>download_all_sets.pl</a>. 
Non-redundant dataset is a subset of the redundant dataset by protein sequence clustering at 90% sequence identity.
For both redundant and non-redundant datasets, the following files are provided:
<br><br>
<li>Receptor: 3D structure of proteins interacting with at least one biologically relevent ligand; original PDB residue numbering.
    Due to hard disk space constraints, BioLiP no longer provides separate batch download for receptor structures with residue number starting from 1. However, we provide a C++ program at <a href=script/receptor1.cpp>receptor1.cpp</a> to facilitate the generation of such files by the user.</li>
<li>Ligand: 3D structures of biologically relevant ligands</li>
<li>Annotations: The annotations for each ligand-protein interaction site, as explained in <a href=download/readme.txt>README</a>.
<br><br>

<style>
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
  text-align: center;
}
</style>

<table border=2 width=100%>
<tr><th width=10%>ID</th><th colspan=3 width=45%>Redundant set</th> <th colspan=3 width=45%>Non-redundant set</th></tr>
EOF
;
#foreach my $divided(`ls $rootdir/weekly/|grep -P "BioLiP_\\w+\\.bsr\\.gz"|cut -f1 -d.|cut -f2 -d_`)
foreach my $divided(`ls $rootdir/weekly/|grep -P "receptor_\\w{2}\\.tar\\.bz2"|cut -f1 -d.|cut -f2 -d_`)
{
    chomp($divided);
    $html.=<<EOF

<tr><td>$divided</td>
    <td><a href=weekly/receptor_$divided.tar.bz2>Receptor</a> </td>
    <td><a href=weekly/ligand_$divided.tar.bz2>Ligand</a> </td>
    <td><a href=weekly/BioLiP_$divided.txt>Annotation</a></td>
    <td><a href=weekly/receptor_${divided}_nr.tar.bz2>Receptor</a> </td>
    <td><a href=weekly/ligand_${divided}_nr.tar.bz2>Ligand</a> </td> 
    <td><a href=weekly/BioLiP_${divided}_nr.txt>Annotation</a></td>
</tr>

EOF
;
}
$html.=<<EOF
</table>
</div>
EOF
;
my $html_tail=<<EOF
</body>
</html>
EOF
;

if (scalar @html_items>=2)
{
    $html_head=$html_items[0];
    $html_tail=$html_items[2];
    if (scalar @html_items>=3)
    {
        $html_tail=$html_items[3];
    }
}

open(FP,">$rootdir/weekly.html");
print FP "$html_head$html$html_tail";
close(FP);

exit();
