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

my %csa_dict;
foreach my $line(`cat $rootdir/data/csa.tsv`)
{
    if ($line=~/^(\S+)\t(\S+)\t(\S+)/)
    {
        my $pdbid  ="$1";
        my $asym_id="$2";
        my @res_list=split(/,/,"$3");

        my $divided=substr($pdbid,(length $pdbid)-3,2);
        my $filename="$rootdir/weekly/$divided/receptor/$pdbid$asym_id.pdb";
           $filename="$rootdir/weekly/$divided/receptor_nr/$pdbid$asym_id.pdb" if (!-s "$filename");
        if (!-s "$filename")
        {
            print "no pdb file for $pdbid$asym_id from csa\n";
            next;
        }
        my $csaOrig="";
        my $csaRenu="";
    }
}
$size=keys %csa_dict;
print "$size csa site\n";

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
        my $bindcn =""; # [16] binding affinity from PDBbind-CN
        my $binddb =""; # [17] binding affinity from BindingDB
        my $uniprot=""; # [18] uniprot accession
           $uniprot=$uniprot_dict{$chain} if (exists $uniprot_dict{$chain});
        my $pubmed =""; # [19] pubmed id
           $pubmed =$pubmed_dict{$pdbid} if (exists $uniprot_dict{$pdbid});
        my $sequence=$fasta_dict{$chain};#[20] receptor sequence
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
    system("mkdir -p $rootdir/weekly/$divided/receptor_nr");
    system("mkdir -p $rootdir/weekly/$divided/receptor_nr1");
    system("mkdir -p $rootdir/weekly/$divided/ligand_nr");
    my %filename_dict=map {$_ => 1} @receptor_nr;
    my @filename_list=keys %filename_dict;
    foreach my $filename(@filename_list)
    {
        system("mv $rootdir/weekly/$divided/receptor/$filename $rootdir/weekly/$divided/receptor_nr/$filename");
        system("mv $rootdir/weekly/$divided/receptor1/$filename $rootdir/weekly/$divided/receptor_nr1/$filename");
    }
    %filename_dict=map {$_ => 1} @ligand_nr;
    @filename_list=keys %filename_dict;
    foreach my $filename(@filename_list)
    {
        system("mv $rootdir/weekly/$divided/ligand/$filename $rootdir/weekly/$divided/ligand_nr/$filename");
    }
    system("cd $rootdir/weekly/$divided/; tar -cjvf $rootdir/weekly/receptor1_${divided}_nr.tar.bz2 receptor_nr1/");
    system("cd $rootdir/weekly/$divided/; tar -cjvf $rootdir/weekly/receptor_${divided}_nr.tar.bz2 receptor_nr/");
    system("cd $rootdir/weekly/$divided/; tar -cjvf $rootdir/weekly/ligand_${divided}_nr.tar.bz2 ligand_nr/");
}


exit();
