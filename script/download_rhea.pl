#!/usr/bin/perl
# download uniprot accession, taxonomy, EC, GO and pubmed from sifts
# download ncbi taxonomy and scientific name from ncbi
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download FireDB\n";
system("mkdir -p $rootdir/firedb/") if (!-d "$rootdir/firedb/");
foreach my $filename (("cognate","ambiguous","non_cognate"))
{
    system("wget https://firedb.bioinfo.cnio.es/rest/ligand/list/$filename -O $rootdir/firedb/$filename");
}

print "download Rhea\n";
system("mkdir -p $rootdir/rhea/txt/") if (!-d "$rootdir/rhea/txt/");
system("wget https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz -O $rootdir/rhea/txt/rhea-reactions.txt.gz");
if (!-s "$rootdir/rhea/txt/rhea-reactions.txt.gz")
{
    system("wget ftp://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz -O $rootdir/rhea/txt/rhea-reactions.txt.gz");
}
system("mkdir -p $rootdir/rhea/tsv/") if (!-d "$rootdir/rhea/tsv/");
foreach my $filename(("rhea2ec.tsv","rhea2go.tsv","rhea2uniprot_sprot.tsv","rhea2uniprot_trembl.tsv.gz"))
{
    system("wget https://ftp.expasy.org/databases/rhea/tsv/$filename -O $rootdir/rhea/tsv/$filename");
    if (!-s "$rootdir/rhea/tsv/$filename")
    {
        system("wget ftp://ftp.expasy.org/databases/rhea/tsv/$filename -O $rootdir/rhea/tsv/$filename");
    }
}

print "download ChEBI\n";
system("mkdir -p $rootdir/chebi/Flat_file_tab_delimited") if (!-d "$rootdir/chebi/Flat_file_tab_delimited");
system("wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv -O $rootdir/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv");
if (!-s "$rootdir/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv")
{
    system("wget ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv -O $rootdir/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv");
}

print "read $rootdir/rhea/txt/rhea-reactions.txt.gz\n";
my $rhea="";
my %rhea2chebi;
foreach my $line(`zcat $rootdir/rhea/txt/rhea-reactions.txt.gz`)
{
    if ($line=~/^ENTRY\s+RHEA:(\d+)/)
    {
        $rhea="$1";
    }
    elsif ($line=~/^EQUATION/)
    {
        my $chebi_line="";
        foreach my $chebi($line=~/CHEBI:(\d+)/g)
        {
            $chebi_line.=",$chebi" if ($chebi ne "15377"); # exclude water
        }
        $rhea2chebi{$rhea}=$chebi_line;
    }
}

foreach my $term(("ec","go"))
{
    print "read $rootdir/rhea/tsv/rhea2$term.tsv\n";
    my @ec_list;
    my %ec_dict;
    foreach my $line(`cat $rootdir/rhea/tsv/rhea2$term.tsv|tail -n +2|cut -f1,4`)
    {
        if ($line=~/(\d+)\s(\S+)/)
        {
            $rhea="$1";
            my $ec="$2";

            if (!exists $ec_dict{$ec})
            {
                push(@ec_list,($ec));
                $ec_dict{$ec}="$rhea\t".substr($rhea2chebi{$rhea},1);
            }
            elsif ($ec_dict{$ec}=~/(\S+)\t(\S+)/)
            {
                my $old_rhea="$1";
                my $old_chebi="$2";
                $ec_dict{$ec}="$old_rhea,$rhea\t$old_chebi$rhea2chebi{$rhea}";
            }
        }
    }
    my $txt="#$term\trhea\tchebi\n";
    foreach my $ec(@ec_list)
    {
        if ($ec_dict{$ec}=~/(\S+)\t(\S+)/)
        {
            my %chebi_dict=map {$_ => 0} split(/,/,"$2");
            my $chebi_line=join(",",keys %chebi_dict);
            $txt.="$ec\t$1\t$chebi_line\n";
        }
    }
    open(FP,">$rootdir/rhea/${term}2rhea.tsv");
    print FP $txt;
    close(FP);
}

print "read $rootdir/rhea/tsv/rhea2uniprot_*\n";
my @uniprot_list;
foreach my $uniprot(`zcat $rootdir/data/pdb_all.tsv.gz|tail -n +2|cut -f8|sort|uniq`)
{
    chomp($uniprot);
    push(@uniprot_list,($uniprot));
}
my %uniprot_dict=map {$_ => ""} @uniprot_list;

foreach my $line(`cat $rootdir/rhea/tsv/rhea2uniprot_sprot.tsv|tail -n +2|cut -f1,4`)
{
    if ($line=~/(\d+)\s(\w+)/)
    {
        my $rhea="$1";
        my $uniprot="$2";
        if (exists $uniprot_dict{$uniprot})
        {
            $uniprot_dict{$uniprot}.=",$rhea";
        }
    }
}
foreach my $line(`zcat $rootdir/rhea/tsv/rhea2uniprot_trembl.tsv.gz |tail -n +2|cut -f1,4`)
{
    if ($line=~/(\d+)\s(\w+)/)
    {
        my $rhea="$1";
        my $uniprot="$2";
        if (exists $uniprot_dict{$uniprot})
        {
            $uniprot_dict{$uniprot}.=",$rhea";
        }
    }
}

print "write $rootdir/rhea/uniprot2rhea.tsv\n";
my $txt="#uniprot\trhea\tchebi\n";
foreach my $uniprot(@uniprot_list)
{
    next if ((length $uniprot_dict{$uniprot})==0);
    $uniprot_dict{$uniprot}=substr($uniprot_dict{$uniprot},1);
    my $chebi_line;
    foreach my $rhea(split(/,/,$uniprot_dict{$uniprot}))
    {
        $chebi_line.=$rhea2chebi{$rhea};
    }
    my %chebi_dict;
    my @chebi_list;
    foreach my $chebi(split(/,/,substr($chebi_line,1)))
    {
        if (!exists $chebi_dict{$chebi})
        {
            $chebi_dict{$chebi}=1;
            push(@chebi_list,($chebi))
        }
    }
    $chebi_line=join(",",@chebi_list);
    $txt.="$uniprot\t$uniprot_dict{$uniprot}\t$chebi_line\n";
}

open(FP,">$rootdir/rhea/uniprot2rhea.tsv");
print FP $txt;
close(FP);

print "write $rootdir/chebi/chebiId_inchi.tsv\n";
my @chebi_list=();
foreach my $chebi(`grep -v '^#' $rootdir/rhea/ec2rhea.tsv $rootdir/rhea/go2rhea.tsv $rootdir/rhea/uniprot2rhea.tsv |cut -f3|sed 's/,/\\n/g'|sort -n|uniq`)
{
    chomp($chebi);
    push(@chebi_list,($chebi));
}
my %chebi_dict=map {$_ => ""} @chebi_list;

$txt="";
foreach my $line(`cat $rootdir/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv`)
{
    if ($line=~/(\d+)\s+(\S+)/)
    {
        my $chebi="$1";
        my $inchi="$2";
        if (exists $chebi_dict{$chebi})
        {
            $txt.="$chebi\t$inchi\n";
        }
    }
}
open(FP,">$rootdir/chebi/chebiId_inchi.tsv");
print FP $txt;
close(FP);

my @go_list;
foreach my $line(`cat $rootdir/rhea/go2rhea.tsv |cut -f1`)
{
    if ($line=~/GO:(\d+)/)
    {
        push(@go_list,("$1"));
    }
}
my %go_dict=map {$_ => ""} @go_list;
my %pdb2go_dict;
foreach my $line(`zcat $rootdir/data/pdb_go.tsv.gz`)
{
    if ($line=~/(\w+)\s(\w+)\s(\S+)/)
    {
        my $pdbid="$1";
        my $asym_id="$2";
        my @go_list;
        foreach my $go(split(/,/,"$3"))
        {
            if (exists $go_dict{$go})
            {
                push(@go_list,("$go"))
            }
        }
        if (scalar @go_list)
        {
            $pdb2go_dict{"$pdbid:$asym_id"}=join(",",@go_list);
        }
    }
}
my $size=scalar %pdb2go_dict;
print "$size chain with go mapped to rhea\n";
my %ec2rhea;
foreach my $line(`tail -n +2 $rootdir/rhea/ec2rhea.tsv`)
{
    if ($line=~/(\S+)\s(\S+)\s(\S+)/)
    {
        $ec2rhea{"$1"}="$2\t$3";
    }
}
$size=scalar %ec2rhea;
print "$size ec mapping to rhea\n";
my %go2rhea;
foreach my $line(`tail -n +2 $rootdir/rhea/go2rhea.tsv`)
{
    if ($line=~/GO:(\d+)\s(\S+)\s(\S+)/)
    {
        $go2rhea{"$1"}="$2\t$3";
    }
}
$size=scalar %go2rhea;
print "$size go mapping to rhea\n";
my %uniprot2rhea;
foreach my $line(`tail -n +2 $rootdir/rhea/uniprot2rhea.tsv`)
{
    if ($line=~/(\w+)\s(\S+)\s(\S+)/)
    {
        $uniprot2rhea{"$1"}="$2\t$3";
    }
}
$size=scalar %uniprot2rhea;
print "$size uniprot mapping to rhea\n";

print "write $rootdir/data/pdb_rhea.tsv.gz\n";
my $txt="#pdb\tchain\trhea\tchebi\n";
foreach my $line(`zcat $rootdir/data/pdb_all.tsv.gz|tail -n +2`)
{
    my @items=split(/\t/,$line);
    my $pdbid  =$items[0];
    my $asym_id=$items[1];
    my $ec_line=$items[5];
    my $uniprot_line=$items[7];
    my $rhea_line="";
    my $chebi_line="";
    if (length $uniprot_line)
    {
        foreach my $uniprot(split(/,/,$uniprot_line))
        {
            if (exists $uniprot2rhea{$uniprot} && 
                $uniprot2rhea{$uniprot}=~/(\S+)\t(\S+)/)
            {
                if (length $chebi_line)
                {
                    $rhea_line .=",$1";
                    $chebi_line.=",$2";
                }
                else
                {
                    $rhea_line ="$1";
                    $chebi_line="$2";
                }
            }
        }
    }
    if (length $chebi_line==0 && length $ec_line)
    {
        foreach my $ec(split(/,/,$ec_line))
        {
            if (exists $ec2rhea{$ec} && $ec2rhea{$ec}=~/(\S+)\t(\S+)/)
            {
                if (length $chebi_line)
                {
                    $rhea_line .=",$1";
                    $chebi_line.=",$2";
                }
                else
                {
                    $rhea_line ="$1";
                    $chebi_line="$2";
                }
            }
        }
    }
    if (length $chebi_line==0 && exists $pdb2go_dict{"$pdbid:$asym_id"})
    {
        foreach my $go(split(/,/,$pdb2go_dict{"$pdbid:$asym_id"}))
        {
            if ($go2rhea{$go}=~/(\S+)\t(\S+)/)
            {
                if (length $chebi_line)
                {
                    $rhea_line .=",$1";
                    $chebi_line.=",$2";
                }
                else
                {
                    $rhea_line ="$1";
                    $chebi_line="$2";
                }
            }
        }
    }
    if ($chebi_line)
    {
        $txt.="$pdbid\t$asym_id\t$rhea_line\t$chebi_line\n";
    }
}
open(FP,">$rootdir/data/pdb_rhea.tsv");
print FP "$txt";
close(FP);
&gzipFile("$rootdir/data/pdb_rhea.tsv");

my %ccd2inchi;
foreach my $line(`zcat $rootdir/data/ligand.tsv.gz|tail -n +2|cut -f1,3`)
{
    if ($line=~/^(\w+)\t(InChI=\S+)/)
    {
        $ccd2inchi{"$1"}="$2";
    }
}
$size=scalar %ccd2inchi;
print "$size CCD to InChI mapping\n";
my %chebi2inchi;
foreach my $line(`cat $rootdir/chebi/chebiId_inchi.tsv`)
{
    if ($line=~/^(\d+)\t(\S+)/)
    {
        $chebi2inchi{"$1"}="$2";
    }
}
$size=scalar %chebi2inchi;
print "$size ChEBI to InChI mapping\n";
my %firedb_dict;
foreach my $line(`cat $rootdir/firedb/cognate`)
{
    if ($line=~/<compound ID=\"(\w+)\">/)
    {
        $firedb_dict{"$1"}="cognate";
    }
}
foreach my $line(`cat $rootdir/firedb/ambiguous`)
{
    if ($line=~/<compound ID=\"(\w+)\">/)
    {
        $firedb_dict{"$1"}="ambiguous";
    }
}
my %pdb_rhea;
foreach my $line(`zcat $rootdir/data/pdb_rhea.tsv.gz|tail -n +2`)
{
    chomp($line);
    my @items     =split(/\t/,$line);
    my $pdbid     =$items[0];
    my $asym_id   =$items[1];
    my $chebi_line=$items[3];
    $pdb_rhea{"$pdbid:$asym_id"}=$chebi_line;
}
#my %ccd2chebi;
#foreach my $line(`zcat $rootdir/UniChem/src3src7.txt.gz|tail -n +2`)
#{
    #if ($line=~/(\w+)\t(\d+)/)
    #{
        #$ccd2chebi{$1}=$2;
    #}
#}
my %lig_rhea;
if (-s "$rootdir/data/lig_rhea.tsv.gz")
{
    foreach my $line(`zcat $rootdir/data/lig_rhea.tsv.gz`)
    {
        chomp($line);
        if ($line=~/(\w+\t\w+\t\w+)\t(\d)$/)
        {
            $lig_rhea{"$1"}="$2";
        }
    }
}

print "write $rootdir/data/lig_rhea.tsv.gz\n";
my %ccd_chebi2score;
$txt="#pdb\tchain\tCCD\tscore\n";
foreach my $line(`zcat $rootdir/data/lig_all.tsv.gz|tail -n +2|cut -f1,2,4`)
{
    chomp($line);
    my @items  =split(/\t/,$line);
    my $pdbid  =$items[0];
    my $asym_id=$items[1];
    my $ccd    =$items[2];
    my $score  ="non_cognate";
    if (exists $lig_rhea{"$pdbid\t$asym_id\t$ccd"})
    {
        $score=$lig_rhea{"$pdbid\t$asym_id\t$ccd"};
    }
    elsif (!exists $ccd2inchi{$ccd} || !exists $pdb_rhea{"$pdbid:$asym_id"})
    {
        if (exists $firedb_dict{$ccd})
        {
            $score=$firedb_dict{$ccd};
        }
    }
    else
    {
        my $chebi_line=$pdb_rhea{"$pdbid:$asym_id"};
        if (exists $ccd_chebi2score{"$ccd\t$chebi_line"})
        {
            $score=$ccd_chebi2score{"$ccd\t$chebi_line"};
        }
        else
        {
            my $inchi_line=$ccd2inchi{$ccd};
            my $count=0;
            #my $lig="";
            #if (exists $ccd2chebi{$ccd})
            #{
                #$lig=$ccd2chebi{$ccd};
            #}
            foreach my $chebi(split(/,/,$chebi_line))
            {
                #if ($chebi eq $lig)
                #{
                    #$score="CHEBI:$lig";
                    #$count=0;
                    #last;
                #}
                if (exists($chebi2inchi{$chebi}))
                {
                    $inchi_line.=" $chebi2inchi{$chebi}";
                    $count++;
                }
            }
            if ($count)
            {
                my $cmd="echo '$inchi_line'|sed 's/ /\\n/g'|$bindir/obabel -iinchi -ofpt -xfECFP4 -xN1024 -p 7|grep -ohP ' = [.\\d]+'|sort -nr|head -1";
                if (`$cmd`=~/([.\d]+)$/)
                {
                    $score="$1";
                    print "$pdbid\t$asym_id\t$ccd\t$score => ";
                    $score=int(5*$score);
                    if ($score>5)    {$score=5;}
                    elsif ($score<1) {$score=1;}
                    print "$score\n";
                }
                $ccd_chebi2score{"$ccd\t$chebi_line"}=$score;
            }
        }
    }
    $txt.="$pdbid\t$asym_id\t$ccd\t$score\n";
}
open(FP,">$rootdir/data/lig_rhea.tsv");
print FP $txt;
close(FP);
&gzipFile("$rootdir/data/lig_rhea.tsv");

exit();

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
