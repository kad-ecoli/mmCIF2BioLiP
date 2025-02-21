#!/usr/bin/perl
# parse gene ontology
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

system("$bindir/obo2csv $rootdir/obo/go/go-basic.obo $rootdir/obo/go/is_a.tsv $rootdir/obo/go/go2name.tsv $rootdir/obo/go/alt_id.tsv");
system("zcat $rootdir/data/pdb_all.tsv.gz $rootdir/data/pdb_all.tsv.gz | grep -v '^#'|sort |uniq | cut -f1,2,7 > $rootdir/obo/go/input.txt");
system("$bindir/backpropagate $rootdir/obo/go/input.txt $rootdir/obo/go/is_a.tsv $rootdir/obo/go/alt_id.tsv $rootdir/obo/go/output.txt");
system("sed 's/GO://g' $rootdir/obo/go/output.txt > $rootdir/data/pdb_go.tsv");
&gzipFile("$rootdir/data/pdb_go.tsv");

my @GO_list;
foreach my $GOterm(`cut -f3 $rootdir/obo/go/output.txt|sed 's/,/\\n/g'|sort|uniq`)
{
    chomp($GOterm);
    push(@GO_list,($GOterm));
}
my %GO_dict=map { $_, 1 } @GO_list;


my $txt="";
foreach my $filename(("go2name.tsv","is_a.tsv"))
{
    $txt="";
    foreach my $line(`cat $rootdir/obo/go/$filename`)
    {
        if ($line=~/^(GO:\d+)\t/)
        {
            my $GOterm="$1";
            next if (!exists $GO_dict{$GOterm});
            $txt.="$line";
        }
    }
    open(FP,">$rootdir/data/$filename");
    print FP $txt;
    close(FP);
    &gzipFile("$rootdir/data/$filename");
}

foreach my $filename(("go2name.tsv","is_a.tsv","alt_id.tsv","input.txt","output.txt"))
{
    system("rm $rootdir/obo/go/$filename");
}


my @uniprot_list;
for my $uniprot(`zcat $rootdir/data/pdb_all.tsv.gz |tail -n +2|cut -f8|sort|uniq`)
{
    chomp($uniprot);
    push(@uniprot_list,($uniprot));
}
my %uniprot_dict=map { $_, 1 } @uniprot_list;
$txt="";
for my $line(`zcat $rootdir/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz |grep '^>'`)
{
    chomp($line);
    if ($line=~/^>sp\|(\w+)\|(\w+_\w+\s+[\s\S]+$)/)
    {
        my $uniprot="$1";
        if (exists $uniprot_dict{$uniprot})
        {
            my $name="$2";
            my $GN="";
            if ($name=~/GN=(\S+)/)
            {
                $GN="$1";
            }
            if ($name=~/([\s\S]+?)\s+\w+=/)
            {
                $name="$1";
            }
            $txt.="$uniprot\t$name\t$GN\n";
        }
    }
}
open(FP,">$rootdir/data/uniprot_sprot.tsv");
print FP $txt;
close(FP);
&gzipFile("$rootdir/data/uniprot_sprot.tsv");

print "generating $rootdir/data/index.txt\n";
my $today=`date '+%Y-%m-%d'`;
chomp($today);
my $numProtein  =`zcat $rootdir/data/pdb_all.tsv.gz|wc -l`-1;
my $numRegular  =0;
my $numMetal    =0;
my $numRna      =0;
my $numDna      =0;
my $numPeptide  =0;
my $numBaff     =`zcat $rootdir/data/lig_all.tsv.gz |cut -f9-12|grep -P "\\S+"|wc -l`-1;
my $numManual   =`zcat $rootdir/data/lig_all.tsv.gz |cut -f9 |grep -P "\\S+"|wc -l`-1;
my $numMoad     =`zcat $rootdir/data/lig_all.tsv.gz |cut -f10|grep -P "\\S+"|wc -l`-1;
my $numPdbbind  =`zcat $rootdir/data/lig_all.tsv.gz |cut -f11|grep -P "\\S+"|wc -l`-1;
my $numBindingdb=`zcat $rootdir/data/lig_all.tsv.gz |cut -f12|grep -P "\\S+"|wc -l`-1;
my $numGO       =`zcat $rootdir/data/pdb_go.tsv.gz  |wc -l`+0;
my $numMF       =`zcat $rootdir/data/pdb_go.tsv.gz  |grep 0003674|wc -l`+0;
my $numBP       =`zcat $rootdir/data/pdb_go.tsv.gz  |grep 0008150|wc -l`+0;
my $numCC       =`zcat $rootdir/data/pdb_go.tsv.gz  |grep 0005575|wc -l`+0;
my $numEC       =`zcat $rootdir/data/pdb_all.tsv.gz |cut -f6 |grep -v '^\$'|wc -l`-1;
my @metal_list=();
foreach my $metal(`zcat $rootdir/data/metal.tsv.gz |cut -f1`)
{
    chomp($metal);
    push(@metal_list,($metal));
}
my %metal_dict=map { $_, 1 } @metal_list;
foreach my $ccd(`zcat $rootdir/data/lig_all.tsv.gz|tail -n +2|cut -f4`)
{
    chomp($ccd);
    if ($ccd eq "rna")
    {
        $numRna++;
    }
    elsif ($ccd eq "dna")
    {
        $numDna++;
    }
    elsif ($ccd eq "peptide")
    {
        $numPeptide++;
    }
    elsif (exists $metal_dict{$ccd})
    {
        $numMetal++;
    }
    else
    {
        $numRegular++;
    }
}
my $numEntry=$numRna+$numDna+$numPeptide+$numMetal+$numRegular;
open(FP,">$rootdir/data/index.txt");
print FP <<EOF
<p>
<h1><u>BioLiP in numbers</u></h1>
</p>

BioLiP is updated weekly and the current version ($today) contains:
<li>Number of entries: <a href=qsearch.cgi>$numEntry</a></li>
<li>Number of entries for regular ligands: <a href=qsearch.cgi?lig3=regular>$numRegular</a></li>
<li>Number of entries for metal ligands: <a href=qsearch.cgi?lig3=metal>$numMetal</a></li>
<li>Number of entries for peptide ligands: <a href=qsearch.cgi?lig3=peptide>$numPeptide</a></li>
<li>Number of entries for DNA ligands: <a href=qsearch.cgi?lig3=dna>$numDna</a></li>
<li>Number of entries for RNA ligands: <a href=qsearch.cgi?lig3=rna>$numRna</a></li>
<li>Number of entries with binding affinity data: <a href=qsearch.cgi?baff=baff>$numBaff</a>
(<a href=qsearch.cgi?baff=moad>$numMoad</a> from Binding MOAD,
 <a href=qsearch.cgi?baff=bindingdb>$numBindingdb</a> from BindingDB, and
 <a href=qsearch.cgi?baff=manual>$numManual</a> from manual survey of the original literature)</li>
<li>Number of protein receptors: <a href=qsearch.cgi>$numProtein</a></li>
<li>Number of protein receptors with Enzyme Commission numbers: <a href=qsearch.cgi?ecn=0>$numEC</a></li>
<li>Number of protein receptors with Gene Ontology annotations: <a href=qsearch.cgi?&got=0>$numGO</a>
(<a href=qsearch.cgi?&got=0003674>$numMF</a> with Molecular Function, 
 <a href=qsearch.cgi?&got=0008150>$numBP</a> with Biological Process, and 
 <a href=qsearch.cgi?&got=0005575>$numCC</a> with Cellular Component)</li>
EOF
;
close(FP);

print "generating $rootdir/download/BioLiP*.txt.gz\n";
system("cd $rootdir/weekly; cat `ls BioLiP*_nr.txt`|sort > $rootdir/download/BioLiP_nr.txt");
system("cd $rootdir/weekly; cat `ls $rootdir/weekly/BioLiP*.txt|grep -v _nr.txt`|sort > $rootdir/download/BioLiP.txt");
system("gzip -f $rootdir/download/BioLiP_nr.txt");
system("gzip -f $rootdir/download/BioLiP.txt");

exit(0);

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
