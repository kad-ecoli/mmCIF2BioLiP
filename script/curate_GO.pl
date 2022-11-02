#!/usr/bin/perl
# parse gene ontology
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

system("$bindir/obo2csv $rootdir/obo/go/go-basic.obo $rootdir/obo/go/is_a.tsv $rootdir/obo/go/go2name.tsv $rootdir/obo/go/alt_id.tsv");
system("zcat $rootdir/data/pdb_all.tsv.gz | tail -n +2 | cut -f1,2,7 > $rootdir/obo/go/input.txt");
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
