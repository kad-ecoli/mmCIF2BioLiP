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
system("sed 's/GO://g' $rootdir/obo/go/output.txt|gzip - > $rootdir/data/pdb_go.tsv.gz");

my @GO_list;
foreach my $GOterm(`cut -f3 $rootdir/obo/go/output.txt|sed 's/,/\\n/g'|sort|uniq`)
{
    chomp($GOterm);
    push(@GO_list,($GOterm));
}
my %GO_dict=map { $_, 1 } @GO_list;


foreach my $filename(("go2name.tsv","is_a.tsv"))
{
    my $txt="";
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
    system("gzip -f $rootdir/data/$filename");
}

foreach my $filename(("go2name.tsv","is_a.tsv","alt_id.tsv","input.txt","output.txt"))
{
    system("rm $rootdir/obo/go/$filename");
}
exit(0);
