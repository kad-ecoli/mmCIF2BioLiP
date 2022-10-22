#!/usr/bin/perl
my $docstring=<<EOF
download_pubmed.pl pmid
    download pubmed abstract to pmid.txt
EOF
;

use strict;
if (scalar @ARGV==0)
{
    print "$docstring\n";
    exit(0);
}

for (my $a=0;$a<scalar @ARGV;$a++)
{
    my $pmid=$ARGV[$a];
    next if (-s "$pmid.txt");
    system("wget -q 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$pmid&retmode=abstract&rettype=text' -O $pmid.txt");
    next if (-s "$pmid.txt");
    system("wget -q 'https://pubmed.ncbi.nlm.nih.gov/$pmid/?format=pubmed' -O $pmid.html");
    my $txt="";
    my $pre=0;
    foreach my $line(`cat $pmid.html`)
    {
        if ($line=~/^\s*<pre/) { $pre=1; }
        elsif ($line=~/^\s*<\/pre>/) { $pre=2; }
        elsif ($pre==1) { $txt.="$line"; }
    }
    open(FP,">$pmid.txt");
    print FP "$txt";
    close(FP);
    system("rm $pmid.html");
}
