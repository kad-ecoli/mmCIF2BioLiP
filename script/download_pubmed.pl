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
    #next if (-s "$pmid.txt");
    system("wget -q 'https://pubmed.ncbi.nlm.nih.gov/$pmid/?format=pubmed' -O $pmid.html");
    my $abstract_txt="";
    my $title_txt="";
    my $ab=0;
    my $ti=0;
    foreach my $line(`cat $pmid.html`)
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
    open(FP,">$pmid.txt");
    print FP "$title_txt\n$abstract_txt\n";
    close(FP);
    system("rm $pmid.html");
}
