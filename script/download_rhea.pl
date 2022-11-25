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
exit();
