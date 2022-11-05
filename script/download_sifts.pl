#!/usr/bin/perl
# download uniprot accession, taxonomy, EC, GO and pubmed from sifts
# download ncbi taxonomy and scientific name from ncbi
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download go\n";
system("mkdir -p $rootdir/obo/go/") if (!-d "$rootdir/obo/go/");
system("wget http://purl.obolibrary.org/obo/go/go-basic.obo -O $rootdir/obo/go/go-basic.obo");
system("wget http://current.geneontology.org/ontology/go-basic.obo -O $rootdir/obo/go/go-basic.obo") if (!-s "$rootdir/obo/go/go-basic.obo");

print "download swissprot\n";
system("mkdir -p $rootdir/uniprot/current_release/knowledgebase/complete") if (!-d "$rootdir/uniprot/current_release/knowledgebase/complete");
system("wget http://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O $rootdir/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz");
system("wget  ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O $rootdir/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz") if (!-s "$rootdir/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz");

print "download sifts\n";
system("mkdir -p $rootdir/sifts/flatfiles/tsv");
foreach my $filename(("pdb_chain_uniprot.tsv.gz",
                      "pdb_chain_taxonomy.tsv.gz",
                      "pdb_pubmed.tsv.gz",
                      "pdb_chain_enzyme.tsv.gz",
                      "pdb_chain_go.tsv.gz"))
{
    system("wget http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/$filename -O $rootdir/sifts/flatfiles/tsv/$filename");
    system("wget  ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/$filename -O $rootdir/sifts/flatfiles/tsv/$filename") if (!-s "$rootdir/sifts/flatfiles/tsv/$filename");
}

system("mkdir -p $rootdir/taxonomy");
system("wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O $rootdir/taxonomy/taxdump.tar.gz");
system("wget  ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O $rootdir/taxonomy/taxdump.tar.gz") if (!-s "$rootdir/taxonomy/taxdump.tar.gz");
system("cd $rootdir/taxonomy; tar -xvf taxdump.tar.gz names.dmp");
if (-s "$rootdir/taxonomy/names.dmp")
{
    system("rm $rootdir/taxonomy/taxdump.tar.gz");
    &gzipFile("$rootdir/taxonomy/names.dmp");
}

system("mkdir -p $rootdir/data");
print "parsing taxonomy\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz")
{
    system("zcat $rootdir/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz |tail -n +3|cut -f1-3|uniq > $rootdir/data/chain2taxonomy.tsv");
}
my @taxid_list;
foreach my $taxid(`cat $rootdir/data/chain2taxonomy.tsv|cut -f3|sort|uniq`)
{
    chomp($taxid);
    push(@taxid_list,($taxid));
}

my %taxid_dict = map { $_ => "" } @taxid_list;
foreach my $line(`zcat $rootdir/taxonomy/names.dmp.gz`)
{
    my @items=split(/\t\|\t/,$line);
    my $taxid=$items[0];
    next if (!exists $taxid_dict{$taxid});
    my $name_type=$items[-1];
    next if (length $taxid_dict{$taxid} && $name_type!~/scientific name/);
    $taxid_dict{$taxid}=$items[1];
}
my $txt="";
foreach my $taxid(@taxid_list)
{
    $txt.="$taxid\t$taxid_dict{$taxid}\n";
}
if (scalar @taxid_list)
{
    open(FP,">$rootdir/data/taxid2name.tsv");
    print FP $txt;
    close(FP);
    &gzipFile("$rootdir/data/taxid2name.tsv");
}
&concatenate_sift("cat $rootdir/data/chain2taxonomy.tsv","$rootdir/data/chain2taxonomy.tsv");

print "parsing uniprot accession\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz")
{
    &concatenate_sift("zcat $rootdir/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz |tail -n +3|cut -f1-3|sort|uniq","$rootdir/data/chain2uniprot.tsv");
}

print "parsing EC\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_chain_enzyme.tsv.gz")
{
    &concatenate_sift("zcat $rootdir/sifts/flatfiles/tsv/pdb_chain_enzyme.tsv.gz |tail -n +3|cut -f1,2,4|sort|uniq","$rootdir/data/chain2ec.tsv");
}

print "parsing GO\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_chain_go.tsv.gz")
{
    &concatenate_sift("zcat $rootdir/sifts/flatfiles/tsv/pdb_chain_go.tsv.gz |tail -n +3|cut -f1,2,6|sort|uniq","$rootdir/data/chain2go.tsv");
}

print "parsing pubmed\n";
if (-s "$rootdir/sifts/flatfiles/tsv/pdb_pubmed.tsv.gz")
{
    system("zcat $rootdir/sifts/flatfiles/tsv/pdb_pubmed.tsv.gz |tail -n +3|grep -P '\\t0\\t'|cut -f1,3|sort|uniq > $rootdir/data/pdb2pubmed.tsv");
    &gzipFile("$rootdir/data/pdb2pubmed.tsv");
}

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

sub concatenate_sift
{
    my ($cmd,$outfile)=@_;
    my @chain_list;
    my %chain_dict;
    foreach my $line(`$cmd`)
    {
        chomp($line);
        my @items=split("\t",$line);
        my $chain="$items[0]\t$items[1]";
        if (exists $chain_dict{$chain})
        {
            $chain_dict{$chain}.=",$items[2]";
        }
        else
        {
            push(@chain_list,($chain));
            $chain_dict{$chain}="$items[2]";
        }
    }
    my $txt="";
    foreach my $chain(@chain_list)
    {
        $txt.="$chain\t$chain_dict{$chain}\n";
    }
    open(FP,">$outfile");
    print FP $txt;
    close(FP);
    &gzipFile($outfile);
}
