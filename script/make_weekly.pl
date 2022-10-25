#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

system("mkdir -p $rootdir/weekly") if (!-d "$rootdir/weekly");
foreach my $divided(`ls $rootdir/interim/`)
{
    chomp($divided);
    next if (`ls $rootdir/interim/$divided/|grep bsr|wc -l`+0==0);
    print "prepare BioLiP_$divided.txt ligand_$divided.tar.bz2 receptor_$divided.tar.bz2\n";
    foreach my $moltype(("protein","peptide","rna","dna"))
    {
        system("cat $rootdir/interim/$divided/*bsr| grep -A1 '^>' --no-group-separator |grep -PA1 '\\t$moltype\\t' --no-group-separator |gzip - > $rootdir/weekly/${moltype}_$divided.fasta.gz");
    }
    system("rm -rf $rootdir/weekly/$divided/") if (-d "$rootdir/weekly/$divided");
    system("mkdir -p $rootdir/weekly/$divided/receptor/");
    system("mkdir -p $rootdir/weekly/$divided/receptor1/");
    system("mkdir -p $rootdir/weekly/$divided/ligand/");
    foreach my $filename(`ls $rootdir/interim/$divided/|grep .tar.bz2`)
    {
        chomp($filename);
        system("cd $rootdir/weekly/$divided; tar -xf $rootdir/interim/$divided/$filename; mv *_*_*.pdb $rootdir/weekly/$divided/ligand/; mv *.pdb $rootdir/weekly/$divided/receptor/");
    }
    system("ls $rootdir/weekly/$divided/receptor/ | $bindir/receptor1 $rootdir/weekly/$divided/receptor/ $rootdir/weekly/$divided/receptor1/ -");
    foreach my $moltype(("receptor","receptor1","ligand"))
    {
        print "$rootdir/weekly/${moltype}_$divided.tar.bz2\n";
        system("cd $rootdir/weekly/$divided; tar -cjvf $rootdir/weekly/${moltype}_$divided.tar.bz2 $moltype/");
    }
    system("cd $rootdir/interim/$divided/; sed -n '/^#pdb/{ :a; n; p; ba; }' *.bsr|gzip - > $rootdir/weekly/BioLiP_$divided.bsr.gz");
}

print "combine data\n";
system("mkdir -p $rootdir/data") if (!-d "$rootdir/data");
foreach my $moltype(("protein","peptide","rna","dna"))
{
    system("cd $rootdir/weekly; zcat ${moltype}_*.fasta.gz > $rootdir/data/$moltype.fasta");
    if ($moltype eq "protein")
    {
        system("$bindir/cd-hit -i $rootdir/data/$moltype.fasta -o $rootdir/data/${moltype}_nr.fasta");
    }
    else
    {
        system("$bindir/fasta2nr $rootdir/data/$moltype.fasta $rootdir/data/${moltype}_nr.fasta");
    }
    system("gzip -f $rootdir/data/$moltype.fasta");
    system("gzip -f $rootdir/data/${moltype}_nr.fasta");
    system("cd $rootdir/weekly; rm ${moltype}_*.fasta.gz");
}

exit();
