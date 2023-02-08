#!/usr/bin/perl -w
use strict;


`mkdir -p BioLiP_updated_set`;
#`rm -fr BioLiP_updated_set/*`;
chdir "BioLiP_updated_set";

my $head= "https://zhanggroup.org/BioLiP2/weekly";
my $address="$head.html";
system("wget -o log -c $address") == 0 or die "System call failed: $!";


my @rst=`cat weekly.html`;
my @all=();
foreach my $r(@rst)
{
    if($r =~ /\<tr\>\<td\>(\S+)\<\/td\>/)
    {
	#print "$1\n";
	push(@all, $1);
    }
}
my $tot=@all;
print "\n====================================================\n";
print "In total, there are $tot weeks to update.\n\n";

my $annotation="BioLiP_UP.txt";
open(OUT, ">$annotation");
close(OUT);
my $annotation1="BioLiP_UP_nr.txt";
open(OUT, ">$annotation1");
close(OUT);

foreach my $r(@all)
{        
    my $rec="receptor_$r.tar.bz2";
    my $lig="ligand_$r.tar.bz2";

	print "Dowload redundant set for the week $r...\n";    
	system("wget -o log -c $head/$rec") == 0 or die "System call failed: $!";
	system("tar -xvf $rec >log")== 0 or die "System call failed: $!";
	system("wget -o log -c $head/$lig") == 0 or die "System call failed: $!";
	system("tar -xvf $lig >log")== 0 or die "System call failed: $!";
	
	
	print "Dowload non-redundant set for the week $r...\n";
	$rec="receptor_$r\_nr.tar.bz2";
	$lig="ligand_$r\_nr.tar.bz2";
	system("wget -o log -c $head/$rec") == 0 or die "System call failed: $!";
	system("tar -xvf $rec >log")== 0 or die "System call failed: $!";
	system("wget -o log -c $head/$lig") == 0 or die "System call failed: $!";
	system("tar -xvf $lig >log")== 0 or die "System call failed: $!";
}


print "Cheers! All updates are done.\n";
print "====================================================\n\n";
";
