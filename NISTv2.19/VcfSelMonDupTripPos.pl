#!/usr/bin/perl -w
#
# VcfSelMonDupTripPos.pl - separate 1st, 2nd, and 3rd lines into different files at sites with duplicate lines in vcf file since GATK won't call them
#
# Version 1.0.0 - Jan 10, 2013
#

use strict;
use warnings;



# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 0) {
 print "usage: VcfSelDupPos.pl input.vcf (without .vcf)\n";
    print "example: perl VcfSelDupPos.pl /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130114_1a\n";
    exit;
}

# Assign the inputs
my $infile = "$ARGV[0]";


# open files
unless ( open(INUG, "$infile.vcf")) {
    print "\nCannot open the file: $infile.vcf! \n\n";
    exit;
}
unless ( open(OUT1, ">${infile}_1st.vcf")) {
    print "\nCannot open the file: ${infile} to write to! \n\n";
    exit;
}
unless ( open(OUT2, ">${infile}_2nd.vcf")) {
    print "\nCannot open the file: ${infile} to write to! \n\n";
    exit;
}
unless ( open(OUT3, ">${infile}_3rd.vcf")) {
    print "\nCannot open the file: ${infile} to write to! \n\n";
    exit;
}

my $line;

my $x=0;
while ($line = <INUG>) {
    if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, with last one being the line starting with #
    print OUT1 $line;
    print OUT2 $line;
    print OUT3 $line;
    if ($x==1) {last;}
}


my @fields;
my $chromprev="";
my $posprev=0;
my $chromprev2="";
my $posprev2=0;
my $trip=0;
while($line=<INUG>) {

    
    # Split up the line into an array
    @fields = split( "\t", $line);

    if ($fields[0] eq $chromprev2 && $fields[1]==$posprev2) {
		if ($trip==0) {
			print OUT3 "$line";
		} else {print "$line"; } #four records with same position!
    	$trip=1;
    } elsif ($fields[0] eq $chromprev && $fields[1]==$posprev) {
		print OUT2 "$line";
    } else {
		print OUT1 "$line";
		$trip=0;
    }
    $chromprev2=$chromprev;
    $posprev2=$posprev;
    $chromprev=$fields[0];
    $posprev=$fields[1];

}


close OUT1;
close OUT2;
close OUT3;
close INUG;

