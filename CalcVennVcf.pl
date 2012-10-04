#!/usr/bin/perl -w
#
# CalcVennVcf.pl - Calculate venn diagram overlaps of vcf calls for 3 methods
#                  from output of GATK VariantsToTable and grep, converted to 
#                  csv in Excel
#
# Version 1.0.0 - Dec 16, 2011
#

use strict;
use warnings;

my ($cnt, $cats);

my $cat1cnt = 0;
my $cat2cnt = 0;
my $cat3cnt = 0;
my $cat12cnt = 0;
my $cat13cnt = 0;
my $cat23cnt = 0;
my $cat123cnt = 0;
my $catothcnt = 0;
my $cat1omnicnt = 0;
my $cat2omnicnt = 0;
my $cat3omnicnt = 0;
my $cat12omnicnt = 0;
my $cat13omnicnt = 0;
my $cat23omnicnt = 0;
my $cat123omnicnt = 0;
my $catothomnicnt = 0;

my $line;

# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 4) {
 print "usage: CalcVennVcf.pl inputfile outputfile cat1 cat2 cat3\n";
 exit;
}

# Assign the two inputs
my $pileupfilename = $ARGV[0];
my $outputfile = $ARGV[1];
my $cat1 = $ARGV[2];
my $cat2 = $ARGV[3];
my $cat3 = $ARGV[4];


# Open the assembly file
unless ( open(PILEUPFILE, $pileupfilename) ) {
    print "\nCould not open file: $pileupfilename!\n\n";
    exit;
}

# Output file is opened (Format ready for further use in R)
unless ( open(OUTPUT, ">$outputfile")) {
    print "\nCannot open the file: $outputfile to write to! \n\n";
    exit;
}

# Write header for output file
print OUTPUT "${cat1}, ${cat2}, ${cat3}, ${cat1}-${cat2}, ${cat1}-${cat3}, ${cat2}-${cat3}, Intersect, Other\n";



while(defined($line = <PILEUPFILE>)) {

    
###############################################################################
# Parse the assembly info

    # Split up the assembly info into an array 
    my @fields = split( ",", $line);
    
    # Get the basic information contig, base number, reference base, coverage
    $cnt = shift @fields;
    $cats = shift @fields;
    $/ = "\r\n"; #make chomp remove windows end of line char
    chomp($cats);
#    print $cats . ";";
    my @catss = split( "-", $cats);
    my $dimcatss = @catss;
    my $catint=0;
    my $catomni=0;
    my $cat1y=0;
    my $cat2y=0;
    my $cat3y=0;
    my $cat1yf=0;
    my $cat2yf=0;
    my $cat3yf=0;

    for(my $i=0; $i<$dimcatss; $i++) {
	if ($catss[$i] eq "Intersection") {
	    $catint=1;
	    $catomni=1;
	} elsif ($catss[$i] eq "OMNI") {
	    $catomni=1;
	} elsif ($catss[$i] eq $cat1) {
	    $cat1y=1;
	} elsif ($catss[$i] eq $cat2) {
	    $cat2y=1;
	} elsif ($catss[$i] eq $cat3) {
	    $cat3y=1;
	} elsif ($catss[$i] eq "filterIn$cat1") {
	    $cat1yf=1;
	} elsif ($catss[$i] eq "filterIn$cat2") {
	    $cat2yf=1;
	} elsif ($catss[$i] eq "filterIn$cat3") {
	    $cat3yf=1;
	} 
#	print $catss[$i] . "," . $cat1y . "," . $cat2y . "," . $cat3y . "\n";
    }

    if ($catint || ($cat1y && $cat2y && $cat3y)) {
	$cat123cnt+=$cnt;
	if ($catomni) { $cat123omnicnt+=$cnt; }
    } elsif ($cat1y && $cat2y) {
	$cat12cnt += $cnt;
	if ($catomni) { $cat12omnicnt+=$cnt; }
    } elsif ($cat1y && $cat3y) {
	$cat13cnt += $cnt;
	if ($catomni) { $cat13omnicnt+=$cnt; }
    } elsif ($cat3y && $cat2y) {
	$cat23cnt += $cnt;
	if ($catomni) { $cat23omnicnt+=$cnt; }
    } elsif ($cat1y) {
	$cat1cnt += $cnt;
	if ($catomni) { $cat1omnicnt+=$cnt; }
    } elsif ($cat2y) {
	$cat2cnt += $cnt;
	if ($catomni) { $cat2omnicnt+=$cnt; }
    } elsif ($cat3y) {
	$cat3cnt += $cnt;
	if ($catomni) { $cat3omnicnt+=$cnt; }
    } else {
	$catothcnt += $cnt;
	if ($catomni) { $catothomnicnt+=$cnt; }
    }
#    exit;
}


    # Print the results in the output file
    printf OUTPUT "%u, %u, %u, %u, %u, %u, %u, %u\n%u, %u, %u, %u, %u, %u, %u, %u\n",
$cat1cnt, 
$cat2cnt, 
$cat3cnt, 
$cat12cnt, 
$cat13cnt, 
$cat23cnt, 
$cat123cnt, 
$catothcnt,        
$cat1omnicnt, 
$cat2omnicnt, 
$cat3omnicnt, 
$cat12omnicnt, 
$cat13omnicnt, 
$cat23omnicnt, 
$cat123omnicnt, 
$catothomnicnt;        
    printf "%u, %u, %u, %u, %u, %u, %u, %u\n%u, %u, %u, %u, %u, %u, %u, %u\n",
$cat1cnt, 
$cat2cnt, 
$cat3cnt, 
$cat12cnt, 
$cat13cnt, 
$cat23cnt, 
$cat123cnt, 
$catothcnt,        
$cat1omnicnt, 
$cat2omnicnt, 
$cat3omnicnt, 
$cat12omnicnt, 
$cat13omnicnt, 
$cat23omnicnt, 
$cat123omnicnt, 
$catothomnicnt;      


close PILEUPFILE;

close OUTPUT;


exit;
