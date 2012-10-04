#!/usr/bin/perl -w
#
# FosmidSelAllMap.pl - Select fosmids that have less than 2000 bases unmapped to avoid difficult regions of the genome and SVs
#
#
# Version 1.0.0 - May 24, 2012
#

use strict;
use warnings;


my $line;
my $cigar;

# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 0) {
 print "usage: FosmidSelAllMap.pl inputsamfile(without .sam)\n";
 exit;
}
my $infile = $ARGV[0];

# Input file is opened 
unless ( open(INPUT, "${infile}.sam")) {
    print "\nCannot open the file: ${infile}.sam \n\n";
    exit;
}

# Output file is opened 
unless ( open(OUTPUT, ">${infile}_allmap4000.sam")) {
    print "\nCannot open the file: ${infile}_allmap4000.sam to write to! \n\n";
    exit;
}
#print "vcfnt:$vcfcnt\n";


while ($line = <INPUT>) {
    if ($line =~ /^\#/ || $line =~ /^\@/) { #skip header lines
	 print OUTPUT "${line}"; #print header 
        next;
    }
    my @fields;
    # Split up the line into an array 
    @fields = split( "\t", $line);
#    if ($fields[0] =~ /GL/) { next; } #ignore GL... contigs

    $cigar = $fields[5]; #extract cigar field

    my $skipcnt = 0; #number of bases skipped on both ends of fosmid
#print "$info\n";
    if ($cigar =~ /^(\d+)S/) { #find number of bases skipped at beginning
	$skipcnt += $1;
    } 
    if ($cigar =~ /(\d+)S$/) { #find number of bases skipped at end
	$skipcnt += $1;
    } 

    if ($skipcnt < 4001 && $fields[4]>100) { #keep line if <4001 bases skipped and MapQ>100
	 print OUTPUT "${line}"; #print header 
	print "$skipcnt\t";
    }
#last;
}
 
close INPUT;

close OUTPUT;

exit;
