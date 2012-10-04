#!/usr/bin/perl -w
#
# removeMNPs.pl - remove MNPs from vcf
#
#
# Version 1.0.0 - Jul 6, 2012
#

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use FileHandle;
use English;
use File::Basename;
#use Try::Tiny;

my $line;

# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 0) {
 print "usage: sqlinsertVcfParams.pl inputvcffile(without .vcf)\n";
 exit;
}
my $infile = $ARGV[0];

# Input file is opened 
unless ( open(INVCF, "${infile}.vcf")) {
    print "\nCannot open the file: ${infile}.vcf to read! \n\n";
    exit;
}
# Output file is opened 
unless ( open(OUTVCF, ">${infile}_noMNP.vcf")) {
    print "\nCannot open the file: ${infile}_noMNP.vcf to read! \n\n";
    exit;
}

my $i=0;
while ($line = <INVCF>) {
    if ($line =~ /^\#/) {
	print OUTVCF "${line}"; #print header lines to output
	next;
    } 
    if ($line =~ /CHROM/) {
	print OUTVCF "${line}"; #print header lines to output
	next;
    } 
    my @fields;
    # Split up the line into an array 
    @fields = split( "\t", $line);
    if ($fields[3] =~ /../) {$i++; next;} #don't output MNP lines
	print OUTVCF "${line}"; #print header lines to output
}
print "$i records filtered\n";
close OUTVCF; 
close INVCF;
