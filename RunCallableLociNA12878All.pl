#!/usr/bin/perl -w
#
# RunCallableLociNA12878All.pl - Run GATK Callable Loci on all NA12878 bam files using NA12878datasets.csv
#
#
# Version 1.0.0 - Sep 19, 2012
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

if ($#ARGV != 2) {
 print "$#ARGV arguments\n usage: RunCallableLociNA12878All.pl inputdatasetscsvfile beginline endline\n";
 exit;
}
my $infile = $ARGV[0];
my $line1 = $ARGV[1];
my $line2 = $ARGV[2];

# Input file is opened 
unless ( open(IN, "${infile}")) {
    print "\nCannot open the file: ${infile} to read! \n\n";
    exit;
}
my $i = 0;
my $res=0;

while ($line = <IN>) {
    $i++;
    if ($i<$line1 || $i>$line2) {next;}
    my @fields;
    # Split up the line into an array 
    @fields = split( ",", $line);
     if ($res == 0) {
	 my $samCmd = "java -jar -Xmx2g  ~/bioinfo/GenomeAnalysisTK-2.1-8-g5efb575/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T CallableLoci -I $fields[0] -summary $fields[0].callable.summary -o $fields[0].callable.bed -frlmq 1 -mlmq 0 -mbq 10 -minDepth 5 -mmq 0 &";
	  print "$samCmd\n";
	 $res = system($samCmd);
    }   

}
