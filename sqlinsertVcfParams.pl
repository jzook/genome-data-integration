#!/usr/bin/perl -w
#
# sqlinsertVcfParams.pl - Convert log10 allele balance probability distributions (from GATK AlleleBalanceDistribution) to normalized decimal probability for forward, reverse, and both strands
#
#
# Version 1.0.0 - Mar 1, 2012
#

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use FileHandle;
use English;
use File::Basename;
#use Try::Tiny;
use DBI;

my $line;
my $line2;

# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 3) {
 print "$#ARGV arguments\n usage: sqlinsertVcfParams.pl inputcsvfile(without _params.csv or _CI.csv) platform dataset bam\n";
 exit;
}
my $infile = $ARGV[0];
my $plat = $ARGV[1];
my $dataset = $ARGV[2];
my $bam = $ARGV[3];

# Input file is opened 
unless ( open(INPARAMS, "${infile}_params.csv")) {
    print "\nCannot open the file: ${infile}_params.csv to read! \n\n";
    exit;
}
unless ( open(INCI, "${infile}_CI.csv")) {
    print "\nCannot open the file: ${infile}_CI.csv to read! \n\n";
    exit;
}

# Connect to database file
my $dbh = DBI->connect("dbi:SQLite:dbname=/data/results/justin/SRA/NA12878/NA12878vcfs.db", "", "", 
	{PrintError => 1, RaiseError =>1, AutoCommit => 0});

# Prepare INSERT statements to lookup r3 and f3 reads based on id from SAM file
  my $sthInsAnnot = $dbh->prepare("INSERT INTO vcfannot (chrompos, platform, dataset, bam, DP, FS, MQ, MQ0, MQ0frac, AB, BaseQRankSum, GC, QD, ReadPosRankSum, MapQRankSum, ACI5, ACI50, ACI95, CCI5, CCI50, CCI95, GCI5, GCI50, GCI95, TCI5, TCI50, TCI95, RefCI5, RefCI50, RefCI95 ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");



my $Log_file="${infile}_sqlins.log";
open (OUTLOG,">$Log_file") or die "Cannot open the log file submitted\n";
my @timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
my $i = 0;

while ($line = <INPARAMS>) {
    if ($line =~ /^((CHROM)|(chrM)|(MT)|(X)|(Y)|(chrX)|(chrY)|(GL))/) { #skip header lines
	$line2=<INCI>;
        next;
    }
    my @fields;
    # Split up the line into an array 
    @fields = split( ",", $line);
    my $chrom;
    if ($fields[0] =~ /chr(.*)/) {
	$chrom = $1;
    } else { $chrom = $fields[0] }
    my $chrompos = $chrom*1e9 + $fields[1];

	$line2=<INCI>;
    my @fields2;
    # Split up the line into an array 
    @fields2 = split( ",", $line2);
    my $endfield=$fields[12];
    if ($fields[12] =~ /(.*)\n/) {  $endfield=$1; }
    my $endfield2=$fields2[14];
    if ($fields2[14] =~ /(.*)\n/) {  $endfield2=$1; }
    my $ref=$fields2[2];
    my $RefCI5=-1;
    my $RefCI50=-1;
    my $RefCI95=-1;
	if ($ref eq "A") { 
	    $RefCI5=$fields2[3];
	    $RefCI50=$fields2[4];
	    $RefCI95=$fields2[5];
	} 
	if ($ref eq "C") { 
	    $RefCI5=$fields2[6];
	    $RefCI50=$fields2[7];
	    $RefCI95=$fields2[8];
	} 
	if ($ref eq "G") { 
	    $RefCI5=$fields2[9];
	    $RefCI50=$fields2[10];
	    $RefCI95=$fields2[11];
	} 
	if ($ref eq "T") { 
	    $RefCI5=$fields2[12];
	    $RefCI50=$fields2[13];
	    $RefCI95=$endfield2;
	} 

    $sthInsAnnot->execute($chrompos,$plat,$dataset,$bam,$fields[2],$fields[3],$fields[4],$fields[5],$fields[6],$fields[7],$fields[8],$fields[9],$fields[10],$fields[11],$endfield,$fields2[3],$fields2[4],$fields2[5],$fields2[6],$fields2[7],$fields2[8],$fields2[9],$fields2[10],$fields2[11],$fields2[12],$fields2[13],$endfield2,$RefCI5,$RefCI50,$RefCI95);
  if (($i % 40000) == 0) {
    $dbh->commit(); #commit inserts into database every 20000 records to save time
#    print $i/40000 . " ";
  }
  $i++;
  #last;
}
    $dbh->commit();
print OUTLOG "done:";
@timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
close OUTLOG; 
close INPARAMS;
close INCI;
$dbh->disconnect;


