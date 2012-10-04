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

if ($#ARGV != 4) {
 print "$#ARGV arguments\n usage: sqlinsertVcfParams.pl inputcsvfile(without _Hom.txt) platform dataset bam vcf\n";
 exit;
}
my $infile = $ARGV[0];
my $plat = $ARGV[1];
my $dataset = $ARGV[2];
my $bam = $ARGV[3];
my $vcf = $ARGV[4];

# Input file is opened 
unless ( open(INHOM, "${infile}_Hom.csv")) {
    print "\nCannot open the file: ${infile}_Hom.csv to read! \n\n";
    exit;
}
unless ( open(INHET, "${infile}_Het.csv")) {
    print "\nCannot open the file: ${infile}_Het.csv to read! \n\n";
    exit;
}
unless ( open(INOTH, "${infile}_Oth.csv")) {
    print "\nCannot open the file: ${infile}_Oth.csv to read! \n\n";
    exit;
}

# Connect to database file
my $dbh = DBI->connect("dbi:SQLite:dbname=/data/results/justin/SRA/NA12878/NA12878vcfs.db", "", "", 
	{PrintError => 1, RaiseError =>1, AutoCommit => 0});

# Prepare INSERT statements to lookup r3 and f3 reads based on id from SAM file
  my $sthInsAnnot = $dbh->prepare("INSERT INTO vcfhomhet (chrompos, platform, dataset, bam, vcf, geno, alt, qual, gq, pl1, pl2, pl3 ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)");



my $Log_file="${infile}_sqlins.log";
open (OUTLOG,">$Log_file") or die "Cannot open the log file submitted\n";
my @timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
my $i = 0;

while ($line = <INHOM>) {
    if ($line =~ /^((CHROM)|(chrM)|(MT)|(X)|(Y)|(chrX)|(chrY)|(GL))/) { #skip header lines
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


    $sthInsAnnot->execute($chrompos,$plat,$dataset,$dataset,$vcf,3,$fields[2],$fields[3],$fields[4],$fields[5],$fields[6],$fields[7]);
  if (($i % 40000) == 0) {
    $dbh->commit(); #commit inserts into database every 20000 records to save time
#    print $i/40000 . " ";
  }
  $i++;
  #last;
}
    $dbh->commit();
print OUTLOG "Hom done:";
@timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";


while ($line = <INHET>) {
    if ($line =~ /^((CHROM)|(chrM)|(MT)|(X)|(Y)|(chrX)|(chrY)|(GL))/) { #skip header lines
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


    $sthInsAnnot->execute($chrompos,$plat,$dataset,$dataset,$vcf,2,$fields[2],$fields[3],$fields[4],$fields[5],$fields[6],$fields[7]);
  if (($i % 40000) == 0) {
    $dbh->commit(); #commit inserts into database every 20000 records to save time
#    print $i/40000 . " ";
  }
  $i++;
  #last;
}
    $dbh->commit();
print OUTLOG "Het done:";
@timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";


while ($line = <INOTH>) {
    if ($line =~ /^((CHROM)|(chrM)|(MT)|(X)|(Y)|(chrX)|(chrY)|(GL))/) { #skip header lines
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
    $fields[3] = chomp($fields[3]);


    $sthInsAnnot->execute($chrompos,$plat,$dataset,$dataset,$vcf,0,$fields[2],$fields[3],$fields[4],$fields[5],$fields[6],$fields[7]);
  if (($i % 40000) == 0) {
    $dbh->commit(); #commit inserts into database every 20000 records to save time
#    print $i/40000 . " ";
  }
  $i++;
  #last;
}
    $dbh->commit();
print OUTLOG "Oth done:";
@timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";


close OUTLOG; 
close INHOM;
close INHET;
close INOTH;
$dbh->disconnect;

#CREATE TABLE vcfannot (chrompos INTEGER NOT NULL, platform INTEGER NOT NULL, dataset INTEGER NOT NULL, bam INTEGER NOT NULL, vcf INTEGER NOT NULL, geno INTEGER NOT NULL, alt TEXT NOT NULL, qual REAL NOT NULL DEFAULT 0)

