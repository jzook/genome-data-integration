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

# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 0) {
 print "usage: sqlinsertVcfParams.pl inputcsvfile(without _Hom20HetCI.csv)\n";
 exit;
}
my $infile = $ARGV[0];

# Input file is opened 
unless ( open(INPARAMS, "${infile}_Hom20HetCI.csv")) {
    print "\nCannot open the file: ${infile}_Hom20HetCI.csv to read! \n\n";
    exit;
}
#unless ( open(INCI, "${infile}_CI.csv")) {
#    print "\nCannot open the file: ${infile}_CI.csv to read! \n\n";
#    exit;
#}
#print "vcfnt:$vcfcnt\n";

# Connect to database file
my $dbh = DBI->connect("dbi:SQLite:dbname=/data/results/justin/SRA/NA12878/NA12878vcfs.db", "", "", 
	{PrintError => 1, RaiseError =>1, AutoCommit => 0});

# Prepare INSERT statements to lookup r3 and f3 reads based on id from SAM file
  my $sthInsAll = $dbh->prepare("INSERT INTO chrposall (chrompos, geno, ref, AHom,AHet,CHom,CHet,GHom,GHet,THom,THet, ACI5, ACI50, ACI95, CCI5, CCI50, CCI95, GCI5, GCI50, GCI95, TCI5, TCI50, TCI95  ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
#  my $sthInsAnnot = $dbh->prepare("INSERT INTO vcfannot (chrompos, platform, dataset, bam, DP, FS, MQ, MQ0, MQ0frac, AB, BaseQRankSum, GC, QD, ReadPosRankSum, MapQRankSum, ACI5, ACI50, ACI95, CCI5, CCI50, CCI95, GCI5, GCI50, GCI95, TCI5, TCI50, TCI95 ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");



my $Log_file="${infile}_sqlins.log";
open (OUTLOG,">$Log_file") or die "Cannot open the log file submitted\n";
my @timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
my $i = 0;

while ($line = <INPARAMS>) {
    if ($line =~ /^((CHROM)|(chrM)|(MT)|(X)|(Y)|(chrX)|(chrY)|(GL)|(,))/) { #skip header lines and non-chrom 1-22
        #$i++; 
	next;
    }
   if ($i<40001) {$i++; next;}
    my @fields;
    # Split up the line into an array 
    @fields = split( ",", $line);
    my $chrompos = $fields[0]*1e9 + $fields[1];
    my $ref=$fields[2];
    my $geno=0;
    if ((($ref eq "A") && $fields[3]>=0.02) || (($ref eq "C") && $fields[5]>=0.02) || (($ref eq "G") && $fields[7]>=0.02) || (($ref eq "T") && $fields[9]>=0.02)) {
	$geno = 1; #HomRef
    } elsif ((($ref eq "A") && $fields[4]>=1) || (($ref eq "C") && $fields[6]>=1) || (($ref eq "G") && $fields[8]>=1) || (($ref eq "T") && $fields[10]>=1)) {
	$geno = 2; #HetRef
    } elsif ((!($ref eq "A") && $fields[3]>=0.02) || (!($ref eq "C") && $fields[5]>=0.02) || (!($ref eq "G") && $fields[7]>=0.02) || (!($ref eq "T") && $fields[9]>=0.02)) {
	$geno = 3; #HomVar
    }  

    $sthInsAll->execute($chrompos,$geno,$ref,$fields[3],$fields[4],$fields[5],$fields[6],$fields[7],$fields[8],$fields[9],$fields[10],$fields[11],$fields[12],$fields[13],$fields[14],$fields[15],$fields[16],$fields[17],$fields[18],$fields[19],$fields[20],$fields[21],$fields[22]);
  if (($i % 40000) == 0) {
    $dbh->commit(); #commit inserts into database every 20000 records to save time
    print $i/40000 . " ";
  }
  $i++;
  #last;
}
    $dbh->commit();
print OUTLOG "done with F3:";
@timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
close OUTLOG; 
close INPARAMS;
$dbh->disconnect;

exit;
#create indexes in database
 my $SampleCmd = "/home/jzook_local/sqlite3 /data/results/justin/SRA/NA12878/NA12878vcfs.db 'CREATE INDEX chrompos_idx on chrposall (chrompos)'";
 print "$SampleCmd\n";
 
 my $res = system($SampleCmd);   
 if ($res != 0) {
  print "creating index Failed: $res\n";
   exit(1);
 }

 $SampleCmd = "/home/jzook_local/sqlite3 /data/results/justin/SRA/NA12878/NA12878vcfs.db 'CREATE INDEX geno_idx on chrposall (geno)'";
 print "$SampleCmd\n";
 
 $res = system($SampleCmd);   
 if ($res != 0) {
  print "creating index Failed: $res\n";
   exit(1);
 }



