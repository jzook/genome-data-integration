#!/usr/bin/perl -w
#
# sqlselectVcf.pl - create vcfs from consensus genotypes in database
#
# Version 1.0.0 - Jul 31, 2012
# Version 1.2 - Aug 28, 2012 - select genotypes based on v1.2, allowing 1/6 datasets to disagree, 
# 								and not filtering by mapping bias for hom calls
#
#
print "test";
use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use FileHandle;
use English;
use File::Basename;
#use Try::Tiny;
use DBI;
use POSIX;

my $line;
my $line2;

# Check if the right number of parameters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 0) {
 print "$#ARGV arguments\n usage: sqlselectVcf_v1.2.pl vcffilestart\n";
 exit;
}
my $filen = $ARGV[0];

# Connect to database file
my $dbh = DBI->connect("dbi:SQLite:dbname=/Volumes/SSD960/workspace/data/NA12878vcfCompare_120315/NA12878vcfs.db", "", "",
	{PrintError => 1, RaiseError =>1, AutoCommit => 0});


# Prepare SELECT statement to create vcf file based on genotypes

  my $sthSelChrPos = $dbh->prepare("SELECT chrompos, genoAll,genoSeqV,(GenoMapGood+GenoMapGoodSeqV),TrancheSSEmin2,TrancheABQDmin2,TrancheAlignmin2,TrancheMapmin2,ref,alt FROM chrposall ORDER BY chrompos");

open (VARONLY,">${filen}_VarOnly.vcf") or die "Cannot open ${filen}_VarOnly.vcf\n";
open (VARUNCERT,">${filen}_VarUncert.vcf") or die "Cannot open ${filen}_VarUncert.vcf\n";
open (ALL,">${filen}_All.vcf") or die "Cannot open ${filen}_All.vcf\n";


my $Log_file="sqlselectVcf.log";
open (OUTLOG,">$Log_file") or die "Cannot open the log file submitted\n";
my @timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
my $i = 0;



  $sthSelChrPos->execute();
my $jj=-1; #ignore datasets>$bams since one dataset is in twice
while (my $row = $sthSelChrPos->fetchrow_arrayref) {
    my $chr = floor(@$row[0] / 1000000000);
    my $pos = @$row[0] % 1000000000;
    my $ref = @$row[8];
    my $alt = @$row[9];
    if ($alt eq "N") { next; }
    my $filt="";
    my $geno=0;
    my $datasets=0;
    my $arbitration=0;
    my $gt="0/0";
    if ((@$row[1]+@$row[2]) % 10 == 1) { #HomRef
        $geno=1; $gt="0/0";
        if (@$row[3]<3) {
#            $filt = $filt . "Map3";
        } 
        if ((@$row[1]+@$row[2]) % 1000 <30) {
            $filt = $filt . "Twodatasets";
        } 
        if (@$row[6]>=99) {
            $filt = $filt . "TrancheAlign";
        } 
        if (@$row[7]>=99) {
#            $filt = $filt . "TrancheMap";
        }
    } elsif ((@$row[1]+@$row[2]) % 10 == 3) { #HomVar
        $geno=3; $gt="1/1";
        if (@$row[3]<3) {
#            $filt = $filt . "Map3";
        } 
        if ((@$row[1]+@$row[2]) % 1000 <30) {
            $filt = $filt . "Twodatasets";
        } 
        if (@$row[5]>=99) {
            $filt = $filt . "TrancheABQD";
        } 
        if (@$row[6]>=99) {
            $filt = $filt . "TrancheAlign";
        } 
        if (@$row[7]>=99) {
#            $filt = $filt . "TrancheMap";
        }
    } elsif ((@$row[1]+@$row[2]) % 10 == 2) { #Het
        $geno=2; $gt="0/1";
        if (@$row[3]<3) {
            $filt = $filt . "Map3";
        } 
        if (@$row[4]>=95) {
            $filt = $filt . "TrancheSSE";
        } 
        if (@$row[5]>=95) {
            $filt = $filt . "TrancheABQD";
        } 
        if (@$row[6]>=95) {
            $filt = $filt . "TrancheAlign";
        } 
        if (@$row[7]>=95) {
            $filt = $filt . "TrancheMap";
        }
    } elsif ((@$row[1]+@$row[2]) % 10 == 0) { #unfiltered uncertain
         $filt="Uncertain";
    }
    my $qual=0;
    if ($filt eq "") {
        $filt = "PASS";
        if ($geno>1) { $qual=100; } 
    } else {
        $filt = "filtered" . $filt;
        $qual=10;
    }

    if ($geno>0) {
        $datasets = floor(((@$row[1]+@$row[2]) % 1000)/10);
        if (@$row[1]>0) { $arbitration=@$row[1]; }
        if (@$row[2]>0) { $arbitration=@$row[2]+10000; }
    }
    if ($geno>1) {
        print VARONLY "$chr\t$pos\t.\t$ref\t$alt\t$qual\t$filt\tdatasets=$datasets;arbitration=$arbitration\tGT\t$gt\n";
    }
    if ($geno!=1 || !($filt eq "PASS")) {
        print VARUNCERT "$chr\t$pos\t.\t$ref\t$alt\t$qual\t$filt\tdatasets=$datasets;arbitration=$arbitration\tGT\t$gt\n";
    }
        print ALL "$chr\t$pos\t.\t$ref\t$alt\t$qual\t$filt\tdatasets=$datasets;arbitration=$arbitration\tGT\t$gt\n";


#  last;
}

#run updates for last record


    $dbh->commit();
print OUTLOG "done:";
@timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
close OUTLOG;
close VARONLY;
close VARUNCERT;
close ALL;

$dbh->disconnect;


