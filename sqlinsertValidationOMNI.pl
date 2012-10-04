#!/usr/bin/perl -w
#
# sqlinsertValidationOMNI.pl - Insert OMNI data into validation table
#
#
# Version 1.0.0 - May 26, 2012
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

if ($#ARGV != 1) {
 print "usage: sqlinsertValidationOMNI.pl inputvcffile(without vcf) dataset\#\n";
 exit;
}
my $infile = $ARGV[0];
my $dataset = $ARGV[1];

# Input file is opened 
unless ( open(INPARAMS, "${infile}.vcf")) {
    print "\nCannot open the file: ${infile}.vcf to read! \n\n";
    exit;
}

# Connect to database file
my $dbh = DBI->connect("dbi:SQLite:dbname=/data/results/justin/SRA/NA12878/NA12878vcfs.db", "", "", 
	{PrintError => 1, RaiseError =>1, AutoCommit => 0});

# Prepare INSERT statement
  my $sthInsAll = $dbh->prepare("INSERT INTO validation (chrompos, dataset, ACount, CCount, GCount, TCount, geno  ) VALUES (?,?,?,?,?,?,?)");



my $Log_file="${infile}_omni_sqlins.log";
open (OUTLOG,">$Log_file") or die "Cannot open the log file submitted\n";
my @timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
my $i = 0;

while ($line = <INPARAMS>) {
    if ($line =~ /^((\#)|(CHROM)|(chrM)|(MT)|(X)|(Y)|(chrX)|(chrY)|(GL)|(,))/) { #skip header lines and non-chrom 1-22
        #$i++; 
	next;
    }
    my @fields;
    my @BCs;
    # Split up the line into an array 
    @fields = split( "\t", $line);
    if (!($fields[6] eq "PASS") && !($fields[6] eq ".")) { next; } #only include variants that passed the filter
    my $chrompos = $fields[0]*1e9 + $fields[1];
    my $ref = $fields[3]; #extract REF field
    my $info = $fields[7]; #extract INFO field

    my $genono=0; # genotype number for table
    my $genoform = $fields[8]; #extract genotype format field
    my $geno = $fields[9]; #extract genotype field
    my $alt =$fields[4]; #extract ALT field


#print "\n$info\n";
#last;
    my @genoformf = split( ":", $genoform);
    my @genof = split( ":", $geno);
    my $genotype; my $GQ=0; my @PL=((0) x 3); my $alt1; my $alt2;
    my $i=0;
    my $hetnoref = 0;
	    if ($alt =~ /^(.),(.)/) {$alt1=$1; $alt2=$2; $hetnoref=1; #Het non-ref SNP
	    } elsif ($alt =~ /^../) {next; #skip non-SNP genotypes
	    } else {$alt1=$alt;}
    while ($i<=$#genoformf) {
	if ($genoformf[$i] =~ /GT/) {$genotype=$genof[$i];}
	if ($genoformf[$i] =~ /GQ/) {if ($genof[$i] =~ /(.*)/) {$GQ=$1;}} #remove newline
	if ($genoformf[$i] =~ /PL/) {
	    if ($hetnoref==0) {if ($genof[$i] =~ /^(.*?),(.*?),(.*)/) {$PL[0]=$1;$PL[1]=$2;$PL[2]=$3;}}
	    if ($hetnoref==1) {if ($genof[$i] =~ /^(.*?),(.*?),(.*?),/) {$PL[0]=$1;$PL[1]=$2;$PL[2]=$3;}}
	}
	$i++;
    }
    my @BCs = ((0) x 4);
    if ($alt1 eq "A") { $BCs[0]=1;
    } elsif ($alt1 eq "C") { $BCs[1]=1;
    } elsif ($alt1 eq "G") { $BCs[2]=1;
    } elsif ($alt1 eq "T") { $BCs[3]=1; }
    if ($hetnoref==1) { 
    if ($alt2 eq "A") { $BCs[0]=1;
    } elsif ($alt2 eq "C") { $BCs[1]=1;
    } elsif ($alt2 eq "G") { $BCs[2]=1;
    } elsif ($alt2 eq "T") { $BCs[3]=1; }
	$sthInsAll->execute($chrompos,$dataset,$BCs[0],$BCs[1],$BCs[2],$BCs[3],4);
    } elsif ($genotype =~ /^((0(\||\/)1)|(1(\||\/)0))/) { #if genotype is 0/1 or 0|1 or 1/0 or 1|0 then it's heterozygous
	$sthInsAll->execute($chrompos,$dataset,$BCs[0],$BCs[1],$BCs[2],$BCs[3],2);
    } elsif ($genotype =~ /^1(\||\/)1/) { #if genotype is 1/1 or 1|1 then it's homozygous variant
	$sthInsAll->execute($chrompos,$dataset,2*$BCs[0],2*$BCs[1],2*$BCs[2],2*$BCs[3],3);
    } elsif ($genotype =~ /^0(\||\/)0/) { #if genotype is 0/0 or 0|0 then it's homozygous ref
	$sthInsAll->execute($chrompos,$dataset,2*$BCs[0],2*$BCs[1],2*$BCs[2],2*$BCs[3],1);
    } else { 
	$sthInsAll->execute($chrompos,$dataset,$BCs[0],$BCs[1],$BCs[2],$BCs[3],4);
    }


  if (($i % 40000) == 0) {
    $dbh->commit(); #commit inserts into database every 40000 records to save time
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


