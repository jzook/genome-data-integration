#!/usr/bin/perl -w
#
# sqlinsertValidationFosmids.pl - Insert Eichler fosmid data into validation table
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
 print "usage: sqlinsertValidationFosmids.pl inputvcffile(without vcf) dataset#\n";
 exit;
}
my $infile = $ARGV[0];
my $dataset = $ARGV[1];

# Input file is opened 
unless ( open(INPARAMS, "${infile}.vcf")) {
    print "\nCannot open the file: ${infile}.vcf to read! \n\n";
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

# Prepare INSERT statement
  my $sthInsAll = $dbh->prepare("INSERT INTO validation (chrompos, dataset, ACount, CCount, GCount, TCount, geno  ) VALUES (?,?,?,?,?,?,?)");



my $Log_file="${infile}_sqlins.log";
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
    my $chrompos = $fields[0]*1e9 + $fields[1];
    my $ref = $fields[3]; #extract REF field
    my $info = $fields[7]; #extract INFO field
    if ($info =~ /BaseCounts=(.*?);/) { #find BaseCounts
	@BCs = split(",", $1);
    } else { next; }

    my $geno=0; #find genotype
    if ($ref eq "A") {
	if ($BCs[0]>0 && ($BCs[1]+$BCs[2]+$BCs[3]==0)) { $geno=5 #one allele is ref, other unknown
	} elsif ($BCs[0]>0 && ($BCs[1]+$BCs[2]+$BCs[3]>0)) { $geno=2 #HetRef
	} elsif (($BCs[1]+$BCs[2]+$BCs[3]==1)) { $geno=6 ##one allele is non-ref, other unknown
	} elsif (($BCs[1]+$BCs[2]+$BCs[3]==1)) { $geno=4 ##both alleles are non-ref
	}
    } elsif ($ref eq "C") {
	if ($BCs[1]>0 && ($BCs[0]+$BCs[2]+$BCs[3]==0)) { $geno=5 #one allele is ref, other unknown
	} elsif ($BCs[1]>0 && ($BCs[0]+$BCs[2]+$BCs[3]>0)) { $geno=2 #HetRef
	} elsif (($BCs[0]+$BCs[2]+$BCs[3]==1)) { $geno=6 ##one allele is non-ref, other unknown
	} elsif (($BCs[0]+$BCs[2]+$BCs[3]==1)) { $geno=4 ##both alleles are non-ref
	}
    } elsif ($ref eq "G") {
	if ($BCs[2]>0 && ($BCs[1]+$BCs[0]+$BCs[3]==0)) { $geno=5 #one allele is ref, other unknown
	} elsif ($BCs[2]>0 && ($BCs[1]+$BCs[0]+$BCs[3]>0)) { $geno=2 #HetRef
	} elsif (($BCs[1]+$BCs[0]+$BCs[3]==1)) { $geno=6 ##one allele is non-ref, other unknown
	} elsif (($BCs[1]+$BCs[0]+$BCs[3]==1)) { $geno=4 ##both alleles are non-ref
	}
    } elsif ($ref eq "T") {
	if ($BCs[3]>0 && ($BCs[1]+$BCs[2]+$BCs[0]==0)) { $geno=5 #one allele is ref, other unknown
	} elsif ($BCs[3]>0 && ($BCs[1]+$BCs[2]+$BCs[0]>0)) { $geno=2 #HetRef
	} elsif (($BCs[1]+$BCs[2]+$BCs[0]==1)) { $geno=6 ##one allele is non-ref, other unknown
	} elsif (($BCs[1]+$BCs[2]+$BCs[0]==1)) { $geno=4 ##both alleles are non-ref
	}
    } 

    $sthInsAll->execute($chrompos,$dataset,$BCs[0],$BCs[1],$BCs[2],$BCs[3],$geno);
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


