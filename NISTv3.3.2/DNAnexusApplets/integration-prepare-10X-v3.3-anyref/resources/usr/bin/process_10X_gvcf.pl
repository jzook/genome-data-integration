#!/usr/bin/perl -w
#
# process_10X_gvcf.pl - find 10X calls to filter when calls are from a single haplotype
#      - require DP>5 and homozygous with PL for other homozygous genotypes >=99

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use FileHandle;
use English;
use File::Basename;

my $line;



# Check if the right number of parameters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 1) {
    print "$#ARGV arguments\n usage: process_10X_gvcf.pl 10xgvcf.vcf maxcov\n";
    print "example: perl process_10X_gvcf.pl 10xgvcf.vcf 50\n";
    print "Finds 10X calls to filter when calls are from a single haplotype. \n";
    print "Requires DP>5 and DP<maxcov, and homozygous with PL for non-called homozygous genotypes >=99. \n";
     exit;
}



my $infile = "$ARGV[0]";
my $maxcov = "$ARGV[1]";


unless ( open(VCF, "${infile}")) {
	print "\nCannot open the file: ${infile}! \n\n";
	exit;
}

my $infilestart = $infile;
if ( $infile =~ /(.*).vcf/ ) {$infilestart = $1;}
# Output files are opened
unless ( open(OUTPUT, ">${infilestart}_notcallable.bed")) {
    print "\nCannot open the file: ${infilestart}_notcallable.bed to write to! \n\n";
    exit;
}

#Read in information about each callset
my @fields;
#print "@platform\n@dataset\n@callset\n$callsets\n";
#skip header lines in individual files and write first one to output
    my $x=0;
    while ($line = <VCF>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
        if ($x==1) {last;}
    }




my $bedchrom="";
my $bedstart=-1;
my $bedend=-1;
while($line=<VCF>) {
             
	# Split up the line into an array
	@fields = split( "\t", $line);
	
	my $chrom = $fields[0];
	my $pos = $fields[1];
	my $id = $fields[2];
	my $ref = $fields[3];
	my $alt = $fields[4];
	my $qual = $fields[5];
	my $info = $fields[7];
	my $format = $fields[8];
	my @formats = split( ":", $format);
	my $gt="";
	my $gq=-1;
	my $pl;	my @pls;
	my $dp=0;

	my $filt=0;
		 
	#find GT, GQ, DP, and PL for this line to determine if it should be filtered
	my @chars = split( ":", $fields[9]);
	for (my $k=0; $k<=$#chars; $k++) {
		if ($formats[$k] eq "GT") {
			if ($chars[$k] =~ /(.*)\n/) {$gt=$1;}
			else { $gt=$chars[$k];}
			if ($gt =~ /(.*)\|(.*)/) {
				$gt="$1/$2";
			}
			if ($gt =~ /(\d)\/(\d)/) {if ($1>$2) {$gt="$2/$1"; }}
			if ($gt =~ /\./) {$gt="."; $filt=1;} #filter if no GT
			if ($gt eq "1") {$gt="1/1";}
		}
		if ($formats[$k] eq "GQ") {
			if ($chars[$k] =~ /(.*)\n/) {$gq=$1;}
			else { $gq=$chars[$k];}
			if ($gq eq ".") { $gq=-1; }
		}
		if ($formats[$k] eq "PL") {
			if ($chars[$k] =~ /(.*)\n/) {$pl=$1;}
			else { $pl=$chars[$k];}
			
		}
		if ($formats[$k] eq "DP") {
			if ($chars[$k] =~ /(.*)\n/) {$dp=$1;}
			else { $dp=$chars[$k];}
		}
	}
	#if ($gq == 0) { $filt=1; } #filter if GQ == 0, too conservative, so don't do this
	if ($dp < 6 || $dp > $maxcov) { $filt=1; } #filter if coverage < 6 or >max cov specified
	if ($gt eq "0/1" || $gt eq "1/2" || $gt eq "0/2") { $filt=1; } # filter if heterozygous
	#filter if homozygous and any other homozygous PL < 99
	if ($gt eq "0/0") { 
		my @pls = split( ",", $pl);
		if ($#pls>1 && $pls[2]<99) { $filt=1; } #if at least 3 entries and 1/1 PL <99, then filter
		if ($#pls>4 && $pls[5]<99) { $filt=1; } #if at least 6 entries
		if ($#pls>8 && $pls[9]<99) { $filt=1; } #if at least 10 entries
	} elsif ($gt eq "1/1") {
		my @pls = split( ",", $pl);
		if ($#pls>1 && $pls[0]<99) { $filt=1; } #if at least 3 entries
		if ($#pls>4 && $pls[5]<99) { $filt=1; } #if at least 6 entries
		if ($#pls>8 && $pls[9]<99) { $filt=1; } #if at least 10 entries
	} elsif ($gt eq "2/2") {
		my @pls = split( ",", $pl);
		if ($#pls>1 && $pls[0]<99) { $filt=1; } #if at least 3 entries
		if ($#pls>1 && $pls[2]<99) { $filt=1; } #if at least 3 entries
		if ($#pls>8 && $pls[9]<99) { $filt=1; } #if at least 10 entries
	}

	if ($filt==1) {
		my $lineend=$fields[1]+length($fields[3])+50;
		if ($gt eq "0/0") { #for homozygous reference regions, get end coordinate from INFO
			my @infoflds=split( ";", $info);
			for (my $k=0;$k<=$#infoflds;$k++) {
				my @infofld=split( "=", $infoflds[$k]);
				if ($#infofld==1 && $infofld[0] eq "END") { 
					$lineend=$infofld[1]+50;
				}
			}
		}			
		if ($bedend==-1) {
			$bedchrom = $fields[0];
			$bedstart = $fields[1]-1-50; #subtract 1 to create 0-based bed and add 50bp padding
			if ($fields[1]<51) { $bedstart = 0; } #don't allow negative starts
			$bedend = $lineend; #add number of reference bases to uncertain region + 50bp padding
			
		} elsif ($fields[1]-1-50>$bedend || !($fields[0] eq $bedchrom)) { #print previous record if far from current row 
			print OUTPUT "$bedchrom\t$bedstart\t$bedend\n";
			$bedchrom = $fields[0];
			$bedstart = $fields[1]-1-50; #subtract 1 to create 0-based bed and add 50bp padding
			$bedend = $lineend; #add number of reference bases to uncertain region + 50bp padding
		} else { #extend previous record
			if ($bedend < $lineend) {$bedend=$lineend;} #add number of reference bases to uncertain region + 50bp padding
		}
	}
}

print OUTPUT "$bedchrom\t$bedstart\t$bedend\n"; #print last record



close VCF;
close OUTPUT;

