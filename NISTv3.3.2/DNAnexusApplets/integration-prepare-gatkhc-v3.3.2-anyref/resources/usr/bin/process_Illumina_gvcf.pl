#!/usr/bin/perl -w
#
# process_Illumina_gvcf.pl - find Illumina calls to filter from GATKHC gvcf

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

if ($#ARGV != 0) {
    print "$#ARGV arguments\n usage: process_Illumina_gvcf.pl Illuminagvcf.vcf \n";
    print "example: perl process_Illumina_gvcf.pl Illuminagvcf.vcf \n";
    print "Finds Illumina regions to filter from GATKHC gvcf. \n";
     exit;
}



my $infile = "$ARGV[0]";


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
my $indel=-1;
my $lastfiltvar=-1;
while($line=<VCF>) {
             
	# Split up the line into an array
	@fields = split( "\t", $line);
	
	my $chrom = $fields[0];
	my $pos = $fields[1];
	my $id = $fields[2];
	my $ref = $fields[3];
	my $alt = $fields[4];
	if ($alt eq "<NON_REF>") {
		$alt="";
	} elsif ($alt =~ /(.*),\<NON_REF\>/) {
		$alt=$1;
	}
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
	if ((length($ref)>1 || length($alt)>1) && $gq>60 && !($gt eq "0/0")) {
		if ($bedstart>$pos-60 && $lastfiltvar<$pos-50) { $bedend=-1;} #ignore low confidence reference assertions that are within 10bp of this good indel, since they tend to inaccurately exclude good indels
		$indel=$pos+length($alt)+length($ref)-1;
	}
		

	if ($gq<=60) {
		my $lineend=$fields[1]+length($fields[3])+50;
		my $skiprow=0;
		if ($gt eq "0/0") { #for homozygous reference regions, get end coordinate from INFO
			my @infoflds=split( ";", $info);
			for (my $k=0;$k<=$#infoflds;$k++) {
				my @infofld=split( "=", $infoflds[$k]);
				if ($#infofld==1 && $infofld[0] eq "END") { 
					if ($infofld[1]<$indel+10) {$skiprow=1;} #ignore low confidence reference assertions that are within 10bp of this good indel, since they tend to inaccurately exclude good indels 
					$lineend=$infofld[1]+50;
				}
			}
		} else {
			$lastfiltvar=$pos;
		}	
		if ($skiprow==1) { #don't make filtered bed region for low confidence reference assertions that are within 10bp of this good indel, since they tend to inaccurately exclude good indels 
		} elsif ($bedend==-1) {
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
	#print "$alt;$gq;$indel;$bedstart;$bedend\n";
	#if ($pos>70000) {exit;}
}

print OUTPUT "$bedchrom\t$bedstart\t$bedend\n"; #print last record



close VCF;
close OUTPUT;

