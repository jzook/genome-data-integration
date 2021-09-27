#!/usr/bin/perl
use strict;

######################################################
# vcfsorter.pl
#
# Copyright (C) 2011 German Gaston Leparc
#
# sorts VCF by reference genome
#
# usage:
#
# vcfsorter.pl genome.dict myvcf.file > mynewvcf.file
#
######################################################

my $usage = <<EOF;
sorts VCF by reference genome

usage:

vcfsorter.pl genome.dict myvcf > mynewvcf.file 2>STDERR
EOF


my $dict_file = @ARGV[0];
my $vcf_file = @ARGV[1];

die "\nERROR: missing an argument!\n\n$usage" if (@ARGV < 2);


#---------------------------------------- LOAD IN FASTA DICT INTO MEMORY
open(DICT,$dict_file) or die "Can't open $dict_file!\n";
my @contig_order;
my $c=0;
while(<DICT>)
{
if($_=~ /\@SQ/)
	{
	my ($contig) = $_ =~ /SN:(\S+)/;
	$contig_order[$c]=$contig;
	++$c; 
	#print $contig,"\n";
	}
}
close(DICT);

#---------------------------------------- PARSE VCF FILE & OUTPUT SORTED VCF

open(VCF,$vcf_file) or die "Can't open $vcf_file!\n";

my %vcf_hash;
my $header;

while(<VCF>)
{
if($_=~/^#/){ $header .= $_; } # store header and comment fields
chomp($_);

my @data = split(/\t/,$_);
my $contig = $data[0];
my $start = $data[1];
my $variant = $data[4]."to".$data[5];
my $line = $_;

#print $contig,":",$start,"\n";

$vcf_hash{$contig}{$start}{$variant}=$line;

}
close(VCF);

#------------------ print out the VCF in the order of the reference genome

#print standard VCF header
print $header;

foreach my $contig (@contig_order) # sort by contig order
	{
	#print $contig,"\n";
	foreach my $start (sort {$a <=> $b} keys %{$vcf_hash{$contig}}) # sort numerically by coordinates
		{
		#print $start,"\n";
		foreach my $variant (keys %{$vcf_hash{$contig}{$start}}) # if overlapping mutation, print each variant
			{
			print $vcf_hash{$contig}{$start}{$variant},"\n";
			}	
		}
		
	}
