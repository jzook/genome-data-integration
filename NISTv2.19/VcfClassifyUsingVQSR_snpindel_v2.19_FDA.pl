#!/usr/bin/perl -w
#
# VcfClassifyUsingVQSR_snpindel_v2.13.pl - arbitrate uncertain genotype calls using VQSR tranches for SSEs, AB/QD, Alignment, and Mapping artifacts
# - similar to sqlupdateVcfUsingVQSR_v1.2.pl but adds indel arbitration, and it uses and outputs vcfs instead of a sqlite db
#
#
# Version 2.0.0 - Feb 8, 2013 - adapted from sqlupdateVcfUsingVQSR_v1.2.pl to adds indel arbitration, and use and output vcfs instead of a sqlite db
# Version 2.1 - Mar 4, 2013 - changed to use Hom VQSR (where HomRef and HomVar sites are combined, since HomRef VQSR often doesn't work
# Version 2.11 - Mar 12, 2013 - fixed bugs in code detecting end of files
# Version 2.12 - Apr 10, 2013 - add annotation indicating number of different platforms calling the genotype
# Version 2.13 - Apr 19, 2013 - fix parentheses for "for loop" in last tranche
# Version 2.14 - May 3, 2013 - output GQ field and correct QUAL field
# Version 2.15 - May 30, 2013 - change align VQSR arbitration to only filter if ReadPosEndDist exists as an annotation, since some real indels were being filtered just because the annotations don't exist for HaplotypeCaller calls.  ALso, If REF overlaps next variant and it isn't homref, then consider other vars uncertain
# Version 2.17 - Aug 21, 2013 - fix platforms and platformcount annotations; call as HomRef sites that are newly called in only 1 datasest when triggering haplotypecaller since most are 454 homopolymer errors; add annotation of platform-specific errors; fix VQSRSSE/AB/Align/Map getting out of sync, and correct arbitration to require SSE or Align rather than both
#   fix FILTER field and PLfiltered annotations so that header lines can conform to vcf spec
# Version 2.18 - Oct 29, 2013 - fix problem with CG call taking precedence when it has a different representation of a complex variant, as at 1:89923327
#	expand window around HC calls from 20bp to 35bp
#   add filter of indels that aren't called by HC in >1 datasets
# Version 2.19 - Feb 10, 2014 - fix HapCallVar to look for maximum among datasets instead of just taking first one since sometimes one dataset would have 0


use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use FileHandle;
use English;
use File::Basename;
use List::Util qw[min max];
#use Try::Tiny;
#use DBI;

my $MINDATASETS=2; #minimum number of datasets making a genotype call for it to be "highly confident"
my $MINYESNORATIO=5; #minimum number of datasets making a genotype call for it to be "highly confident"


my $MINPLSUMHOMSNP = 120; #call homozygous SNP genotype if difference between PL of most likely genotype and next most likely genotype for sum of all datasets > MINPLSUMHOMSNP
my $MINPLBYDPSUMHOMSNP = 1.6; #call homozygous SNP genotype if difference between PL of most likely genotype and next most likely genotype divided by DP for sum of all datasets > MINPLBYDPSUMHOMSNP
my $MINPLSUMHETSNP = 200; #call heterozygous SNP genotype if difference between PL of most likely genotype and next most likely genotype for sum of all datasets > MINPLSUMHETSNP
my $MINPLBYDPSUMHETSNP = 6.8; #call heterozygous SNP genotype if difference between PL of most likely genotype and next most likely genotype divided by DP for sum of all datasets > MINPLBYDPSUMHETSNP

my $MINPLSUMHOMINDEL = 80; #call homozygous INDEL genotype if difference between PL of most likely genotype and next most likely genotype for sum of all datasets > MINPLSUMHOMINDEL
my $MINPLBYDPSUMHOMINDEL = 0.8; #call homozygous INDEL genotype if difference between PL of most likely genotype and next most likely genotype divided by DP for sum of all datasets > MINPLBYDPSUMHOMINDEL
my $MINPLSUMHETINDEL = 100; #call heterozygous INDEL genotype if difference between PL of most likely genotype and next most likely genotype for sum of all datasets > MINPLSUMHETINDEL
my $MINPLBYDPSUMHETINDEL = 3.4; #call heterozygous INDEL genotype if difference between PL of most likely genotype and next most likely genotype divided by DP for sum of all datasets > MINPLBYDPSUMHETINDEL

my $MAPGOODMIN=3; #minimum datasets with mapping tranche < 99 to make a Het call

#define arrays with minimum values above corresponding to Hom and Het genotypes (Hom in PL positions 0,2,5,9,14)
my @MinPLSumSNP = ($MINPLSUMHOMSNP,$MINPLSUMHETSNP,$MINPLSUMHOMSNP,$MINPLSUMHETSNP,$MINPLSUMHETSNP,$MINPLSUMHOMSNP,$MINPLSUMHETSNP,$MINPLSUMHETSNP,$MINPLSUMHETSNP,$MINPLSUMHOMSNP,$MINPLSUMHETSNP,$MINPLSUMHETSNP,$MINPLSUMHETSNP,$MINPLSUMHETSNP,$MINPLSUMHOMSNP);
my @MinPLbyDPSumSNP = ($MINPLBYDPSUMHOMSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHOMSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHOMSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHOMSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHETSNP,$MINPLBYDPSUMHOMSNP);
my @MinPLSumINDEL = ($MINPLSUMHOMINDEL,$MINPLSUMHETINDEL,$MINPLSUMHOMINDEL,$MINPLSUMHETINDEL,$MINPLSUMHETINDEL,$MINPLSUMHOMINDEL,$MINPLSUMHETINDEL,$MINPLSUMHETINDEL,$MINPLSUMHETINDEL,$MINPLSUMHOMINDEL,$MINPLSUMHETINDEL,$MINPLSUMHETINDEL,$MINPLSUMHETINDEL,$MINPLSUMHETINDEL,$MINPLSUMHOMINDEL);
my @MinPLbyDPSumINDEL = ($MINPLBYDPSUMHOMINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHOMINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHOMINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHOMINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHETINDEL,$MINPLBYDPSUMHOMINDEL);

my @ALTPL=(3,6,10,15); #number of PLs for each number of ALT alleles

#F(j/k) = (k*(k+1)/2)+j -> formula for determining PL order from genotype
my @GTs=("0/0","0/1","1/1","0/2","1/2","2/2","0/3","1/3","2/3","3/3","0/4","1/4","2/4","3/4","4/4");
my $line;
my $line2;

#print OUTVARALL "$chromVar[0]\t$posVar[0]\t$idVar[0]\t$refVarLong\t$altVarLong\t$PL0diffsum\t$filtfield\t${filterReason}geno=$genosum;platforms=$platformcount;platformnames=$platformtxt;platformbias=$platbias;varType=$varType;YesPLtot=$YesPLTot[$genosum-1];NoPLTot=$NoPLTot[$genosum-1];HapNoVar=$HapNoVar;PLminsum=$PLminsum;DPSum=$DPSum;PLminsumOverDP=$PLbyDPsum;genoMapGood=$mapgood[$genosum-1];TrancheSSEmin2=$TrancheSSEmin2;TrancheABQDmin2=$TrancheABQDmin2;TrancheAlignmin2=$TrancheAlignmin2;TrancheMapmin2=$TrancheMapmin2;datasetcalls=$datasetcalls;allalts=$altsall$allPLannot$hruntxt$RUtxt$RPAtxt\tGT:DP:GQ:PL\t$GTs[$genosum-1]:$DPSum:$PLminsum:$PLTottxt\n";

#create header
my $header = "##fileformat=VCFv4.1\
##FILTER=<ID=Uncertain,Description=\"Uncertain genotype due to reason in filter INFO field\">\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth summed across all datasets, excluding MQ0 reads\">\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Net Genotype quality across all datasets, defined as difference between most likely and next most likely genotype likelihoods\">\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Net Genotype across all datasets\">\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods summed across all unfiltered datasets for genotypes as defined in the VCF specification\">\
##INFO=<ID=allalts,Number=1,Type=Integer,Description=\"All ALT alleles originally considered at this position\">\
##INFO=<ID=datasetcalls,Number=1,Type=Integer,Description=\"Number of datasets with any genotype call at this position\">\
##INFO=<ID=DPSum,Number=1,Type=Integer,Description=\"Total read depth summed across all datasets, excluding MQ0 reads\">\
##INFO=<ID=Entropy,Number=1,Type=Float,Description=\"Shannon entropy of variant flanking regions, 12bp on both sides\">\
##INFO=<ID=filter,Number=1,Type=String,Description=\"Reason for filtering this genotype as uncertain\">\
##INFO=<ID=geno,Number=1,Type=Integer,Description=\"Most probable genotype, corresponding to the minimum entry in the PL field (e.g., 1=0/0,2=0/1,3=1/1,4=0/2,etc)\">\
##INFO=<ID=genoMapGood,Number=1,Type=Integer,Description=\"Number of datasets calling this genotype with VQSR mapping tranche <= 95\">\
##INFO=<ID=HapNoVar,Number=1,Type=Integer,Description=\"Number of datasets for which HaplotypeCaller called a variant within 35bp and did not call a variant at this location\">\
##INFO=<ID=HRun,Number=1,Type=Integer,Description=\"Largest Contiguous Homopolymer Run of Variant Allele In Either Direction\">\
##INFO=<ID=NoCG,Number=0,Type=Flag,Description=\"Present if no consensus was reached for arbitration of all datasets, so we looked at all datasets except Complete Genomics since it may have a different representation of complex variants\">\
##INFO=<ID=NoPLTot,Number=1,Type=Integer,Description=\"Number of datasets with likelihood ratio > 20 for a genotype different from the called genotype\">\
##INFO=<ID=platforms,Number=1,Type=Integer,Description=\"Number of different platforms that called this genotype\">\
##INFO=<ID=platformbias,Number=.,Type=String,Description=\"Names of platforms that have at more than twice as many incorrect than correct genotypes at this location, indicating platform-specific bias (ill=Illumina,sol=SOLiD,454=454,ion=Ion Torrent,cg=Complete Genomics)\">\
##INFO=<ID=platformnames,Number=.,Type=String,Description=\"Names of platforms that called this genotype (ill=Illumina,sol=SOLiD,454=454,ion=Ion Torrent,cg=Complete Genomics)\">\
##INFO=<ID=PL454WG,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~16x 454 whole genome sequencing from 1000 Genomes Project, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLCG,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~73x Complete Genomics whole genome sequencing, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLHSWEx,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~66x 2x100bp Illumina exome sequencing from Broad Institute, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLHSWG,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~68x 2x100bp Illumina whole genome sequencing from Broad Institute, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLILL250,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~50x 2x250bp Illumina PCR-free whole genome sequencing from Broad Institute, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLILLCLIA,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~80x 2x100bp Illumina whole genome sequencing from Illumina CLIA lab, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLIllPCRFree,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~56x 2x100bp Illumina PCR-free whole genome sequencing from Illumina Platinum Genomes Project, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLILLWEx,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~30x 2x54bp Illumina exome sequencing from Broad Institute, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLILLWG,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~39x 2x44bp Illumina whole genome sequencing from Broad Institute, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLIonEx,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~80x mean 237bp Ion Torrent exome sequencing from Life Technologies, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLPlatGen,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~190x 2x100bp Illumina PCR-free whole genome sequencing from Illumina Platinum Genomes Project, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLXIll,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~37x 2x100bp Illumina whole genome sequencing from X Prize, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLXPSolWGLS,Number=.,Type=String,Description=\"Genotype likelihoods (PL) for ~24x 50bpx35bp SOLiD whole genome sequencing from X Prize, preceded by filtering info if this dataset was not used due to evidence of bias\">\
##INFO=<ID=PLminsum,Number=1,Type=Integer,Description=\"Net Genotype quality across all datasets, defined as difference between most likely and next most likely genotype likelihoods\">\
##INFO=<ID=PLminsumOverDP,Number=1,Type=Float,Description=\"Net Genotype quality across all datasets, defined as difference between most likely and next most likely genotype likelihoods, divided by the depth of coverage\">\
##INFO=<ID=RU,Number=1,Type=String,Description=\"Tandem repeat unit (bases)\">\
##INFO=<ID=RPA,Number=.,Type=Integer,Description=\"Number of times tandem repeat unit is repeated, for each allele (including reference)\">\
##INFO=<ID=TrancheABQDmin2,Number=1,Type=Float,Description=\"2nd lowest VQSR tranche for the called genotype for annotations associated with abnormal allele balance (AB and QD)\">\
##INFO=<ID=TrancheAlignmin2,Number=1,Type=Float,Description=\"2nd lowest VQSR tranche for the called genotype for annotations associated with local alignment errors (distance from the end of the read and clipping)\">\
##INFO=<ID=TrancheMapmin2,Number=1,Type=Float,Description=\"2nd lowest VQSR tranche for the called genotype for annotations associated with mapping errors (mapping quality and depth of coverage)\">\
##INFO=<ID=TrancheSSEmin2,Number=1,Type=Float,Description=\"2nd lowest VQSR tranche for the called genotype for annotations associated with systematic sequencing errors (strand bias and neighboring base quality)\">\
##INFO=<ID=varType,Number=1,Type=String,Description=\"Type of variant\">\
##INFO=<ID=YesPLtot,Number=1,Type=Integer,Description=\"Number of datasets with likelihood ratio > 20 for the called genotype\">\
##contig=<ID=1,length=249250621,assembly=b37>\
##contig=<ID=2,length=243199373,assembly=b37>\
##contig=<ID=3,length=198022430,assembly=b37>\
##contig=<ID=4,length=191154276,assembly=b37>\
##contig=<ID=5,length=180915260,assembly=b37>\
##contig=<ID=6,length=171115067,assembly=b37>\
##contig=<ID=7,length=159138663,assembly=b37>\
##contig=<ID=8,length=146364022,assembly=b37>\
##contig=<ID=9,length=141213431,assembly=b37>\
##contig=<ID=10,length=135534747,assembly=b37>\
##contig=<ID=11,length=135006516,assembly=b37>\
##contig=<ID=12,length=133851895,assembly=b37>\
##contig=<ID=13,length=115169878,assembly=b37>\
##contig=<ID=14,length=107349540,assembly=b37>\
##contig=<ID=15,length=102531392,assembly=b37>\
##contig=<ID=16,length=90354753,assembly=b37>\
##contig=<ID=17,length=81195210,assembly=b37>\
##contig=<ID=18,length=78077248,assembly=b37>\
##contig=<ID=19,length=59128983,assembly=b37>\
##contig=<ID=20,length=63025520,assembly=b37>\
##contig=<ID=21,length=48129895,assembly=b37>\
##contig=<ID=22,length=51304566,assembly=b37>\
##contig=<ID=X,length=155270560,assembly=b37>\
##contig=<ID=Y,length=59373566,assembly=b37>\
##contig=<ID=MT,length=16569,assembly=b37>\
##fileDate=20130719\
##phasing=none\
##reference=human_g1k_v37.fasta\
##variants_justified=left\
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878\n";


# Check if the right number of parameters are given at start up
#
# The following form is needed at start up:

if ($#ARGV < 5) {
    print "$#ARGV arguments\n usage: VcfClassifyUsingVQSR_snpindel_v2.18.pl inputvcfUGHapMergefilenamestart chromosomeName datasetNames(>=2 separated by spaces and alternating with platform identifier)\n";
    print "example: perl /Applications/bioinfo/perl/VcfClassifyUsingVQSR_snpindel_v2.18.pl AllFDAdatasets_130517_ all HSWG ill ILLWG ill XIll ill 454WG 454 ILLCLIA ill IllPCRFree ill XPSolWGLS sol IonEx ion HSWEx ill ILLWEx ill CG cg PlatGen ill ILL250 ill\n";
 exit;
}



my $vcfcnt = 2;

my $infilestart = "$ARGV[0]";
my $chrom = "$ARGV[1]";
my @infiles;
my @infilehandlesDistHomAB;
my @infilehandlesDistHomAlign;
my @infilehandlesDistHomSSE;
my @infilehandlesDistHomMap;
my @infilehandlesDistHetAB;
my @infilehandlesDistHetAlign;
my @infilehandlesDistHetSSE;
my @infilehandlesDistHetMap;
my @platforms=((0) x (($#ARGV-1/2)));
my @platnames=(("") x 10);
my $plats=0;

while ($vcfcnt <= $#ARGV) {
	$infiles[($vcfcnt-2)/2] = $ARGV[$vcfcnt];
	$platforms[($vcfcnt-2)/2] = $ARGV[$vcfcnt+1];
    
    my $j=0; #add platform to platnames if it's new
    for (my $i=0; $i<$plats; $i++) {
        if ($platforms[($vcfcnt-2)/2] eq $platnames[$i]) { $j=1; }
    }
    if ($j==0) {
        $platnames[$plats] = $platforms[($vcfcnt-2)/2];
        $plats++;
    }
    
	local *FILEIN;
	open (FILEIN, "${infilestart}${infiles[($vcfcnt-2)/2]}call_UGHapMerge_Hom_AB_recal90_${chrom}.vcf") || die;
	push(@infilehandlesDistHomAB, *FILEIN);
	local *FILEIN2;
	open (FILEIN2, "${infilestart}${infiles[($vcfcnt-2)/2]}call_UGHapMerge_Hom_Align_recal90_${chrom}.vcf") || die;
	push(@infilehandlesDistHomAlign, *FILEIN2);
	local *FILEIN3;
	open (FILEIN3, "${infilestart}${infiles[($vcfcnt-2)/2]}call_UGHapMerge_Hom_SSE_recal90_${chrom}.vcf") || die;
	push(@infilehandlesDistHomSSE, *FILEIN3);
	local *FILEIN4;
	open (FILEIN4, "${infilestart}${infiles[($vcfcnt-2)/2]}call_UGHapMerge_Hom_map_recal90_${chrom}.vcf") || die;
	push(@infilehandlesDistHomMap, *FILEIN4);
	local *FILEIN5;
	open (FILEIN5, "${infilestart}${infiles[($vcfcnt-2)/2]}call_UGHapMerge_Het_AB_recal90_${chrom}.vcf") || die;
	push(@infilehandlesDistHetAB, *FILEIN5);
	local *FILEIN6;
	open (FILEIN6, "${infilestart}${infiles[($vcfcnt-2)/2]}call_UGHapMerge_Het_Align_recal90_${chrom}.vcf") || die;
	push(@infilehandlesDistHetAlign, *FILEIN6);
	local *FILEIN7;
	open (FILEIN7, "${infilestart}${infiles[($vcfcnt-2)/2]}call_UGHapMerge_Het_SSE_recal90_${chrom}.vcf") || die;
	push(@infilehandlesDistHetSSE, *FILEIN7);
	local *FILEIN8;
	open (FILEIN8, "${infilestart}${infiles[($vcfcnt-2)/2]}call_UGHapMerge_Het_map_recal90_${chrom}.vcf") || die;
	push(@infilehandlesDistHetMap, *FILEIN8);
	$vcfcnt+=2;
}
$vcfcnt -= 2;


# Output file is opened
unless ( open(OUTPUT, ">${infilestart}allcall_UGHapMerge_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTHOMREF, ">${infilestart}allcall_UGHapMerge_HomRef_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_HomRef_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTHET, ">${infilestart}allcall_UGHapMerge_Het_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_Het_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTHOMVAR, ">${infilestart}allcall_UGHapMerge_HomVar_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_HomVar_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTVARPASS, ">${infilestart}allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTVARALL, ">${infilestart}allcall_UGHapMerge_HetHomVarAll_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_HetHomVarAll_VQSRv2.19_${MINDATASETS}mindatasets_${MINYESNORATIO}minYesNoRatio_${chrom}.vcf to write to! \n\n";
    exit;
}
#skip header lines in individual files and write first one to output
my $fhn = 0;
foreach my $fh (@infilehandlesDistHomAB) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
        if ($x==1) {last;}
    }
    $fhn++;
    
}
      print OUTPUT $header;
		print OUTHOMREF $header;
		print OUTHET $header;
		print OUTHOMVAR $header;
		print OUTVARPASS $header;
		print OUTVARALL $header;


foreach my $fh (@infilehandlesDistHomAlign) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
        if ($x==1) {last;}
    }
    $fhn++;    
}
foreach my $fh (@infilehandlesDistHomSSE) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
            if ($x==1) {last;}
    }
    $fhn++;
}
foreach my $fh (@infilehandlesDistHomMap) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
            if ($x==1) {last;}
    }
    $fhn++;
}

foreach my $fh (@infilehandlesDistHetAB) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
        if ($x==1) {last;}
    }
    $fhn++;
}

foreach my $fh (@infilehandlesDistHetAlign) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
            if ($x==1) {last;}
    }
    $fhn++;
}
foreach my $fh (@infilehandlesDistHetSSE) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
            if ($x==1) {last;}
    }
    $fhn++;
}
foreach my $fh (@infilehandlesDistHetMap) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
            if ($x==1) {last;}
    }
    $fhn++;
}


my @fields;



my @VQSRHomAB=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomSSE=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomAlign=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomMap=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHetAB=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHetSSE=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHetAlign=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHetMap=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomVarAB=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomVarSSE=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomVarAlign=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomVarMap=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomABprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomSSEprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomAlignprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomMapprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHetABprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHetSSEprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHetAlignprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHetMapprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomVarABprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomVarSSEprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomVarAlignprev=((0) x ($#infilehandlesDistHomAB+1));
my @VQSRHomVarMapprev=((0) x ($#infilehandlesDistHomAB+1));

my @chromposHom=((0) x ($#infilehandlesDistHomAB+1));
my @chromposHet=((0) x ($#infilehandlesDistHomAB+1));
my @chromposHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @chromposHomprev=((0) x ($#infilehandlesDistHomAB+1));
my @chromposHetprev=((0) x ($#infilehandlesDistHomAB+1));
my @chromposHomVarprev=((0) x ($#infilehandlesDistHomAB+1));
my @chromHom=((0) x ($#infilehandlesDistHomAB+1));
my @chromHet=((0) x ($#infilehandlesDistHomAB+1));
my @chromHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @posHom=((0) x ($#infilehandlesDistHomAB+1));
my @posHet=((0) x ($#infilehandlesDistHomAB+1));
my @posHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @idHom=((0) x ($#infilehandlesDistHomAB+1));
my @idHet=((0) x ($#infilehandlesDistHomAB+1));
my @idHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @refHom=((0) x ($#infilehandlesDistHomAB+1));
my @refHet=((0) x ($#infilehandlesDistHomAB+1));
my @refHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @altHom=((0) x ($#infilehandlesDistHomAB+1));
my @altHet=((0) x ($#infilehandlesDistHomAB+1));
my @altHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @qualHom=((0) x ($#infilehandlesDistHomAB+1));
my @qualHet=((0) x ($#infilehandlesDistHomAB+1));
my @qualHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @infoHom=((0) x ($#infilehandlesDistHomAB+1));
my @infoHet=((0) x ($#infilehandlesDistHomAB+1));
my @infoHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @formatHom=((0) x ($#infilehandlesDistHomAB+1));
my @formatHet=((0) x ($#infilehandlesDistHomAB+1));
my @formatHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @charHom=((0) x ($#infilehandlesDistHomAB+1));
my @charHet=((0) x ($#infilehandlesDistHomAB+1));
my @charHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @PLtextHom=(("") x ($#infilehandlesDistHomAB+1));
my @PLtextHet=(("") x ($#infilehandlesDistHomAB+1));
my @PLtextHomVar=(("") x ($#infilehandlesDistHomAB+1));
my @genoHom=((0) x ($#infilehandlesDistHomAB+1));
my @genoHet=((0) x ($#infilehandlesDistHomAB+1));
my @genoHomVar=((0) x ($#infilehandlesDistHomAB+1));
my @genoHomprev=((0) x ($#infilehandlesDistHomAB+1));
my @genoHetprev=((0) x ($#infilehandlesDistHomAB+1));
my @genoHomVarprev=((0) x ($#infilehandlesDistHomAB+1));
my @infoHomprev=((0) x ($#infilehandlesDistHomAB+1));
my @infoHetprev=((0) x ($#infilehandlesDistHomAB+1));
my @infoHomVarprev=((0) x ($#infilehandlesDistHomAB+1));
my @PLtextHomprev=(("") x ($#infilehandlesDistHomAB+1));
my @PLtextHetprev=(("") x ($#infilehandlesDistHomAB+1));
my @PLtextHomVarprev=(("") x ($#infilehandlesDistHomAB+1));
my @qualHetprev=((0) x ($#infilehandlesDistHomAB+1));
my @qualHomVarprev=((0) x ($#infilehandlesDistHomAB+1));
my $chromposdelregion=0; #specify region each confident variant REF field covers so that we consider other variants inside it to be uncertain


my $minchrompos=0; #previous minimum chrompos
my $newline=1; #printed out a new line, so ready for next one

my @endlineHom = ((1) x ($#infilehandlesDistHomAB+1));
my @endlineHet = ((1) x ($#infilehandlesDistHomAB+1));
my @endlineHomVar = ((1) x ($#infilehandlesDistHomAB+1));
my $endlineall=1;
while($endlineall==1) {

    my $PLno=0;
    my @HapNoVarHom=((0) x ($#infilehandlesDistHomAB+1));
    my @HapNoVarHet=((0) x ($#infilehandlesDistHomAB+1));
    my @HapNoVarHomVar=((0) x ($#infilehandlesDistHomAB+1));
    
    ##########
    #Load Hom vcfs
    ##########
    my $i=0;
    my @readnewline=((0) x ($#infilehandlesDistHomAB+1));
    foreach my $fh (@infilehandlesDistHomAB) {
		if ($endlineHom[$i]==1 && $chromposHom[$i]==$minchrompos) {
            #print "\n$i $endline[$i] $chrompos[$i] $minchrompos\n";
            #readline($i);
            my $dataset = $i;
            if ($genoHom[$i]>1) {
                $chromposHomprev[$i]=$chromposHom[$i];
                $infoHomprev[$i] = $infoHom[$i];
                $PLtextHomprev[$i] = $PLtextHom[$i];
                $genoHomprev[$i] = $genoHom[$i];
                $VQSRHomABprev[$i]=$VQSRHomAB[$i];
                $VQSRHomAlignprev[$i]=$VQSRHomAlign[$i];
                $VQSRHomSSEprev[$i]=$VQSRHomSSE[$i];
                $VQSRHomMapprev[$i]=$VQSRHomMap[$i];
            }
            
            my $hapref=1; #continue reading next lines while
            while ($hapref==1) {
                my $line = <$fh>;
                $readnewline[$i]+=1; #read new lines in SSE, Align, and map
                #print $line[$i];
                if (defined($line)) {
                } else { #quit loop if reach end of file
                    $endlineHom[$i] = 0;
                    $chromposHom[$i]=1000000000000; 
                    last;
                }
                
                # Split up the line into an array
                @fields = split( "\t", $line);
                my $chrno=$fields[0];
                if ($chrno eq "M" || $chrno eq "MT") {
                    $chrno=23;
                } elsif ($chrno eq "X") {$chrno=24;
                } elsif ($chrno eq "Y") {$chrno=25;
                } elsif ($chrno =~ /^\d*/) {
                } else { $endlineHom[$i] = 0; $chromposHom[$i]=1000000000000; last;}
                $chromposHom[$i] = $chrno*1000000000+$fields[1];
                #print "chrpos=$chrompos[$i]\n";
                
                
                $chromHom[$i] = $fields[0];
                $posHom[$i] = $fields[1];
                $idHom[$i] = $fields[2];
                $refHom[$i] = $fields[3];
                $altHom[$i] = $fields[4];
                $qualHom[$i] = $fields[5];
                $infoHom[$i] = $fields[7];
                $formatHom[$i] = $fields[8];
                $charHom[$i] = $fields[9];
                
                if ($infoHom[$i] =~ /UGHap=HapRefUGVar/) {
                    $HapNoVarHom[$i]=1;
                    next;
                } else { $HapNoVarHom[$i]=0; $hapref=0; }
                
                
                my $filter=$fields[6];
                if ($filter =~ /99.90to/) { $VQSRHomAB[$i]=99.9;
                } elsif ($filter =~ /99.50to/) { $VQSRHomAB[$i]=99.5;
                } elsif ($filter =~ /99.00to/) { $VQSRHomAB[$i]=99;
                } elsif ($filter =~ /95.00to/) { $VQSRHomAB[$i]=95;
                } elsif ($filter =~ /90.00to/) { $VQSRHomAB[$i]=90;
                } else { $VQSRHomAB[$i]=0; }
                
                my @formats = split( ":", $formatHom[$i]);
                my @chars = split( ":", $charHom[$i]);
                $PLtextHom[$i]="";
                for (my $j=0; $j<=$#formats; $j++) {
                    #if ($formats[$i] eq "GT") { $gt[$i]=$chars[$j]; }
                    if ($formats[$j] eq "PL") {
                        if ($chars[$j] =~ /(.*)\n/) {$PLtextHom[$i]=$1;}
                        else { $PLtextHom[$i]=$chars[$j];}
                    }
                }
                if (!($PLtextHom[$i] eq "")) {
                    my @PLs = split( ",", $PLtextHom[$i]);
                    $PLno=$#PLs;
                    my $j=0; my $PL0=0; my $PLmin=1000000;
                    foreach my $PL (@PLs) {
                        if ($PL==0) { $genoHom[$i]=$j+1; }
                        $j++;
                    }
                }
                
            }
		}
        $i++;
    }

    $i=0;
    foreach my $fh (@infilehandlesDistHomSSE) {
		for (my $j=0;$j<$readnewline[$i];$j++) {
            #SSE
            my $line = <$fh>;
            @fields = split( "\t", $line);
            my $filter=$fields[6];
            if ($filter =~ /99.90to/) { $VQSRHomSSE[$i]=99.9;
            } elsif ($filter =~ /99.50to/) { $VQSRHomSSE[$i]=99.5;
            } elsif ($filter =~ /99.00to/) { $VQSRHomSSE[$i]=99;
            } elsif ($filter =~ /95.00to/) { $VQSRHomSSE[$i]=95;
            } elsif ($filter =~ /90.00to/) { $VQSRHomSSE[$i]=90;
            } else { $VQSRHomSSE[$i]=0; }
        }
        $i++;
    }
    $i=0;
    foreach my $fh (@infilehandlesDistHomAlign) {
		for (my $j=0;$j<$readnewline[$i];$j++) {
            #Align
            my $line = <$fh>;
            @fields = split( "\t", $line);
            my $filter=$fields[6];
            $VQSRHomAlign[$i]=0;
            if ($fields[7] =~ /ReadPosEndDist/) { #only use VQSRAlign if ReadPosEndDist has a value
                if ($filter =~ /99.90to/) { $VQSRHomAlign[$i]=99.9;
                } elsif ($filter =~ /99.50to/) { $VQSRHomAlign[$i]=99.5;
                } elsif ($filter =~ /99.00to/) { $VQSRHomAlign[$i]=99;
                } elsif ($filter =~ /95.00to/) { $VQSRHomAlign[$i]=95;
                } elsif ($filter =~ /90.00to/) { $VQSRHomAlign[$i]=90;
                }
            }
        }
        $i++;
    }
    $i=0;
    foreach my $fh (@infilehandlesDistHomMap) {
		for (my $j=0;$j<$readnewline[$i];$j++) {
            #Map
            my $line = <$fh>;
            @fields = split( "\t", $line);
            my $filter=$fields[6];
            if ($filter =~ /99.90to/) { $VQSRHomMap[$i]=99.9;
            } elsif ($filter =~ /99.50to/) { $VQSRHomMap[$i]=99.5;
            } elsif ($filter =~ /99.00to/) { $VQSRHomMap[$i]=99;
            } elsif ($filter =~ /95.00to/) { $VQSRHomMap[$i]=95;
            } elsif ($filter =~ /90.00to/) { $VQSRHomMap[$i]=90;
            } else { $VQSRHomMap[$i]=0; }
        }
        $i++;
    }

    
    ##########
    #Load Het vcfs
    ##########
    $i=0;
    @readnewline=((0) x ($#infilehandlesDistHetAB+1));
    foreach my $fh (@infilehandlesDistHetAB) {
		if ($endlineHet[$i]==1 && $chromposHet[$i]==$minchrompos) {
            #print "\n$i $endline[$i] $chrompos[$i] $minchrompos\n";
            #readline($i);
            my $dataset = $i;
            if ($genoHet[$i]>1) {
                $chromposHetprev[$i]=$chromposHet[$i];
                $infoHetprev[$i] = $infoHet[$i];
                $PLtextHetprev[$i] = $PLtextHet[$i];
                $genoHetprev[$i] = $genoHet[$i];
                $VQSRHetABprev[$i]=$VQSRHetAB[$i];
                $VQSRHetAlignprev[$i]=$VQSRHetAlign[$i];
                $VQSRHetSSEprev[$i]=$VQSRHetSSE[$i];
                $VQSRHetMapprev[$i]=$VQSRHetMap[$i];
            }
            
            my $hapref=1; #continue reading next lines while
            while ($hapref==1) {
                my $line = <$fh>;
                $readnewline[$i]+=1; #read new lines in SSE, Align, and map
                #print $line[$i];
                if (defined($line)) {
                } else { #quit loop if reach end of file
                    $endlineHet[$i] = 0;
                    $chromposHet[$i]=1000000000000; 
                    last;
                }
                
                # Split up the line into an array
                @fields = split( "\t", $line);
                my $chrno=$fields[0];
                if ($chrno eq "M" || $chrno eq "MT") {
                    $chrno=23;
                } elsif ($chrno eq "X") {$chrno=24;
                } elsif ($chrno eq "Y") {$chrno=25;
                } elsif ($chrno =~ /^\d*/) {
                } else { $endlineHet[$i] = 0; $chromposHet[$i]=1000000000000; last;}
                $chromposHet[$i] = $chrno*1000000000+$fields[1];
                #print "chrpos=$chrompos[$i]\n";
                
                
                $chromHet[$i] = $fields[0];
                $posHet[$i] = $fields[1];
                $idHet[$i] = $fields[2];
                $refHet[$i] = $fields[3];
                $altHet[$i] = $fields[4];
                $qualHet[$i] = $fields[5];
                $infoHet[$i] = $fields[7];
                $formatHet[$i] = $fields[8];
                $charHet[$i] = $fields[9];
                
                if ($infoHet[$i] =~ /UGHap=HapRefUGVar/) {
                    $HapNoVarHet[$i]=1;
                    next;
                } else { $HapNoVarHet[$i]=0; $hapref=0; }
                
                
                my $filter=$fields[6];
                if ($filter =~ /99.90to/) { $VQSRHetAB[$i]=99.9;
                } elsif ($filter =~ /99.50to/) { $VQSRHetAB[$i]=99.5;
                } elsif ($filter =~ /99.00to/) { $VQSRHetAB[$i]=99;
                } elsif ($filter =~ /95.00to/) { $VQSRHetAB[$i]=95;
                } elsif ($filter =~ /90.00to/) { $VQSRHetAB[$i]=90;
                } else { $VQSRHetAB[$i]=0; }
                
                my @formats = split( ":", $formatHet[$i]);
                my @chars = split( ":", $charHet[$i]);
                $PLtextHet[$i]="";
                for (my $j=0; $j<=$#formats; $j++) {
                    #if ($formats[$i] eq "GT") { $gt[$i]=$chars[$j]; }
                    if ($formats[$j] eq "PL") {
                        if ($chars[$j] =~ /(.*)\n/) {$PLtextHet[$i]=$1;}
                        else { $PLtextHet[$i]=$chars[$j];}
                    }
                }
                if (!($PLtextHet[$i] eq "")) {
                    my @PLs = split( ",", $PLtextHet[$i]);
                    $PLno=$#PLs;
                    my $j=0; my $PL0=0; my $PLmin=1000000;
                    foreach my $PL (@PLs) {
                        if ($PL==0) { $genoHet[$i]=$j+1; }
                        $j++;
                    }
                }
                
            }
		}
        $i++;
    }
    
    $i=0;
    foreach my $fh (@infilehandlesDistHetSSE) {
		for (my $j=0;$j<$readnewline[$i];$j++) {
            #SSE
            my $line = <$fh>;
            @fields = split( "\t", $line);
            my $filter=$fields[6];
            if ($filter =~ /99.90to/) { $VQSRHetSSE[$i]=99.9;
            } elsif ($filter =~ /99.50to/) { $VQSRHetSSE[$i]=99.5;
            } elsif ($filter =~ /99.00to/) { $VQSRHetSSE[$i]=99;
            } elsif ($filter =~ /95.00to/) { $VQSRHetSSE[$i]=95;
            } elsif ($filter =~ /90.00to/) { $VQSRHetSSE[$i]=90;
            } else { $VQSRHetSSE[$i]=0; }
        }
        $i++;
    }
    $i=0;
    foreach my $fh (@infilehandlesDistHetAlign) {
		for (my $j=0;$j<$readnewline[$i];$j++) {
            #Align
            my $line = <$fh>;
            @fields = split( "\t", $line);
            my $filter=$fields[6];
            $VQSRHetAlign[$i]=0;
            if ($fields[7] =~ /ReadPosEndDist/) { #only use VQSRAlign if ReadPosEndDist has a value
                if ($filter =~ /99.90to/) { $VQSRHetAlign[$i]=99.9;
                } elsif ($filter =~ /99.50to/) { $VQSRHetAlign[$i]=99.5;
                } elsif ($filter =~ /99.00to/) { $VQSRHetAlign[$i]=99;
                } elsif ($filter =~ /95.00to/) { $VQSRHetAlign[$i]=95;
                } elsif ($filter =~ /90.00to/) { $VQSRHetAlign[$i]=90;
                }
            }
        }
        $i++;
    }
    $i=0;
    foreach my $fh (@infilehandlesDistHetMap) {
		for (my $j=0;$j<$readnewline[$i];$j++) {
            #Map
            my $line = <$fh>;
            @fields = split( "\t", $line);
            my $filter=$fields[6];
            if ($filter =~ /99.90to/) { $VQSRHetMap[$i]=99.9;
            } elsif ($filter =~ /99.50to/) { $VQSRHetMap[$i]=99.5;
            } elsif ($filter =~ /99.00to/) { $VQSRHetMap[$i]=99;
            } elsif ($filter =~ /95.00to/) { $VQSRHetMap[$i]=95;
            } elsif ($filter =~ /90.00to/) { $VQSRHetMap[$i]=90;
            } else { $VQSRHetMap[$i]=0; }
        }
        $i++;
    }

    
    
    
    if (max(@endlineHom)+max(@endlineHet)==0) { last; }
    my $minchromposprev=$minchrompos;
    $minchrompos = min(min(@chromposHom),min(@chromposHet));
    if ($minchrompos==1000000000000) { last; }
    if (int($minchrompos/1000000000)>int($minchromposprev/1000000000)) {
    	$chromposdelregion=0; #reset deletion region if new chromosome
    }
    if ($minchrompos==1000948929) {
    	print "before vqsr\n";
    }
    
    my @VQSRABHapNoVar=((0) x ($#infilehandlesDistHomAB+1));
    my @VQSRSSEHapNoVar=((0) x ($#infilehandlesDistHomAB+1));
    my @VQSRAlignHapNoVar=((0) x ($#infilehandlesDistHomAB+1));
    my @VQSRMapHapNoVar=((0) x ($#infilehandlesDistHomAB+1));
    my @VQSRABVar=((0) x ($#infilehandlesDistHomAB+1));
    my @VQSRSSEVar=((0) x ($#infilehandlesDistHomAB+1));
    my @VQSRAlignVar=((0) x ($#infilehandlesDistHomAB+1));
    my @VQSRMapVar=((0) x ($#infilehandlesDistHomAB+1));
    
    
    my @genoHapNoVar=((0) x ($#infilehandlesDistHomAB+1));
    
    my @chromVar=((0) x ($#infilehandlesDistHomAB+1));
    my @posVar=((0) x ($#infilehandlesDistHomAB+1));
    my @idVar=((0) x ($#infilehandlesDistHomAB+1));
    my @refVar=((0) x ($#infilehandlesDistHomAB+1));
    my @altVar=((0) x ($#infilehandlesDistHomAB+1));
    my @DPVar=((0) x ($#infilehandlesDistHomAB+1));
    my @PLtextVar=(("") x ($#infilehandlesDistHomAB+1));
    my @genoVar=((0) x ($#infilehandlesDistHomAB+1));
    my @datasetVar=((0) x ($#infilehandlesDistHomAB+1));
    my @platformsVar=(("") x ($#infilehandlesDistHomAB+1));
    
    my $HapCallVar=-1;
    my $hrun=-1;
    my $RU="";
    my $RPA="";
    if ($minchrompos==1000948929) {
    	print "before PLs\n";
    }

    ###########
    # Find the PLs and genotype for each dataset that has a row at this position
    ###########
    my $datasetcalls=0; my $HapNoVar=0;
    for (my $i=0;$i<=$#infilehandlesDistHomAB;$i++) {
		if ($chromposHom[$i]!=$minchrompos && $chromposHet[$i]!=$minchrompos) {
            if ($chromposHet[$i]<=$minchrompos+35) {
                if (!($qualHet[$i] eq ".") && $qualHet[$i]>40 && ($infoHet[$i] =~ /UGHap=Hap/ || $infoHet[$i] =~ /UGHap=Both/ || $infoHet[$i] =~ /UGHap=both/)) {
                    if (!($PLtextHet[$i] eq "")) {
                        my @PLs = split( ",", $PLtextHet[$i]);
                        $PLno=$#PLs;
                        my $j=0; my $PL0=0; my $PLmin=1000000;
                        foreach my $PL (@PLs) {
                            if ($PL==0) { $genoHapNoVar[$HapNoVar]=$j; }
                            $j++;
                        }
                    }
                    $VQSRABHapNoVar[$HapNoVar]=$VQSRHetAB[$i];
                    $VQSRSSEHapNoVar[$HapNoVar]=$VQSRHetSSE[$i];
                    $VQSRAlignHapNoVar[$HapNoVar]=$VQSRHetAlign[$i];
                    $VQSRMapHapNoVar[$HapNoVar]=$VQSRHetMap[$i];
                    $HapNoVar++;
                }
            } elsif ($chromposHetprev[$i]>=$minchrompos-35) {
                if (!($qualHetprev[$i] eq ".") && $qualHetprev[$i]>40 && ($infoHetprev[$i] =~ /UGHap=Hap/ || $infoHetprev[$i] =~ /UGHap=Both/ || $infoHetprev[$i] =~ /UGHap=both/)) {
                    $genoHapNoVar[$HapNoVar]=$genoHetprev[$i];
                    $VQSRABHapNoVar[$HapNoVar]=$VQSRHetABprev[$i];
                    $VQSRSSEHapNoVar[$HapNoVar]=$VQSRHetSSEprev[$i];
                    $VQSRAlignHapNoVar[$HapNoVar]=$VQSRHetAlignprev[$i];
                    $VQSRMapHapNoVar[$HapNoVar]=$VQSRHetMapprev[$i];
                    $HapNoVar++;
                }
            }
		} else { #dataset has a genotype call at this position
        
            if ($infoHom[$i] =~ /HapCallVar=(.*?);/) {
                if ($1>$HapCallVar) {$HapCallVar=$1;}
            }
            if ($infoHet[$i] =~ /HapCallVar=(.*?);/) {
                if ($1>$HapCallVar) {$HapCallVar=$1;}
            }

            if ($hrun==-1 && $infoHom[$i] =~ /HRun=(.*?);/) {
                $hrun=$1;
            }
            if ($hrun==-1 && $infoHet[$i] =~ /HRun=(.*?);/) {
                $hrun=$1;
            }

            if ($RU eq "" && $infoHom[$i] =~ /RU=(.*?);/) {
                $RU=$1;
            }
            if ($RU eq "" && $infoHet[$i] =~ /RU=(.*?);/) {
                $RU=$1;
            }
            
            if ($RPA eq "" && $infoHom[$i] =~ /RPA=(.*?);/) {
                $RPA=$1;
            }
            if ($RPA eq "" && $infoHet[$i] =~ /RPA=(.*?);/) {
                $RPA=$1;
            }

            if ($chromposHom[$i]==$minchrompos && ($genoHom[$i]==1 || $genoHom[$i]==3 || $genoHom[$i]==6 || $genoHom[$i]==10 || $genoHom[$i]==15)) {
                
                #if this dataset has a haplotypecaller variant call within 35 bases but no variant at this position, then something went wrong
                if ($infoHom[$i] =~ /UGHap=HapRefUGVar/) {
                    print "HapRefUGVar got through";
                    exit;
                }
                
                if (!($PLtextHom[$i] eq "")) {
                    $datasetVar[$datasetcalls]=$i;
                    $chromVar[$datasetcalls] = $chromHom[$i];
                    $posVar[$datasetcalls] = $posHom[$i];
                    $idVar[$datasetcalls] = $idHom[$i];
                    $refVar[$datasetcalls] = $refHom[$i];
                    $altVar[$datasetcalls] = $altHom[$i];
                    $PLtextVar[$datasetcalls]=$PLtextHom[$i];
                    if ($infoHom[$i] =~ /DP=(.*?);/) { $DPVar[$datasetcalls] = $1; }
                    $genoVar[$datasetcalls]=$genoHom[$i];
                    $VQSRABVar[$datasetcalls]=$VQSRHomAB[$i];
                    $VQSRSSEVar[$datasetcalls]=$VQSRHomSSE[$i];
                    $VQSRAlignVar[$datasetcalls]=$VQSRHomAlign[$i];
                    $VQSRMapVar[$datasetcalls]=$VQSRHomMap[$i];
                    $platformsVar[$datasetcalls]=$platforms[$i];
                    $datasetcalls++;
                }

            } elsif ($chromposHet[$i]==$minchrompos && ($genoHet[$i]!=1 && $genoHet[$i]!=3 && $genoHet[$i]!=6 && $genoHet[$i]!=10 && $genoHet[$i]!=15)) {
                
                #if this dataset has a haplotypecaller variant call within 35 bases but no variant at this position, then something went wrong
                if ($infoHet[$i] =~ /UGHap=HapRefUGVar/) {
                    print "HapRefUGVar got through";
                    exit;
                }
                
                if (!($PLtextHet[$i] eq "")) {
                    $datasetVar[$datasetcalls]=$i;
                    $PLtextVar[$datasetcalls]=$PLtextHet[$i];
                    if ($infoHet[$i] =~ /DP=(.*?);/) { $DPVar[$datasetcalls] = $1; }
                    $chromVar[$datasetcalls] = $chromHet[$i];
                    $posVar[$datasetcalls] = $posHet[$i];
                    $idVar[$datasetcalls] = $idHet[$i];
                    $refVar[$datasetcalls] = $refHet[$i];
                    $altVar[$datasetcalls] = $altHet[$i];
                    $genoVar[$datasetcalls]=$genoHet[$i];
                    $VQSRABVar[$datasetcalls]=$VQSRHetAB[$i];
                    $VQSRSSEVar[$datasetcalls]=$VQSRHetSSE[$i];
                    $VQSRAlignVar[$datasetcalls]=$VQSRHetAlign[$i];
                    $VQSRMapVar[$datasetcalls]=$VQSRHetMap[$i];
                    $platformsVar[$datasetcalls]=$platforms[$i];
                    $datasetcalls++;
                }
                
            }
        }
        
    }
    
    if ($datasetcalls<2) { next; } #skip this record if only in one dataset, since it had to be missed by UG and HaplotypeCaller initially in all datasets and only called by HaplotypeCaller in one dataset when triggered
    
    ###TODO: take into account multiple Ref/Alt combinations
    my @altallVar=(("") x 6);
    my $varType="SNP";
    my $altno=0; my $altVarLong=""; my $refVarLong=""; my $multiplealts=0;
    for (my $kk=0; $kk<$datasetcalls; $kk++) {
        
        if (length($refVar[$kk])>1) {$varType="INDEL";}
        my @alts = split( ",", $altVar[$kk]);
        my $alt4="";
        my $j=0; 
        foreach my $alt (@alts) {
            if (length($alt)>1) {$varType="INDEL";}
            if ($j==0) {
            	$alt4=$alt;
            } elsif ($j<=$#ALTPL) { #if there are more than 4 alts, only include the first 4
            	$alt4="$alt4,$alt"; 
            } else {
            	$j--;
            }
            $j++;
        }
        if ($kk>0 && !($alt4 eq $altVarLong)) { $multiplealts=1;}
        if ($j>$altno) {$altno=$j; $altVarLong=$alt4; $refVarLong=$refVar[$kk];}
    }
    
    if ($minchrompos==1000948929) {
    	print "before out\n";
    }
    
	my $DPSum=0;
	
    my @PLTot=((0) x 15);    #allows up to 4 ALT alleles
    my @NoPLTot=((0) x 15);
    my @YesPLTot=((0) x 15);
    my @filteredgenos=((0) x 15);
    my @mapgood=((0) x 15);
    my @platformgenotxt=((",") x 15);
    my $genosum=0;
    
    my $allPLannot="";
    my $filtergeno="filtered";


    my $TrancheSSEmin1=100.0; my $TrancheSSEmin2=100.0; my $TrancheABQDmin1=100.0; my $TrancheABQDmin2=100.0; my $TrancheAlignmin1=100.0; my $TrancheAlignmin2=100.0; my $TrancheMapmin1=100.0; my $TrancheMapmin2=100.0;
    
    #first, check if there is a consensus genotype without filtering
    for (my $kk=0; $kk<$datasetcalls; $kk++) {
        #print "$chromVar[$kk],$posVar[$kk],$idVar[$kk],$refVar[$kk],$altVar[$kk],$genoVar[$kk],$VQSRABVar[$kk],$VQSRSSEVar[$kk],$VQSRAlignVar[$kk],$VQSRMapVar[$kk],$datasetcalls,$PLtextVar[$kk]\n";
        #exit;
        $allPLannot="$allPLannot;PL${infiles[$datasetVar[$kk]]}=$PLtextVar[$kk]";
        if ($PLtextVar[$kk] eq "") { next; } #no coverage
        my @PLs = split( ",", $PLtextVar[$kk]);
        my $j=0; my $PL0=0; my $PLmin=1000000;
        foreach my $PL (@PLs) {
            $PLTot[$j] += $PL;
            if ($PL>20) { #confidently not this genotype in dataset
                $NoPLTot[$j]+=1;
            } 
            if ($PL>0 && $PL<$PLmin) { $PLmin=$PL; } #find the minimum PL>0, since this is the likelihood ratio of the most likely to the next most likely genotype
            $j++;
        }
        if ($PLmin>20) {
            #add this platform if it is a new platform calling this genotype
            if ($platformgenotxt[$genoVar[$kk]-1] =~ /,$platformsVar[$kk],/) {} else {
                $platformgenotxt[$genoVar[$kk]-1] = "$platformgenotxt[$genoVar[$kk]-1]$platformsVar[$kk],";
            }
            
            $YesPLTot[$genoVar[$kk]-1]++;
            if ($VQSRMapVar[$kk]<=95) {$mapgood[$genoVar[$kk]-1]++;}
            if ($VQSRSSEVar[$kk]<=$TrancheSSEmin1 && $VQSRSSEVar[$kk]>=0) {
                $TrancheSSEmin2=$TrancheSSEmin1;
                $TrancheSSEmin1=$VQSRSSEVar[$kk];
            } elsif ($VQSRSSEVar[$kk]<$TrancheSSEmin2 && $VQSRSSEVar[$kk]>=0) {
                $TrancheSSEmin2=$VQSRSSEVar[$kk];
            }
            if ($VQSRABVar[$kk]<=$TrancheABQDmin1 && $VQSRABVar[$kk]>=0) {
                $TrancheABQDmin2=$TrancheABQDmin1;
                $TrancheABQDmin1=$VQSRABVar[$kk];
            } elsif ($VQSRABVar[$kk]<$TrancheABQDmin2 && $VQSRABVar[$kk]>=0) {
                $TrancheABQDmin2=$VQSRABVar[$kk];
            }
            if ($VQSRAlignVar[$kk]<=$TrancheAlignmin1 && $VQSRAlignVar[$kk]>=0) {
                $TrancheAlignmin2=$TrancheAlignmin1;
                $TrancheAlignmin1=$VQSRAlignVar[$kk];
            } elsif ($VQSRAlignVar[$kk]<$TrancheAlignmin2 && $VQSRAlignVar[$kk]>=0) {
                $TrancheAlignmin2=$VQSRAlignVar[$kk];
            }
            if ($VQSRMapVar[$kk]<=$TrancheMapmin1 && $VQSRMapVar[$kk]>=0) {
                $TrancheMapmin2=$TrancheMapmin1;
                $TrancheMapmin1=$VQSRMapVar[$kk];
            } elsif ($VQSRMapVar[$kk]<$TrancheMapmin2 && $VQSRMapVar[$kk]>=0) {
                $TrancheMapmin2=$VQSRMapVar[$kk];
            }
        }
        $DPSum+=$DPVar[$kk];
    }
    
    my $j=0; my $PL0=0; my $PLminsum1=1000000; my $PLminsum2=1000000; my $PLTottxt="";
    foreach my $PL (@PLTot) {
        if ($j>=$ALTPL[$altno-1]) {last;}
        if ($j==0) {
            $PLTottxt=$PL;
        } else {
            $PLTottxt="$PLTottxt,$PL";
        }
        if ($PL<$PLminsum1) { #most likely net genotype
            $genosum=$j+1;
            $PLminsum2=$PLminsum1;
            $PLminsum1=$PL;
        } elsif ($PL<$PLminsum2) { $PLminsum2=$PL; }
        $j++;
    }
    my $PLminsum=$PLminsum2-$PLminsum1;
    my $PLbyDPsum=0;
    if ($DPSum>0) {$PLbyDPsum=$PLminsum/($DPSum+0.0);}
    my $DPSumall=$DPSum;
    my $PLminsumall=$PLminsum;
    my $PLbyDPsumall=$PLbyDPsum;
    my $PLTottxtall=$PLTottxt;
    my $PL0diffsumall=$PLTot[0]-$PLTot[$genosum-1];
    my $genosumall=$genosum;
    #if ($minchrompos>1000948929) { exit; }
    if ($minchrompos==1000948929) { 
    	print "All:$varType,$DPSum,$PLminsum,$PLbyDPsum,$MinPLbyDPSumSNP[$genosum-1],$YesPLTot[$genosum-1],$NoPLTot[$genosum-1],$HapCallVar\n";
    	if ($varType eq "SNP" && $PLminsum>$MinPLSumSNP[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumSNP[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)) { print "true\n"; }
    	if (!($varType eq "SNP") && $PLminsum>$MinPLSumINDEL[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumINDEL[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1))  { print "true\n"; }
    }
    if ($DPSum>0 && (($varType eq "SNP" && $PLminsum>$MinPLSumSNP[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumSNP[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)) || (!($varType eq "SNP") && $PLminsum>$MinPLSumINDEL[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumINDEL[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)))) {
        $filtergeno="PASS";
        if ($minchrompos==1000010285) { print "PASS\n"; }
    } else { #genotype cannot be determined from combined data, so try filtering methods
	        
        my @tranches = (99.5, 99.0, 95.0, 90.0);
        my @readlens = (40,60,80,90);

            #first filtering method starts with highest tranche, filters SSEs then alignment then readlength, and then tries lower tranches, but tests if datasets agree after each filter
        my $tt=0;
            my $nexttranche=1;
        while ($tt<4 && $nexttranche==1) {
            $nexttranche=0;
            @filteredgenos=((0) x 15); @mapgood=((0) x 15);
            $DPSum=0; $allPLannot="";
            @PLTot=((0) x 15); @NoPLTot=((0) x 15); @YesPLTot=((0) x 15);
            $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
            for (my $kk=0; $kk<$datasetcalls; $kk++) {
                if ($genoVar[$kk]>1 && $VQSRSSEVar[$kk]>$tranches[$tt]) {
                    $allPLannot="$allPLannot;PL${infiles[$datasetVar[$kk]]}=filteredSSE${VQSRSSEVar[$kk]},$PLtextVar[$kk]";
                    my @PLs = split( ",", $PLtextVar[$kk]);
                    my $j=0; my $PL0=0; my $PLmin=1000000;
                    foreach my $PL (@PLs) {
                        if ($PL>0 && $PL<$PLmin) { $PLmin=$PL; }
                        $j++;
                    }
                    if ($PLmin>20) { $filteredgenos[$genoVar[$kk]-1]++;}
                } else {
                
                    $allPLannot="$allPLannot;PL${infiles[$datasetVar[$kk]]}=$PLtextVar[$kk]";
                    
                    if ($PLtextVar[$kk] eq "") { next; } #no coverage
                    my @PLs = split( ",", $PLtextVar[$kk]);
                    my $j=0; my $PL0=0; my $PLmin=1000000;
                    foreach my $PL (@PLs) {
                        $PLTot[$j] += $PL;
                        if ($PL>20) { #confidently not this genotype in dataset
                            $NoPLTot[$j]+=1;
                        }
                        if ($PL>0 && $PL<$PLmin) { $PLmin=$PL; }
                        $j++;
                    }
                    if ($PLmin>20) {
                        $YesPLTot[$genoVar[$kk]-1]++;
                        if ($VQSRMapVar[$kk]<=95) {$mapgood[$genoVar[$kk]-1]++;}
                        if ($VQSRSSEVar[$kk]<=$TrancheSSEmin1 && $VQSRSSEVar[$kk]>=0) {
                            $TrancheSSEmin2=$TrancheSSEmin1;
                            $TrancheSSEmin1=$VQSRSSEVar[$kk];
                        } elsif ($VQSRSSEVar[$kk]<$TrancheSSEmin2 && $VQSRSSEVar[$kk]>=0) {
                            $TrancheSSEmin2=$VQSRSSEVar[$kk];
                        }
                        if ($VQSRABVar[$kk]<=$TrancheABQDmin1 && $VQSRABVar[$kk]>=0) {
                            $TrancheABQDmin2=$TrancheABQDmin1;
                            $TrancheABQDmin1=$VQSRABVar[$kk];
                        } elsif ($VQSRABVar[$kk]<$TrancheABQDmin2 && $VQSRABVar[$kk]>=0) {
                            $TrancheABQDmin2=$VQSRABVar[$kk];
                        }
                        if ($VQSRAlignVar[$kk]<=$TrancheAlignmin1 && $VQSRAlignVar[$kk]>=0) {
                            $TrancheAlignmin2=$TrancheAlignmin1;
                            $TrancheAlignmin1=$VQSRAlignVar[$kk];
                        } elsif ($VQSRAlignVar[$kk]<$TrancheAlignmin2 && $VQSRAlignVar[$kk]>=0) {
                            $TrancheAlignmin2=$VQSRAlignVar[$kk];
                        }
                        if ($VQSRMapVar[$kk]<=$TrancheMapmin1 && $VQSRMapVar[$kk]>=0) {
                            $TrancheMapmin2=$TrancheMapmin1;
                            $TrancheMapmin1=$VQSRMapVar[$kk];
                        } elsif ($VQSRMapVar[$kk]<$TrancheMapmin2 && $VQSRMapVar[$kk]>=0) {
                            $TrancheMapmin2=$VQSRMapVar[$kk];
                        }
                    }
                    $DPSum+=$DPVar[$kk];
                }
            }
            #if ($posVar[0]==28593) { print "$allPLannot\t$tt\n";}
                
             $j=0;  $PL0=0;  $PLminsum1=1000000;  $PLminsum2=1000000;  $PLTottxt="";
            foreach my $PL (@PLTot) {
                if ($j>=$ALTPL[$altno-1]) {last;}
                if ($j==0) {
                    $PLTottxt=$PL;
                } else {
                    $PLTottxt="$PLTottxt,$PL";
                }
                if ($PL<$PLminsum1) { #most likely net genotype
                    $genosum=$j+1;
                    $PLminsum2=$PLminsum1;
                    $PLminsum1=$PL;
                } elsif ($PL<$PLminsum2) { $PLminsum2=$PL; }
                $j++;
            }
             $PLminsum=$PLminsum2-$PLminsum1;
            my $PLbyDPsum=0;
            if ($DPSum>0) {$PLbyDPsum=$PLminsum/($DPSum+0.0);}
    if ($minchrompos==1000010285) { print "SSE,$tt:$DPSum,$PLminsum,$PLbyDPsum,$MinPLbyDPSumSNP[$genosum-1]\n"; }
            if ($DPSum>0 && (($varType eq "SNP" && $PLminsum>$MinPLSumSNP[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumSNP[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)) || (!($varType eq "SNP") && $PLminsum>$MinPLSumINDEL[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumINDEL[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)))) {
                    $filtergeno="PASS";
                } else { #genotype cannot be determined from filtering SSEs, so also filter alignment artifacts
                    @filteredgenos=((0) x 15); @mapgood=((0) x 15);
                    $DPSum=0; $allPLannot="";
                    @PLTot=((0) x 15); @NoPLTot=((0) x 15); @YesPLTot=((0) x 15);
                    $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
                    for (my $kk=0; $kk<$datasetcalls; $kk++) {
                        if (($genoVar[$kk]>1 && $VQSRSSEVar[$kk]>$tranches[$tt]) || $VQSRAlignVar[$kk]>$tranches[$tt]) {
                            $allPLannot="$allPLannot;PL${infiles[$datasetVar[$kk]]}=filteredSSE${VQSRSSEVar[$kk]}Align${VQSRAlignVar[$kk]},$PLtextVar[$kk]";
                            my @PLs = split( ",", $PLtextVar[$kk]);
                            my $j=0; my $PL0=0; my $PLmin=1000000;
                            foreach my $PL (@PLs) {
                                if ($PL>0 && $PL<$PLmin) { $PLmin=$PL; }
                                $j++;
                            }
                            if ($PLmin>20) { $filteredgenos[$genoVar[$kk]-1]++;}
                        } else {
                            
                            $allPLannot="$allPLannot;PL${infiles[$datasetVar[$kk]]}=$PLtextVar[$kk]";
                            
                            if ($PLtextVar[$kk] eq "") { next; } #no coverage
                            my @PLs = split( ",", $PLtextVar[$kk]);
                            my $j=0; my $PL0=0; my $PLmin=1000000;
                            foreach my $PL (@PLs) {
                                $PLTot[$j] += $PL;
                                if ($PL>20) { #confidently not this genotype in dataset
                                    $NoPLTot[$j]+=1;
                                }
                                if ($PL>0 && $PL<$PLmin) { $PLmin=$PL; }
                                $j++;
                            }
                            if ($PLmin>20) {
                                $YesPLTot[$genoVar[$kk]-1]++;
                                if ($VQSRMapVar[$kk]<=95) {$mapgood[$genoVar[$kk]-1]++;}
                                if ($VQSRSSEVar[$kk]<=$TrancheSSEmin1 && $VQSRSSEVar[$kk]>=0) {
                                    $TrancheSSEmin2=$TrancheSSEmin1;
                                    $TrancheSSEmin1=$VQSRSSEVar[$kk];
                                } elsif ($VQSRSSEVar[$kk]<$TrancheSSEmin2 && $VQSRSSEVar[$kk]>=0) {
                                    $TrancheSSEmin2=$VQSRSSEVar[$kk];
                                }
                                if ($VQSRABVar[$kk]<=$TrancheABQDmin1 && $VQSRABVar[$kk]>=0) {
                                    $TrancheABQDmin2=$TrancheABQDmin1;
                                    $TrancheABQDmin1=$VQSRABVar[$kk];
                                } elsif ($VQSRABVar[$kk]<$TrancheABQDmin2 && $VQSRABVar[$kk]>=0) {
                                    $TrancheABQDmin2=$VQSRABVar[$kk];
                                }
                                if ($VQSRAlignVar[$kk]<=$TrancheAlignmin1 && $VQSRAlignVar[$kk]>=0) {
                                    $TrancheAlignmin2=$TrancheAlignmin1;
                                    $TrancheAlignmin1=$VQSRAlignVar[$kk];
                                } elsif ($VQSRAlignVar[$kk]<$TrancheAlignmin2 && $VQSRAlignVar[$kk]>=0) {
                                    $TrancheAlignmin2=$VQSRAlignVar[$kk];
                                }
                                if ($VQSRMapVar[$kk]<=$TrancheMapmin1 && $VQSRMapVar[$kk]>=0) {
                                    $TrancheMapmin2=$TrancheMapmin1;
                                    $TrancheMapmin1=$VQSRMapVar[$kk];
                                } elsif ($VQSRMapVar[$kk]<$TrancheMapmin2 && $VQSRMapVar[$kk]>=0) {
                                    $TrancheMapmin2=$VQSRMapVar[$kk];
                                }
                            }
                            $DPSum+=$DPVar[$kk];
                        }
                    }
                        #if ($posVar[0]==28593) { print "$allPLannot\t$tt\n";}
                        
                    $j=0;  $PL0=0;  $PLminsum1=1000000;  $PLminsum2=1000000;  $PLTottxt="";
                    foreach my $PL (@PLTot) {
                        if ($j>=$ALTPL[$altno-1]) {last;}
                        if ($j==0) {
                            $PLTottxt=$PL;
                        } else {
                            $PLTottxt="$PLTottxt,$PL";
                        }
                        if ($PL<$PLminsum1) { #most likely net genotype
                            $genosum=$j+1;
                            $PLminsum2=$PLminsum1;
                            $PLminsum1=$PL;
                        } elsif ($PL<$PLminsum2) { $PLminsum2=$PL; }
                        $j++;
                    }
                    $PLminsum=$PLminsum2-$PLminsum1;
                    my $PLbyDPsum=0;
                    if ($DPSum>0) {$PLbyDPsum=$PLminsum/($DPSum+0.0);}
    if ($minchrompos==1000010285) { print "Align,$tt:$DPSum,$PLminsum,$PLbyDPsum,$MinPLbyDPSumSNP[$genosum-1]\n"; }
                    if ($DPSum>0 && (($varType eq "SNP" && $PLminsum>$MinPLSumSNP[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumSNP[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)) || (!($varType eq "SNP") && $PLminsum>$MinPLSumINDEL[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumINDEL[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)))) {
                            $filtergeno="PASS";
                    } else { #genotype cannot be determined from filtering this tranche, so try the next one
                                    $nexttranche=1;
                    }
                    
                }
                
                $tt++;
        }
    }




    #last, if there is no consensus, check if there is a consensus genotype with only filtering CG since it sometimes has a different representation of complex variants
    if (!($filtergeno eq "PASS")) {
            @filteredgenos=((0) x 15); @mapgood=((0) x 15);
            $DPSum=0; $allPLannot=";NoCG";
            @PLTot=((0) x 15); @NoPLTot=((0) x 15); @YesPLTot=((0) x 15);
            $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
    for (my $kk=0; $kk<$datasetcalls; $kk++) {
        #print "$chromVar[$kk],$posVar[$kk],$idVar[$kk],$refVar[$kk],$altVar[$kk],$genoVar[$kk],$VQSRABVar[$kk],$VQSRSSEVar[$kk],$VQSRAlignVar[$kk],$VQSRMapVar[$kk],$datasetcalls,$PLtextVar[$kk]\n";
        #exit;
        $allPLannot="$allPLannot;PL${infiles[$datasetVar[$kk]]}=$PLtextVar[$kk]";
        if ($PLtextVar[$kk] eq "") { next; } #no coverage
        if ($platformsVar[$kk] eq "cg") { next; } #skip CG
        my @PLs = split( ",", $PLtextVar[$kk]);
        my $j=0; my $PL0=0; my $PLmin=1000000;
        foreach my $PL (@PLs) {
            $PLTot[$j] += $PL;
            if ($PL>20) { #confidently not this genotype in dataset
                $NoPLTot[$j]+=1;
            } 
            if ($PL>0 && $PL<$PLmin) { $PLmin=$PL; } #find the minimum PL>0, since this is the likelihood ratio of the most likely to the next most likely genotype
            $j++;
        }
        if ($PLmin>20) {
            #add this platform if it is a new platform calling this genotype
            if ($platformgenotxt[$genoVar[$kk]-1] =~ /,$platformsVar[$kk],/) {} else {
                $platformgenotxt[$genoVar[$kk]-1] = "$platformgenotxt[$genoVar[$kk]-1]$platformsVar[$kk],";
            }
            
            $YesPLTot[$genoVar[$kk]-1]++;
            if ($VQSRMapVar[$kk]<=95) {$mapgood[$genoVar[$kk]-1]++;}
            if ($VQSRSSEVar[$kk]<=$TrancheSSEmin1 && $VQSRSSEVar[$kk]>=0) {
                $TrancheSSEmin2=$TrancheSSEmin1;
                $TrancheSSEmin1=$VQSRSSEVar[$kk];
            } elsif ($VQSRSSEVar[$kk]<$TrancheSSEmin2 && $VQSRSSEVar[$kk]>=0) {
                $TrancheSSEmin2=$VQSRSSEVar[$kk];
            }
            if ($VQSRABVar[$kk]<=$TrancheABQDmin1 && $VQSRABVar[$kk]>=0) {
                $TrancheABQDmin2=$TrancheABQDmin1;
                $TrancheABQDmin1=$VQSRABVar[$kk];
            } elsif ($VQSRABVar[$kk]<$TrancheABQDmin2 && $VQSRABVar[$kk]>=0) {
                $TrancheABQDmin2=$VQSRABVar[$kk];
            }
            if ($VQSRAlignVar[$kk]<=$TrancheAlignmin1 && $VQSRAlignVar[$kk]>=0) {
                $TrancheAlignmin2=$TrancheAlignmin1;
                $TrancheAlignmin1=$VQSRAlignVar[$kk];
            } elsif ($VQSRAlignVar[$kk]<$TrancheAlignmin2 && $VQSRAlignVar[$kk]>=0) {
                $TrancheAlignmin2=$VQSRAlignVar[$kk];
            }
            if ($VQSRMapVar[$kk]<=$TrancheMapmin1 && $VQSRMapVar[$kk]>=0) {
                $TrancheMapmin2=$TrancheMapmin1;
                $TrancheMapmin1=$VQSRMapVar[$kk];
            } elsif ($VQSRMapVar[$kk]<$TrancheMapmin2 && $VQSRMapVar[$kk]>=0) {
                $TrancheMapmin2=$VQSRMapVar[$kk];
            }
        }
        $DPSum+=$DPVar[$kk];
    }
    
    $j=0; $PL0=0; $PLminsum1=1000000; $PLminsum2=1000000; $PLTottxt="";
    foreach my $PL (@PLTot) {
        if ($j>=$ALTPL[$altno-1]) {last;}
        if ($j==0) {
            $PLTottxt=$PL;
        } else {
            $PLTottxt="$PLTottxt,$PL";
        }
        if ($PL<$PLminsum1) { #most likely net genotype
            $genosum=$j+1;
            $PLminsum2=$PLminsum1;
            $PLminsum1=$PL;
        } elsif ($PL<$PLminsum2) { $PLminsum2=$PL; }
        $j++;
    }
    $PLminsum=$PLminsum2-$PLminsum1;
    $PLbyDPsum=0;
    if ($DPSum>0) {$PLbyDPsum=$PLminsum/($DPSum+0.0);}
    $DPSumall=$DPSum;
    $PLminsumall=$PLminsum;
    $PLbyDPsumall=$PLbyDPsum;
    $PLTottxtall=$PLTottxt;
    $PL0diffsumall=$PLTot[0]-$PLTot[$genosum-1];
    $genosumall=$genosum;
    #if ($minchrompos>1000948929) { exit; }
    if ($minchrompos==1000948929) { 
    	print "All:$varType,$DPSum,$PLminsum,$PLbyDPsum,$MinPLbyDPSumSNP[$genosum-1],$YesPLTot[$genosum-1],$NoPLTot[$genosum-1],$HapCallVar\n";
    	print "HapCallVar:$HapCallVar\n";
    	if ($varType eq "SNP" && $PLminsum>$MinPLSumSNP[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumSNP[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)) { print "true\n"; }
    	if (!($varType eq "SNP") && $PLminsum>$MinPLSumINDEL[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumINDEL[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1))  { print "true\n"; }
    }
    if ($DPSum>0 && (($varType eq "SNP" && $PLminsum>$MinPLSumSNP[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumSNP[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)) || (!($varType eq "SNP") && $PLminsum>$MinPLSumINDEL[$genosum-1] && $PLbyDPsum>$MinPLbyDPSumINDEL[$genosum-1] && $YesPLTot[$genosum-1]>$MINDATASETS && ($NoPLTot[$genosum-1]==0 || $YesPLTot[$genosum-1]/($NoPLTot[$genosum-1]+0.0)>=$MINYESNORATIO) && ($HapNoVar<2 || $genosum==1)))) {
        $filtergeno="PASS";
        if ($minchrompos==1000010285) { print "PASS\n"; }
    } 
    
    } #end if for removing CG calls




    #if ($posVar[0]==28593) { exit;}
    
    if ($minchrompos==1000948929) { 
    	print "after arbitrate\n";
	}
    #count number of platforms supporting called genotype
    my @platformcounts=split(",",$platformgenotxt[$genosum-1]);
    my $platformcount=$#platformcounts;
    my $platformtxt="none";
    if ($#platformcounts>0) {$platformtxt="$platformcounts[1]";}
    for (my $jj=2;$jj<=$#platformcounts;$jj++) {
        $platformtxt="$platformtxt,$platformcounts[$jj]";
    }
    
    my $PL0diffsum=$PLTot[0]-$PLTot[$genosum-1];

    if ($DPSum==0) {
        $PLbyDPsum=$PLbyDPsumall;
        $DPSum=$DPSumall;
        $PLTottxt=$PLTottxtall;
        $PLminsum=$PLminsumall;
        $PL0diffsum=$PL0diffsumall;
        $genosum=$genosumall;
    }

    if ($filtergeno eq "PASS") {
        my $filtertxt="";
        if ($multiplealts==1) { $filtertxt="Multiplealts"; }
        if (($genosum!=1 && $genosum!=3 && $genosum!=6 && $genosum!=10 && $genosum!=15) && $mapgood[$genosum-1]<$MAPGOODMIN) {
            $filtertxt="${filtertxt}MapGood$mapgood[$genosum-1]";
        }
        if (($genosum!=1 && $genosum!=3 && $genosum!=6 && $genosum!=10 && $genosum!=15) && $TrancheMapmin2>=95) {
            $filtertxt="${filtertxt}Map$TrancheMapmin2";
        }
        if ((($genosum==3 || $genosum==6 || $genosum==10 || $genosum==15) && $TrancheSSEmin2>=99) || (($genosum!=1 && $genosum!=3 && $genosum!=6 && $genosum!=10 && $genosum!=15) && $TrancheSSEmin2>=95)) {
            $filtertxt="${filtertxt}SSE$TrancheSSEmin2";
        }
        if (($genosum==1 && $TrancheAlignmin2>=99.5) || (($genosum==3 || $genosum==6 || $genosum==10 || $genosum==15 || (!($varType eq "SNP") && $genosum>1)) && $TrancheAlignmin2>=99) || (($genosum!=1 && $genosum!=3 && $genosum!=6 && $genosum!=10 && $genosum!=15 && ($varType eq "SNP")) && $TrancheAlignmin2>=95)) {
            $filtertxt="${filtertxt}Align$TrancheAlignmin2";
        }
        if ((($genosum==3 || $genosum==6 || $genosum==10 || $genosum==15 || !($varType eq "SNP")) && $TrancheABQDmin2>=99) || (($genosum!=1 && $genosum!=3 && $genosum!=6 && $genosum!=10 && $genosum!=15 && ($varType eq "SNP")) && $TrancheABQDmin2>=95)) {
            $filtertxt="${filtertxt}ABQD$TrancheABQDmin2";
        }
        if ($minchrompos<=$chromposdelregion) {
            if ($HapNoVar>=2 || $genosum==1) { # don't filter if it looks like HomRef or HaplotypeCaller doesn't call a variant here, since it's likely HomRef
            	$genosum=1;
            	$filtertxt="";
            } elsif ($HapNoVar>0) {
            	$filtertxt="InsideDel";
            }
        }
        if ($varType eq "INDEL" && $HapCallVar<2 && $genosum>1) { #filter indels that were not called by at least 2 HC datasets
        	$filtertxt = "IndelNoHCcall";
        }
        
        if (!($filtertxt eq "") && !($HapNoVar>=2 && $genosum==1)) {
            $filtergeno="filter=$filtertxt";
        } 
    } elsif ($HapNoVar>=($#infilehandlesDistHomAB/2)) {
        $genosum=1;
        $filtergeno="PASS"; #if haplotypecaller doesn't call a variant in more than half the datasets, then likely homref
    } elsif ($HapNoVar>=2) {
        if ($genosum==1) {
        	$filtergeno="PASS"; #if haplotypecaller doesn't call a variant and the net genotype is homref, then likely homref
        } else {
        	$filtergeno="filter=HapNoVar";
        }
    } elsif ($datasetcalls==2) {
        $filtergeno="filter=${datasetcalls}Dataset(s)";
    } elsif ($DPSum<20) {
        $filtergeno="filter=Cov$DPSum";
    } else {
        $filtergeno="filter=ConflictingGenotypes";
    }
    
    my $filterReason=""; my $filtfield="";
    
    if ($filtergeno eq "PASS") {
        $filtfield="PASS";
    } else {
        $filtfield="Uncertain";
        $filterReason="$filtergeno;";        
    }
    
    #keep only ALT alleles that are used in confident genotypes
    my $altsall=$altVarLong;
    if ($filtergeno eq "PASS") {
        for (my $kk=0; $kk<$datasetcalls; $kk++) {
            if ($PLtextVar[$kk] eq "") { next; } #no coverage o
            my @PLs = split( ",", $PLtextVar[$kk]);
            if ($genosum<4) { #first alt only
            	$PLtextVar[$kk]="$PLs[0],$PLs[1],$PLs[2]";
            } elsif ($genosum==4 || $genosum==6) { #second alt only
            	$PLtextVar[$kk]="$PLs[0],$PLs[3],$PLs[5]";
            } elsif ($genosum==7 || $genosum==10) { #third alt only
            	$PLtextVar[$kk]="$PLs[0],$PLs[6],$PLs[9]";
            } elsif ($genosum==11 || $genosum==15) { #fourth alt only
            	$PLtextVar[$kk]="$PLs[0],$PLs[10],$PLs[14]";
            } 
        }
            my @alts = split( ",", $altVarLong);
            my @PLs = split( ",", $PLTottxt);
            if ($genosum<4) { #first alt only
            	$altVarLong="$alts[0]";
            	$PLTottxt="$PLs[0],$PLs[1],$PLs[2]";
            } elsif ($genosum==4 || $genosum==6) { #second alt only
            	$altVarLong="$alts[1]";
            	$PLTottxt="$PLs[0],$PLs[3],$PLs[5]";
            	if ($genosum==4) {
            		$genosum=2;
            	} else {
            		$genosum=3;
            	}
            } elsif ($genosum==7 || $genosum==10) { #third alt only
            	$altVarLong="$alts[2]";
            	$PLTottxt="$PLs[0],$PLs[6],$PLs[9]";
            	if ($genosum==7) {
            		$genosum=2;
            	} else {
            		$genosum=3;
            	}
            } elsif ($genosum==11 || $genosum==15) { #fourth alt only
            	$altVarLong="$alts[3]";
            	$PLTottxt="$PLs[0],$PLs[10],$PLs[14]";
            	if ($genosum==11) {
            		$genosum=2;
            	} else {
            		$genosum=3;
            	}
            } 
	}

    #add annotation to locations that have platform-specific errors
    my @platgood = ((0) x $plats); #count accurate calls on each platform
    my @platbad = ((0) x $plats); #count inaccurate calls on each platform
    my $platbias = "none";
    if ($filtergeno eq "PASS" && $platformcount>2) {
        for (my $kk=0; $kk<$datasetcalls; $kk++) {
            if ($PLtextVar[$kk] eq "") { next; } #no coverage o
            my @PLs = split( ",", $PLtextVar[$kk]);
            my $j=0; my $PL0=0; my $PLmin=1000000;
            foreach my $PL (@PLs) {
                $PLTot[$j] += $PL;
                if ($PL>20) { #confidently not this genotype in dataset
                    $NoPLTot[$j]+=1;
                }
                if ($PL>0 && $PL<$PLmin) { $PLmin=$PL; } #find the minimum PL>0, since this is the likelihood ratio of the most likely to the next most likely genotype
                $j++;
            }
            if ($PLmin>20) {
                if ($genoVar[$kk] == $genosum) {
                    my $j=0; #find which platform it is
                    for (my $i=0; $i<$plats; $i++) {
                        if ($platformsVar[$kk] eq $platnames[$i]) { $j=$i; }
                    }
                    $platgood[$j]++;
                } else {
                    my $j=0; #find which platform it is
                    for (my $i=0; $i<$plats; $i++) {
                        if ($platformsVar[$kk] eq $platnames[$i]) { $j=$i; }
                    }
                    $platbad[$j]++;

                }
            }
        }
        for (my $i=0; $i<$plats; $i++) {
            if (($platbad[$i]>0 && $platgood[$i]==0) || $platbad[$i]>2) {
                if ($platbias eq "none") {
                    $platbias = $platnames[$i];
                } else {
                    $platbias = "$platbias,$platnames[$i]";
                }
            }
        }

    }
    if ($chromVar[0] eq "0" || $DPSum==0) { next; } #skip output if no datasets exist at this row
    my $hruntxt = ""; my $RUtxt = ""; my $RPAtxt = "";
    if ($hrun>-1) { $hruntxt=";HRun=$hrun"; }
    if (!($RU eq "")) { $RUtxt=";RU=$RU"; }
    if (!($RPA eq "")) { $RPAtxt=";RPA=$RPA"; }

    if ($minchrompos==1000948929) { 
    	print "before out\n";
	}
    $PLbyDPsum=int(($PLminsum/$DPSum*100)+0.5)/100.0;
            print OUTPUT "$chromVar[0]\t$posVar[0]\t$idVar[0]\t$refVarLong\t$altVarLong\t$PL0diffsum\t$filtfield\t${filterReason}geno=$genosum;platforms=$platformcount;platformnames=$platformtxt;platformbias=$platbias;varType=$varType;YesPLtot=$YesPLTot[$genosum-1];NoPLTot=$NoPLTot[$genosum-1];HapNoVar=$HapNoVar;PLminsum=$PLminsum;DPSum=$DPSum;PLminsumOverDP=$PLbyDPsum;genoMapGood=$mapgood[$genosum-1];TrancheSSEmin2=$TrancheSSEmin2;TrancheABQDmin2=$TrancheABQDmin2;TrancheAlignmin2=$TrancheAlignmin2;TrancheMapmin2=$TrancheMapmin2;datasetcalls=$datasetcalls;allalts=$altsall$allPLannot$hruntxt$RUtxt$RPAtxt\tGT:DP:GQ:PL\t$GTs[$genosum-1]:$DPSum:$PLminsum:$PLTottxt\n";
            if ($genosum==1) {
                print OUTHOMREF "$chromVar[0]\t$posVar[0]\t$idVar[0]\t$refVarLong\t$altVarLong\t$PL0diffsum\t$filtfield\t${filterReason}geno=$genosum;platforms=$platformcount;platformnames=$platformtxt;platformbias=$platbias;varType=$varType;YesPLtot=$YesPLTot[$genosum-1];NoPLTot=$NoPLTot[$genosum-1];HapNoVar=$HapNoVar;PLminsum=$PLminsum;DPSum=$DPSum;PLminsumOverDP=$PLbyDPsum;genoMapGood=$mapgood[$genosum-1];TrancheSSEmin2=$TrancheSSEmin2;TrancheABQDmin2=$TrancheABQDmin2;TrancheAlignmin2=$TrancheAlignmin2;TrancheMapmin2=$TrancheMapmin2;datasetcalls=$datasetcalls;allalts=$altsall$allPLannot$hruntxt$RUtxt$RPAtxt\tGT:DP:GQ:PL\t$GTs[$genosum-1]:$DPSum:$PLminsum:$PLTottxt\n";
            } elsif ($genosum==3 || $genosum==6 || $genosum==10 || $genosum==15) {
                print OUTHOMVAR "$chromVar[0]\t$posVar[0]\t$idVar[0]\t$refVarLong\t$altVarLong\t$PL0diffsum\t$filtfield\t${filterReason}geno=$genosum;platforms=$platformcount;platformnames=$platformtxt;platformbias=$platbias;varType=$varType;YesPLtot=$YesPLTot[$genosum-1];NoPLTot=$NoPLTot[$genosum-1];HapNoVar=$HapNoVar;PLminsum=$PLminsum;DPSum=$DPSum;PLminsumOverDP=$PLbyDPsum;genoMapGood=$mapgood[$genosum-1];TrancheSSEmin2=$TrancheSSEmin2;TrancheABQDmin2=$TrancheABQDmin2;TrancheAlignmin2=$TrancheAlignmin2;TrancheMapmin2=$TrancheMapmin2;datasetcalls=$datasetcalls;allalts=$altsall$allPLannot$hruntxt$RUtxt$RPAtxt\tGT:DP:GQ:PL\t$GTs[$genosum-1]:$DPSum:$PLminsum:$PLTottxt\n";
            } else {
                print OUTHET "$chromVar[0]\t$posVar[0]\t$idVar[0]\t$refVarLong\t$altVarLong\t$PL0diffsum\t$filtfield\t${filterReason}geno=$genosum;platforms=$platformcount;platformnames=$platformtxt;platformbias=$platbias;varType=$varType;YesPLtot=$YesPLTot[$genosum-1];NoPLTot=$NoPLTot[$genosum-1];HapNoVar=$HapNoVar;PLminsum=$PLminsum;DPSum=$DPSum;PLminsumOverDP=$PLbyDPsum;genoMapGood=$mapgood[$genosum-1];TrancheSSEmin2=$TrancheSSEmin2;TrancheABQDmin2=$TrancheABQDmin2;TrancheAlignmin2=$TrancheAlignmin2;TrancheMapmin2=$TrancheMapmin2;datasetcalls=$datasetcalls;allalts=$altsall$allPLannot$hruntxt$RUtxt$RPAtxt\tGT:DP:GQ:PL\t$GTs[$genosum-1]:$DPSum:$PLminsum:$PLTottxt\n";
            }
    if ($genosum!=1 || !($filtergeno eq "PASS")) {
        print OUTVARALL "$chromVar[0]\t$posVar[0]\t$idVar[0]\t$refVarLong\t$altVarLong\t$PL0diffsum\t$filtfield\t${filterReason}geno=$genosum;platforms=$platformcount;platformnames=$platformtxt;platformbias=$platbias;varType=$varType;YesPLtot=$YesPLTot[$genosum-1];NoPLTot=$NoPLTot[$genosum-1];HapNoVar=$HapNoVar;PLminsum=$PLminsum;DPSum=$DPSum;PLminsumOverDP=$PLbyDPsum;genoMapGood=$mapgood[$genosum-1];TrancheSSEmin2=$TrancheSSEmin2;TrancheABQDmin2=$TrancheABQDmin2;TrancheAlignmin2=$TrancheAlignmin2;TrancheMapmin2=$TrancheMapmin2;datasetcalls=$datasetcalls;allalts=$altsall$allPLannot$hruntxt$RUtxt$RPAtxt\tGT:DP:GQ:PL\t$GTs[$genosum-1]:$DPSum:$PLminsum:$PLTottxt\n";
        if ($filtergeno eq "PASS") {
            $chromposdelregion=$minchrompos+length($refVarLong)-1; #specify region this REF covers so that we consider other variants inside it to be uncertain
            print OUTVARPASS "$chromVar[0]\t$posVar[0]\t$idVar[0]\t$refVarLong\t$altVarLong\t$PL0diffsum\t$filtfield\t${filterReason}geno=$genosum;platforms=$platformcount;platformnames=$platformtxt;platformbias=$platbias;varType=$varType;YesPLtot=$YesPLTot[$genosum-1];NoPLTot=$NoPLTot[$genosum-1];HapNoVar=$HapNoVar;PLminsum=$PLminsum;DPSum=$DPSum;PLminsumOverDP=$PLbyDPsum;genoMapGood=$mapgood[$genosum-1];TrancheSSEmin2=$TrancheSSEmin2;TrancheABQDmin2=$TrancheABQDmin2;TrancheAlignmin2=$TrancheAlignmin2;TrancheMapmin2=$TrancheMapmin2;datasetcalls=$datasetcalls;allalts=$altsall$allPLannot$hruntxt$RUtxt$RPAtxt\tGT:DP:GQ:PL\t$GTs[$genosum-1]:$DPSum:$PLminsum:$PLTottxt\n";
        }
    }
        
 	$i++;
}




close OUTPUT;
    close OUTHOMREF;
    close OUTHET;
    close OUTHOMVAR;
close OUTVARPASS;
close OUTVARALL;

