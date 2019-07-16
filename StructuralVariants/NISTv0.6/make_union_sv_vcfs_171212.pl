#To load following perl modules
use strict;
use warnings;
use Getopt::Long;


#Declare variables
my $VTPATH="/Applications/bioinfo/vt/";
my $BCFTOOLSPATH="/Applications/bioinfo/bcftools/";
my $BEDTOOLSPATH="/Applications/bioinfo/bedtools2.26.0/";

#create header
my $header = "##fileformat=VCFv4.2\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth summed across all datasets, excluding MQ0 reads\">\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Net Genotype quality across all datasets, calculated from GQ scores of callsets supporting the consensus GT, using only one callset from each dataset\">\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Net allele depths across all unfiltered datasets with called genotype\">\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Consensus Genotype across all datasets with called genotype\">\
##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase set in which this variant falls\">\
##fileDate=20171212\
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AJUNION\n";

#To check input files
if ($#ARGV < 0) {print usage();die "\n\nPlease specify the path of the vcf files.\n\n"}
#Input parameters
my $VCFLIST=$ARGV[0];
my $FILEPATH="";
if ($#ARGV > 0) {$FILEPATH=$ARGV[1];}

#To define usage of the preprocess_combine_vcfs
print "Usage = perl make_union_sv_vcfs.pl $VCFLIST [$FILEPATH] \n\n";

sub usage {
	print "\nusage: perl make_union_sv_vcfs.pl <VCFLIST> [<FILEPATH>] \n";
	print "\tVCFLIST \t\t\tThe input Tab-delimited text file containing 2 columns: vcf file name and ID prefix for ID column \n";
	print "\tFILEPATH\t\t\t\tPath where vcf and bed files are located with trailing \\\n";
	exit 1;
}

my $m = 1; #vcf number in order in VCFLIST
my $vcfs=""; #output variable that is printed that will contain names of processed vcfs for integration


my $delbeds="";
my $delbednames="";

if (1==0) {
	#open output union vcf 
	open (OUTUNION, ">", "${FILEPATH}union_171212.vcf") or die $!;
open (VCFLIST, "<", $VCFLIST) or die "Can't open the input vcf list file. Please check if the file exists\n\n";


#loop through each line in callset table
while (my $line = <VCFLIST>) {
	my $vcfno="$m";
	if ($m<10) {
		$vcfno="0$m";
	}
	chomp($line); #removes new line char
	
	my @line_data = split("\t", $line);
    
    if($line =~ m/^\#/) { #get header line and go to next line
        next;
    }
    
    my $vcfstart=""; #begininning of vcf file name
    if($line_data[0] =~ /(.*)\.vcf.gz/) { $vcfstart=$1; } else { print "vcf does not end in .vcf.gz in input table\n"; exit 1; }
    
    #unzip vcf, keep first 10 columns, removes hom ref lines, removes chr prefix, breaks into "primitive" SNPs and indels, and sorts
    if ($line_data[1] =~ /vcfBeta/) { 
    	`gunzip -c ${FILEPATH}${vcfstart}.vcf.gz | awk '{if((length(\$4)-length(\$5)>19 || length(\$4)-length(\$5)<-19 || (length(\$4)>19 && length(\$5)>19)) && !(\$5 ~ /\\[/) && !(\$5 ~ /\\]/)) print}'  | grep -v '\\./\\.' | grep -v '\\.:\\.\$' > ${FILEPATH}intermediate_171212/${vcfstart}_refalt.vcf`;  #select rows with REF and ALT sequence and likely indel>19bp (not perfect for multi-allelic sites
    } else {
    	`zgrep -v '^MT' ${FILEPATH}${vcfstart}.vcf.gz | awk '{if(length(\$4)-length(\$5)>19 || length(\$4)-length(\$5)<-19 || (length(\$4)>19 && length(\$5)>19)) print}' | grep -v '^#' | grep -v '^NC' | grep -v '^GL' | grep -v '^hs' > ${FILEPATH}intermediate_171212/${vcfstart}_refalt.vcf`;  #select rows with REF and ALT sequence and likely indel>19bp (not perfect for multi-allelic sites
    }
    `cat header.vcf ${FILEPATH}intermediate_171212/${vcfstart}_refalt.vcf | sed 's/^chr//' > ${FILEPATH}intermediate_171212/${vcfstart}_refalt_head.vcf`;
    `${BCFTOOLSPATH}bcftools norm -D -c x -f /Users/jzook/Documents/references/human_g1k_v37.fasta ${FILEPATH}intermediate_171212/${vcfstart}_refalt_head.vcf | ${VTPATH}vt decompose -s - | awk '!(\$5 ~ /\\./)' | ${BCFTOOLSPATH}bcftools norm -D  - > ${FILEPATH}intermediate_171212/${vcfstart}_decomposed_norm.vcf`; #split, deduplicate, and normalize alleles, and exclude incorrect ref alleles
#    `${BCFTOOLSPATH}bcftools norm -D -c x -f /Users/jzook/Documents/references/human_g1k_v37.fasta ${FILEPATH}intermediate_171212/${vcfstart}_decomposed.vcf > ${FILEPATH}intermediate_171212/${vcfstart}_decomposed_norm.vcf`; #normalize alleles and exclude incorrect ref alleles
    `awk '{if(length(\$4)-length(\$5)>19 || length(\$4)-length(\$5)<-19 || (length(\$4)>19 && length(\$5)>19)) print}'  ${FILEPATH}intermediate_171212/${vcfstart}_decomposed_norm.vcf > ${FILEPATH}intermediate_171212/${vcfstart}_refalt_decomposed_norm_indelgt19.vcf`; #select indels>19bp
    `gunzip -c ${FILEPATH}${vcfstart}.vcf.gz | awk '{if((\$5 ~ /DEL/ || \$5 ~ /INS/ || \$5 ~ /MEI/ || \$5 ~ /INV/ || \$5 ~ /DUP/) && (\$7 ~ /\./ || \$7 ~ /PASS/)) print}'  > ${FILEPATH}intermediate_171212/${vcfstart}_imprecise.vcf`;  #select SV rows with inexact sequence


   
    #open formatted vcf from above with REF and ALT
    open (VCFIN, "<", "${FILEPATH}intermediate_171212/${vcfstart}_refalt_decomposed_norm_indelgt19.vcf") or die "Can't open the vcf ${FILEPATH}intermediate_171212/${vcfstart}_refalt_decomposed_norm_indelgt19.vcf. Please check if the file exists\n\n";
	#open output vcf with additional formatting
	open (OUTINDIV, ">", "${FILEPATH}intermediate_171212/${vcfstart}_withID_gt19bp.vcf") or die $!;
    
    my $i=1;
    while (my $line1 = <VCFIN>) {
		chomp($line1);
		my @line_data1 = split("\t", $line1);
	
		if($line1 =~ m/^\#/) { #ignore header lines
		} else { #otherwise, print the revised line with ID and no chr prefix
			my $chrom = $line_data1[0];
			if ($chrom =~ /chr(.*)/ || $chrom =~ /M/) { $chrom=$1; }
			if (length($chrom)>2) {next;} #skip lines that aren't 1-22, X or Y
			print OUTINDIV "$chrom\t$line_data1[1]\t${line_data[1]}_$i\t$line_data1[3]\t$line_data1[4]\t$line_data1[5]\t$line_data1[6]\t$line_data1[7]\t$line_data1[8]\t$line_data1[9]\n";
			print OUTUNION "$chrom\t$line_data1[1]\t${line_data[1]}_$i\t$line_data1[3]\t$line_data1[4]\t$line_data1[5]\t$line_data1[6]\t$line_data1[7]\t$line_data1[8]\t$line_data1[9]\n";
			$i=$i+1;
		}
	}
	close VCFIN;
    #open formatted vcf from above with imprecise SVs
    open (VCFIN, "<", "${FILEPATH}intermediate_171212/${vcfstart}_imprecise.vcf") or die "Can't open the vcf ${FILEPATH}intermediate_171212/${vcfstart}_imprecise.vcf. Please check if the file exists\n\n";
	open (OUTIMPRECISEBED, ">", "${FILEPATH}intermediate_171212/${vcfstart}_imprecise.bed") or die $!;
	open (OUTDELBED, ">", "${FILEPATH}intermediate_171212/${vcfstart}_imprecise_del.bed") or die $!;
    while (my $line1 = <VCFIN>) {
		chomp($line1);
		my @line_data1 = split("\t", $line1);
	
		if($line1 =~ m/^\#/) { #ignore header lines
		} else { #otherwise, print the revised line with ID and no chr prefix
			my $chrom = $line_data1[0];
			if ($chrom =~ /chr(.*)/) { $chrom=$1; }
			if (length($chrom)>2 || $chrom =~ /M/) {next;} #skip lines that aren't 1-22, X or Y
			if ($#line_data1<9) {
				print OUTINDIV "$chrom\t$line_data1[1]\t${line_data[1]}_$i\t$line_data1[3]\t$line_data1[4]\t$line_data1[5]\t$line_data1[6]\t$line_data1[7]\t.\t.\n";
				print OUTUNION "$chrom\t$line_data1[1]\t${line_data[1]}_$i\t$line_data1[3]\t$line_data1[4]\t$line_data1[5]\t$line_data1[6]\t$line_data1[7]\t.\t.\n";
			} else {
				print OUTINDIV "$chrom\t$line_data1[1]\t${line_data[1]}_$i\t$line_data1[3]\t$line_data1[4]\t$line_data1[5]\t$line_data1[6]\t$line_data1[7]\t$line_data1[8]\t$line_data1[9]\n";
				print OUTUNION "$chrom\t$line_data1[1]\t${line_data[1]}_$i\t$line_data1[3]\t$line_data1[4]\t$line_data1[5]\t$line_data1[6]\t$line_data1[7]\t$line_data1[8]\t$line_data1[9]\n";
			}
			my $endpos= $line_data1[1]+1;
			my $svlen= 0;
			my @infoflds=split( ";", $line_data1[7]);
			for (my $k=0;$k<=$#infoflds;$k++) {
				my @infofld=split( "=", $infoflds[$k]);
				if ($#infofld==1 && $infofld[0] eq "END" && $infofld[1]>$endpos) { 
					$endpos=$infofld[1];
				}
				if ($#infofld==1 && $infofld[0] eq "SVLEN" ) { 
					$svlen=$infofld[1];
				}
			}
			$i=$i+1;
			my $startbed=$line_data1[1]-1000;
			if ($startbed<0) {$startbed=0;}
			my $endbed=$endpos+1000;
			if ($endbed-$startbed<300000) {print OUTIMPRECISEBED "$chrom\t$startbed\t$endbed\n";}

			my $startdel=$line_data1[1]-1;
			my $enddel=$endpos;
			if ($svlen>0) {$enddel=$startdel+$svlen;}
			if ($enddel-$startdel<300000 && $enddel-$startdel>19 && ($line1 =~ m/DEL/ || $line1 =~ m/eletion/)) {print OUTDELBED "$chrom\t$startdel\t$enddel\n";}

		}
	}
	close VCFIN;
	close OUTINDIV;
	close OUTIMPRECISEBED;
	close OUTDELBED;

	#create sorted vcf with IDs
	print `sed 's/^X/23/;s/^Y/24/' ${FILEPATH}intermediate_171212/${vcfstart}_withID_gt19bp.vcf | sort -k1,1n -k2,2n -k4,4 -k5,5 | sed 's/^23/X/;s/^24/Y/' > ${FILEPATH}intermediate_171212/${vcfstart}_withID_gt19bp.sort.nohead.vcf`;
	print `cat header.vcf ${FILEPATH}intermediate_171212/${vcfstart}_withID_gt19bp.sort.nohead.vcf  > ${FILEPATH}intermediate_171212/${vcfstart}_withID_gt19bp.sort.vcf`;
	
    #create bed file for input to svrefine
    `awk '{FS="\t";OFS="\t"} {print \$1, \$2-1000, \$2+length(\$4)+1000}' ${FILEPATH}intermediate_171212/${vcfstart}_refalt_decomposed_norm_indelgt19.vcf | sed 's/chr//' > ${FILEPATH}intermediate_171212/${vcfstart}_refalt.bed`;

    #create bed file for deletions
    `awk '{FS="\t";OFS="\t"} {if(length(\$4)-length(\$5)>19) print \$1, \$2, \$2+length(\$4)-length(\$5)}' ${FILEPATH}intermediate_171212/${vcfstart}_refalt.vcf | sed 's/chr//' > ${FILEPATH}intermediate_171212/${vcfstart}_refalt_del.bed`;
	print `cat ${FILEPATH}intermediate_171212/${vcfstart}_refalt_del.bed ${FILEPATH}intermediate_171212/${vcfstart}_imprecise_del.bed | sed 's/^chr//;s/X/23/;s/Y/24/;s/^MT/25/;s/.*://' | grep -v GL | grep -v hap | grep -v ctg | grep -v NC | grep -v hs | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/;s/^25/MT/' | ${BEDTOOLSPATH}bin/mergeBed -i stdin > ${FILEPATH}intermediate_171212/${vcfstart}_all_del.bed`;
	$delbeds="$delbeds ${FILEPATH}intermediate_171212/${vcfstart}_all_del.bed";
	$delbednames="$delbednames ${line_data[1]}";
	$m=$m+1;

}
close VCFLIST;
close OUTUNION;

}
print `sed 's/^X/23/;s/^Y/24/;s/^MT/25/' ${FILEPATH}union_171212.vcf | sort -k1,1n -k2,2n -k4,4 -k5,5 | sed 's/^23/X/;s/^24/Y/;s/^25/MT/' > ${FILEPATH}union_171212.sort.vcf`;  #sort union vcf
print `cat header.vcf ${FILEPATH}union_171212.sort.vcf > ${FILEPATH}union_171212.sort.head.vcf`;
print `awk '{if(\$1 ~ /^#/ || length(\$4)-length(\$5)>49 || length(\$4)-length(\$5)<-49 || (length(\$4)>49 && length(\$5)>49) || \$5 ~ /DEL/ || \$5 ~ /INS/ || \$5 ~ /MEI/ || \$5 ~ /INV/ || \$5 ~ /DUP/) print}' ${FILEPATH}union_171212.sort.head.vcf > ${FILEPATH}union_171212_gt49bp.sort.vcf`;  #SVs >49bp
print `awk '{if(\$1 ~ /^#/ || length(\$4)-length(\$5)>19 || length(\$4)-length(\$5)<-19  || (length(\$4)>19 && length(\$5)>19)) print}' ${FILEPATH}union_171212.sort.head.vcf > ${FILEPATH}union_171212_refalt.sort.vcf`;  #SVs >49bp
print `awk '{if(\$1 ~ /^#/ || length(\$4)-length(\$5)>49 || length(\$4)-length(\$5)<-49 || (length(\$4)>49 && length(\$5)>49) ) print}' ${FILEPATH}union_171212_refalt.sort.vcf > ${FILEPATH}union_171212_refalt_gt49bp.sort.vcf`;  #SVs >49bp
print `awk 'BEGIN {FS=OFS="\t"} {endref=(length(\$4)+\$2) ; if(!(\$8 ~ "SVTYPE")) \$8=\$8";SVTYPE=INS;END="endref; print}' ${FILEPATH}union_171212_refalt.sort.vcf  | awk 'BEGIN {FS=OFS="\t"} {endref=(length(\$4)+\$2) ; if(!(\$8 ~ "END=")) \$8=\$8";END="endref; print}' | sed 's/SVTYPE=OTHER/SVTYPE=INS/;s/SVTYPE=DUP/SVTYPE=INS/;s/SVTYPE=COMPLEX/SVTYPE=INS/;s/SVTYPE=INS_EXTRA/SVTYPE=INS/;s/SVTYPE=DEL_EXTRA/SVTYPE=DEL/' | awk 'BEGIN {FS=OFS="\t"} {ref1=substr(\$4,0,1) ; if((\$5 == ".")) \$5=ref1; print}' | awk 'BEGIN {FS=OFS="\t"} {endref=(length(\$4)+\$2) ; if((\$8 ~ "SVTYPE=DEL" && length(\$4)>1 && length(\$5)>1)) \$8="SVTYPE=INS;END="endref; print}' > ${FILEPATH}union_171212_refalt.sort.formattedforsvviz.vcf`;


#create deletion bed file for input to R
print `${BEDTOOLSPATH}bin/multiIntersectBed -i $delbeds  -header -names $delbednames | sed 's/chrom/#chrom/' > ${FILEPATH}uniondel_171212.bed`;

print `${BEDTOOLSPATH}bin/mergeBed -i ${FILEPATH}uniondel_171212.bed -d 1000 -c 4,6,7,8,9,10 -o max > ${FILEPATH}uniondel_171212_merged1000.bed`;

print `cut -f1-3 ${FILEPATH}uniondel_171212_merged1000.bed > ${FILEPATH}uniondel_171212_merged1000_coordonly.bed`;

print `wc -l ${FILEPATH}uniondel_171212_merged1000.bed > ${FILEPATH}union_171212_linecounts.txt`;

print `awk '\$4>1' ${FILEPATH}uniondel_171212_merged1000.bed | wc -l >> ${FILEPATH}union_171212_linecounts.txt`;

print `${BEDTOOLSPATH}bin/annotateBed -i ${FILEPATH}uniondel_171212_merged1000_coordonly.bed -both -files $delbeds /Users/jzook/Downloads/opt/ga4gh-benchmarking-tools/resources/stratification-bed-files/LowComplexity/AllRepeats_gt95percidentity_slop5.bed.gz /Users/jzook/Downloads/opt/ga4gh-benchmarking-tools/resources/stratification-bed-files/SegmentalDuplications/hg19_self_chain_split_withalts_gt10k.bed.gz /Applications/bioinfo/nist-integration-v3.2.2/resources/example_of_no_ref_regions_input_file_b37.bed -names $delbednames tandemrep segdup refNs  | sed 's/^chr//;s/X/23/;s/Y/24/;s/^MT/25/;s/.*://' |  sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/;s/^25/MT/' | sed 's/^#./#/' > ${FILEPATH}uniondel_171212_callset.bed`;

my $col=($m-1)*2 + 5;
print `awk '{if (\$$col>0.25) sum+=1} END {print sum}' ${FILEPATH}uniondel_171212_callset.bed >> ${FILEPATH}union_171212_linecounts.txt`;
$col=$col+2;
print `awk '{if (\$$col>0.25) sum+=1} END {print sum}' ${FILEPATH}uniondel_171212_callset.bed >> ${FILEPATH}union_171212_linecounts.txt`;
$col=$col+2;
print `awk '{if (\$$col>0.01) sum+=1} END {print sum}' ${FILEPATH}uniondel_171212_callset.bed >> ${FILEPATH}union_171212_linecounts.txt`;
$col=$col-4;
print `awk '{if (\$$col>0.25) print }' ${FILEPATH}uniondel_171212_callset.bed > ${FILEPATH}uniondel_171212_callset_TRs.bed`;
print `awk '{if (\$$col<=0.25 || \$1 ~ /^#/) print }' ${FILEPATH}uniondel_171212_callset.bed > ${FILEPATH}uniondel_171212_callset_noTRs.bed`;

#create bed file for input to svrefine
print `cat ${FILEPATH}intermediate_171212/\*/\*_refalt.bed ${FILEPATH}intermediate_171212/\*/\*_imprecise.bed | sed 's/^chr//;s/X/23/;s/Y/24/;s/^MT/25/;s/.*://' | grep -v GL | grep -v hap | grep -v ctg | grep -v NC | grep -v hs | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/;s/^25/MT/' | ${BEDTOOLSPATH}bin/mergeBed -i stdin > ${FILEPATH}union_171212.merged.sort.bed`;
print `wc -l ${FILEPATH}union_171212.merged.sort.bed >> ${FILEPATH}union_171212_linecounts.txt`;
print `awk '{FS="\t";OFS="\t"} {print \$1,\$2,\$2+length(\$4)}' ${FILEPATH}union_171212.sort.head.vcf | grep -v '^#' > ${FILEPATH}union_171212.sort.head.bed`;
print `${BEDTOOLSPATH}bin/intersectBed -u -a ${FILEPATH}union_171212.merged.sort.bed -b ${FILEPATH}union_171212.sort.head.bed | wc -l >> ${FILEPATH}union_171212_linecounts.txt`;
print `awk '{FS="\t";OFS="\t"} {print \$1,\$2,\$2+length(\$4)}' ${FILEPATH}union_171212_refalt.sort.vcf  | grep -v '^#' > ${FILEPATH}union_171212_refalt.sort.bed`;
print `${BEDTOOLSPATH}bin/intersectBed -u -a ${FILEPATH}union_171212.merged.sort.bed -b ${FILEPATH}union_171212_refalt.sort.bed | wc -l >> ${FILEPATH}union_171212_linecounts.txt`;
print `grep -v refine ${FILEPATH}union_171212_refalt.sort.vcf > ${FILEPATH}union_171212_refalt.sort.norefine.vcf`;
print `awk '{FS="\t";OFS="\t"} {print \$1,\$2,\$2+length(\$4)}' ${FILEPATH}union_171212_refalt.sort.norefine.vcf  | grep -v '^#' > ${FILEPATH}union_171212_refalt.sort.norefine.bed`;
print `${BEDTOOLSPATH}bin/intersectBed -u -a ${FILEPATH}union_171212.merged.sort.bed -b ${FILEPATH}union_171212_refalt.sort.norefine.bed | wc -l >> ${FILEPATH}union_171212_linecounts.txt`;

print `grep -v '^#' ${FILEPATH}union_171212.sort.vcf | wc -l  >> ${FILEPATH}union_171212_linecounts.txt`;
print `grep -v '^#' ${FILEPATH}union_171212_refalt.sort.vcf | wc -l >> ${FILEPATH}union_171212_linecounts.txt`;
print `grep -v '^#' ${FILEPATH}union_171212_refalt.sort.vcf | cut -f1,2,4,5 | uniq | wc -l >> ${FILEPATH}union_171212_linecounts.txt`;
print `grep -v '^#' ${FILEPATH}union_171212_refalt_gt49bp.sort.vcf | wc -l >> ${FILEPATH}union_171212_linecounts.txt`;
print `grep -v '^#' ${FILEPATH}union_171212_refalt_gt49bp.sort.vcf | cut -f1,2,4,5 | uniq | wc -l >> ${FILEPATH}union_171212_linecounts.txt`;
print `grep -v '^#' ${FILEPATH}union_171212_refalt.sort.vcf | cut -f3 | sed 's/.*_\\(.*_.*\\)_.*/\\1/' | sort | uniq -c >> ${FILEPATH}union_171212_linecounts.txt`;


print `${BEDTOOLSPATH}bin/intersectBed -v -a ${FILEPATH}union_171212.merged.sort.bed -b ${FILEPATH}union_171212_refalt.sort.bed > ${FILEPATH}union_171212_norefalt.merged.sort.bed`;
print `${BEDTOOLSPATH}bin/intersectBed -u -a ${FILEPATH}union_171212.sort.head.vcf -b ${FILEPATH}union_171212_norefalt.merged.sort.bed > ${FILEPATH}union_171212_norefalt.sort.vcf`;

print `for v in ${FILEPATH}union_171212*.vcf; do /Applications/bioinfo/tabix/bgzip -f \$v; /Applications/bioinfo/tabix/tabix \$v.gz; done`;
print `for v in ${FILEPATH}intermediate_171212/*/*.vcf; do /Applications/bioinfo/tabix/bgzip -f \$v; /Applications/bioinfo/tabix/tabix \$v.gz; done`;
exit;

