#!/usr/bin/perl -w
#
# VcfClassifyUsingFilters_v3.pl - arbitrate uncertain genotype calls using filters and callable regions, determined differently for each technology
# v3.1 - if genotypes disagree between callsets and the arbitration decides to use 0/1, then make it uncertain, since CG and freebayes have this type of genotyping error
# v3.1 - remove callsettable column giving option to use only to confirm variants, since this can be done simply by using a notcallable bed containing the whole genome
#      - also filter implied homref records if another callset for the same dataset is filtered
# v3.1.1 - filter ion if it doesn't make a call at indels called by other methods, since it misses some
#        - don't report confident homref from any callset if it doesn't make a call at indels >10bp called by other methods, since none are perfectly sensitive to larger indels
# v3.3 - Changed not callable to callable bed file annotations
#      - added ADsum and alleleimbalance FILTER
#      - now take phasing and ID from GATK PGT and PID fields
# v3.3.1 - correct PS type to String in header
#        - enable GRCh38 or GRCh37
#        - filter sites where Ion misses indel even if all others are not callable
# v3.3.2 - filter sites that are within 50bp of another passing call but none of the callsets that support the 2 calls match, because some nearby conflicting calls from different callers were both considered high confidence if another callset from the same dataset was filtered
#        - transfer difficultregion annotation from input to output files
#		 - Fix interpretation of CG AD field for homozygous sites

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use FileHandle;
use English;
use File::Basename;

my $line;

#create header
my $header = "##fileformat=VCFv4.2\
##FILTER=<ID=GQlessthan70,Description=\"Sum of GQ for datasets with this genotype less than 70\">\
##FILTER=<ID=allfilteredanddisagree,Description=\"All callsets have this call filtered or outside the callable regions and they have discordant genotypes or variant calls\">\
##FILTER=<ID=allfilteredbutagree,Description=\"All callsets have this call filtered or outside the callable regions but they have the same genotype\">\
##FILTER=<ID=discordantunfiltered,Description=\"Callsets with unfiltered calls have discordant genotypes or variant calls\">\
##FILTER=<ID=discordanthet,Description=\"Filtered calls where a passing call is het and a high GQ but filtered call is hom var, since often the het is wrong\">\
##FILTER=<ID=questionableindel,Description=\"Filtered calls where some callsets have a filtered indel larger than 10bp and another dataset has an implied homozygous reference call\">\
##FILTER=<ID=cgonly,Description=\"Filtered calls where only Complete Genomics had this call and it was completely missing from any other callset\">\
##FILTER=<ID=alleleimbalance,Description=\"Filtered calls where the net allele balance for unfiltered datasets is <0.2 or >0.8\">\
##FILTER=<ID=overlappingcall,Description=\"Filtered sites that are within 50bp of another passing call but none of the callsets that support the 2 calls match\">\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth summed across all datasets, excluding MQ0 reads\">\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Net Genotype quality across all datasets, calculated from GQ scores of callsets supporting the consensus GT, using only one callset from each dataset\">\
##FORMAT=<ID=ADALL,Number=R,Type=Integer,Description=\"Net allele depths across all datasets\">\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Net allele depths across all unfiltered datasets with called genotype\">\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Consensus Genotype across all datasets with called genotype\">\
##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase set in which this variant falls\">\
##INFO=<ID=DPSum,Number=1,Type=Integer,Description=\"Total read depth summed across all datasets, excluding MQ0 reads\">\
##INFO=<ID=platforms,Number=1,Type=Integer,Description=\"Number of different platforms for which at least one callset called this genotype, whether filtered or not\">\
##INFO=<ID=platformnames,Number=.,Type=String,Description=\"Names of platforms for which at least one callset called this genotype, whether filtered or not\">\
##INFO=<ID=platformbias,Number=.,Type=String,Description=\"Names of platforms that have reads containing a variant at this location, but the high-confidence call is homozygous reference, indicating that there is a potential bias.\">\
##INFO=<ID=datasets,Number=1,Type=Integer,Description=\"Number of different datasets for which at least one callset called this genotype, whether filtered or not\">\
##INFO=<ID=datasetnames,Number=.,Type=String,Description=\"Names of datasets for which at least one callset called this genotype, whether filtered or not\">\
##INFO=<ID=datasetsmissingcall,Number=.,Type=String,Description=\"Names of datasets that are missing a call or have an incorrect call at this location, and the high-confidence call is a variant\">\
##INFO=<ID=callsets,Number=1,Type=Integer,Description=\"Number of different callsets that called this genotype, whether filtered or not\">\
##INFO=<ID=callsetnames,Number=.,Type=String,Description=\"Names of callsets that called this genotype, whether filtered or not\">\
##INFO=<ID=varType,Number=1,Type=String,Description=\"Type of variant\">\
##INFO=<ID=filt,Number=.,Type=String,Description=\"List of callsets that had this call filtered.\">\
##INFO=<ID=callable,Number=.,Type=String,Description=\"List of callsets that had this call in a region with low coverage of high MQ reads.\">\
##INFO=<ID=difficultregion,Number=.,Type=String,Description=\"List of difficult region bed files containing this call.\">\
##INFO=<ID=arbitrated,Number=1,Type=String,Description=\"TRUE if callsets had discordant calls so that arbitration was needed.\">\
##INFO=<ID=callsetwiththisuniqgenopassing,Number=.,Type=String,Description=\"Callset that uniquely calls the PASSing genotype in GT when 2+ PASSing callsets support a different genotype.\">\
##INFO=<ID=callsetwithotheruniqgenopassing,Number=.,Type=String,Description=\"Callset that uniquely calls a PASSing genotype different from GT when 2+ PASSing callsets support the genotype in GT.\">\n";

unless ( open(GENOME, "human.genome")) {
	print "\nCannot open the file: human.genome! \n\n";
	exit;
}
while(my $line=<GENOME>) {
	my @fields = split( "\t", $line);
	my $len1;
	if ($fields[1] =~ /(.*)\n/) {$len1=$1;}
	else { $len1=$fields[1];}
	$header = "${header}##contig=<ID=${fields[0]},length=${len1}>\n";
}
close GENOME;
$header = "${header}##fileDate=20160824\
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	INTEGRATION\n";


# Check if the right number of parameters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 2) {
    print "$#ARGV arguments\n usage: VcfClassifyUsingFilters.pl unionvcf CallsetTable outputfilestart\n";
    print "example: perl /Applications/bioinfo/perl/VcfClassifyUsingFilters.pl NA12878_union.vcf NA12878_RM8398_Datasets.txt NA12878_integration\n";
    print "CallsetTable must have callsets in the same order as the columns in the unionvcf\n"; 
    print "callsets should be in the order GATK, CG, freebayes, Ion TVC, SOLiD for best results.\n"; exit;
}



my $infile = "$ARGV[0]";
my $callsettable = "$ARGV[1]";
my $outfilestart = "$ARGV[2]";


unless ( open(VCFALL, "${infile}")) {
	print "\nCannot open the file: ${infile}! \n\n";
	exit;
}
unless ( open(CALLSETTABLEFILE, "${callsettable}")) {
	print "\nCannot open the file: ${callsettable}! \n\n";
	exit;
}

# Output files are opened
unless ( open(OUTPUT, ">${outfilestart}_ClassifyUsingFilters_allcalls.vcf")) {
    print "\nCannot open the file: ${outfilestart}_ClassifyUsingFilters_allcalls.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTPUTARBITRATED, ">${outfilestart}_ClassifyUsingFilters_arbitrated.vcf")) {
	print "\nCannot open the file: ${outfilestart}_ClassifyUsingFilters_arbitrated.vcf to write to! \n\n";
	exit;
}
unless ( open(OUTPUT2PLATFORMS, ">${outfilestart}_ClassifyUsingFilters_2platforms.vcf")) {
	print "\nCannot open the file: ${outfilestart}_ClassifyUsingFilters_2platforms.vcf to write to! \n\n";
	exit;
}
unless ( open(TESTOUT, ">${outfilestart}_testout.txt")) {
	print "\nCannot open the file: ${outfilestart}_testout.txt to write to! \n\n";
	exit;
}


#Read in information about each callset from callsettable
my @platform=(("") x (20));
my @dataset=(("") x (20));
my @callset=(("") x (20));
my @fields;
$line=<CALLSETTABLEFILE>; #skip header line
my $callsets=0;
while($line=<CALLSETTABLEFILE>) {
	#print "$line\n";
	# Split up the line into an array
	@fields = split( "\t", $line);
	$platform[$callsets]=$fields[0];
	$dataset[$callsets]=$fields[1];
	$callset[$callsets]=$fields[2];
	$callsets++;
}
#print "@platform\n@dataset\n@callset\n$callsets\n";
#skip header lines in individual files and write first one to output
    my $x=0;
    while ($line = <VCFALL>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
        if ($x==1) {last;}
    }

#print header created above to each output
print OUTPUT $header;
		print OUTPUTARBITRATED $header;
		print OUTPUT2PLATFORMS $header;

my $agreecallsetsprev="";
my $prevpasspos=-51;

#Loop through union vcf and determine which calls are high confidence
while($line=<VCFALL>) {
             
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
	my @callable=((0) x ($callsets));
	my @filt=((0) x ($callsets));
	my $callabletxt="";
	my $filttxt="";
	my @gt=(("") x ($callsets));
	my @gq=((-1) x ($callsets));
	my @gl=(("") x ($callsets));
	my @dp=((0) x ($callsets));
	my @ad1=((0) x ($callsets));
	my @ad2=((0) x ($callsets));
	my @ad3=((0) x ($callsets));
	my $gtphased="";
	my $ps=".";

	my $datasetfilttxt="";
	my $difficultregion="";
	
	my @infosplit = split( ";", $info);
	for (my $k=0; $k<=$#infosplit; $k++) {
		if ($infosplit[$k] =~ /difficultregion/) {
			$difficultregion=";$infosplit[$k]";
		} 
	} 

	#Loop through callsets to determine if callable or filtered, and then find FORMAT annotations for each callset like GT, GQ, etc.
	for (my $j=0;$j<$callsets;$j++) {
		#is this call in the callset's callable regions?
		if ($info =~ /CS_${callset[$j]}_callable/) {
			$callable[$j]=1;
			if ($callabletxt eq "") {
				$callabletxt=";callable=CS_${callset[$j]}_callable";
			} else {
				$callabletxt="$callabletxt,CS_${callset[$j]}_callable";
			}
		}
		#is this call in the callset's filtered regions?
		if ($info =~ /CS_${callset[$j]}_filt/) {
			$filt[$j]=1;
			if ($filttxt eq "") {
				$filttxt=";filt=CS_${callset[$j]}_filt";
			} else {
				$filttxt="$filttxt,CS_${callset[$j]}_filt";
			}
			if ($datasetfilttxt eq "") {
				$datasetfilttxt=";datasetfilt=DS_${dataset[$j]}_filt";
			} else {
				$datasetfilttxt="$datasetfilttxt,DS_${dataset[$j]}_filt";
			}
		}
	 
		#find GT and GQ for this callset
		my @chars = split( ":", $fields[$j+9]);
		if ($chars[0] =~ /\./) { #if no genotype, then don't need to find other fields
			$gt[$j]=".";
			$gq[$j]=60;
			$gl[$j]="";
			
		} else {
			$gt[$j]="";
			$gq[$j]=60;
			$gl[$j]="";
			my $revGT=0; #are GT's not in increasing numerical order; if not, reverse AD's
			for (my $k=0; $k<=$#chars; $k++) {
				#For GT field, format to (smaller number)/(larger number) to enable comparisons between callsets
				if ($formats[$k] eq "GT") {
					if ($chars[$k] =~ /(.*)\n/) {$gt[$j]=$1;}
					else { $gt[$j]=$chars[$k];}
					if ($gt[$j] =~ /(.*)\|(.*)/) { # if phased, convert from | to / for below analyses
						if ($j==0) {$gtphased=$gt[$j];} #only keep phasing information from first callset
						$gt[$j]="$1/$2";
					}
					if ($gt[$j] =~ /(\d)\/(\d)/) {if ($1>$2) {$gt[$j]="$2/$1"; $revGT=1;}} #order GT from smallest to largest number
					if ($gt[$j] =~ /\./) {$gt[$j]=".";}
					if ($gt[$j] eq "1") {$gt[$j]="1/1";}
				}
				#GATKHC outputs phased genotype in PGT
				if ($j==0 && $formats[$k] eq "PGT") {
					if ($gt[$j] eq "1/1") { 
						$gtphased="1|1"; #account for bug in GATKHC that outputs 0|1 in PGT at hom ref sites
					} else {
						if ($chars[$k] =~ /(.*)\n/) {$gtphased=$1;}
						else { $gtphased=$chars[$k];}
					}
				}
				#get phase set from first callset (normally in PS but in PID for GATKHC)
				if ($j==0 && ($formats[$k] eq "PS" || $formats[$k] eq "PID")) { 
					if ($chars[$k] =~ /(.*)\n/) {$ps=$1;}
					else { $ps=$chars[$k];}
				}

				#get genotype quality
				if ($formats[$k] eq "GQ") {
					if ($chars[$k] =~ /(.*)\n/) {$gq[$j]=$1;}
					else { $gq[$j]=$chars[$k];}
					if ($gq[$j] eq ".") { $gq[$j]=-1; }
				}
				
				#get coverage from DP
				if ($formats[$k] eq "DP") {
					if ($chars[$k] =~ /(.*)\n/) {$dp[$j]=$1;}
					else { $dp[$j]=$chars[$k];}
				}
				
				#get allele depths, flip the order if we flipped the order of the GT field above
				if ($formats[$k] eq "AD") {
					if ($chars[$k] =~ /.*,.*,.*/) { #if 2 ALT alleles, then 3 allele depths
						if ($chars[$k] =~ /(.*),(.*),(.*)\n/) {$ad1[$j]=$1;$ad2[$j]=$2;$ad3[$j]=$3;}
						elsif ($chars[$k] =~ /(.*),(.*),(.*)/) {$ad1[$j]=$1;$ad2[$j]=$2;$ad3[$j]=$3;}
						if ($revGT==1) {my $x=$ad3[$j]; $ad3[$j]=$ad2[$j]; $ad2[$j]=$x;}
					} else {
						if ($chars[$k] =~ /(.*),(.*)\n/) {$ad1[$j]=$1;$ad2[$j]=$2;}
						elsif ($chars[$k] =~ /(.*),(.*)/) {$ad1[$j]=$1;$ad2[$j]=$2;}
						if ($revGT==1) {my $x=$ad1[$j]; $ad1[$j]=$ad2[$j]; $ad2[$j]=$x;}
					}
				}
			}
		}
	}


	#first, check if all genotypes that are callable agree
	my $gt1=""; #genotype of first callset without low cov
	my $agree=0; #number of callsets that agree with gt1 regardless of filter
	my $disagree=0; #number of callsets that disagree with gt1 regardless of filter 
	my $agreecallsets=""; #list of callsets that agree with gt1 regardless of filter
	my $agreecallsetscnt=0; #number of callsets that agree with gt1 regardless of filter
	my $agreedatasets=""; #list of datasets that agree with gt1 regardless of filter
	my $agreedatasetscnt=0; #number of datasets that agree with gt1 regardless of filter
	my $agreeplatforms=""; #list of platforms that agree with gt1 regardless of filter
	my $agreeplatformscnt=0; #number of platforms that agree with gt1 regardless of filter
	my $notcallablecnt=0; #number of call sets that are not callable 
	my $gqsum=0; #sum GQ's of the datasets with the final genotype
	my $DPSum=0; #sum DP's of the datasets with the final genotype
	my $ADSum1=0; #sum AD's of the datasets with the final genotype and not filtered
	my $ADSum2=0; #sum AD's of the datasets with the final genotype and not filtered
	my $ADSum3=0; #sum AD's of the datasets with the final genotype and not filtered
	my $ADallSum1=0; #sum AD's of the datasets with the final genotype
	my $ADallSum2=0; #sum AD's of the datasets with the final genotype
	my $ADallSum3=0; #sum AD's of the datasets with the final genotype
	#Loop through all callsets, check if they are callable, and check if callable callsets have the same genotype
	for (my $j=0;$j<$callsets;$j++) {
		if ($gt1 eq "") { #Have we found any callable callsets yet?
			if ($callable[$j]==0 || (!($gt[$j] eq ".") && $gq[$j]<20)) { #is this not callable or (not homref and GQ<20)?
				$notcallablecnt++;
			} else { #otherwise, it's callable
				$gt1=$gt[$j];
				$agree++;
				$agreecallsets=$callset[$j];
				$agreecallsetscnt++;
				$agreedatasets=$dataset[$j];
				$agreedatasetscnt++;
				$agreeplatforms=$platform[$j];
				$agreeplatformscnt++;
				$gqsum=$gqsum+$gq[$j];
			}
		} else { #if already found first callable callset, then compare to first callable callset's genotype
			if ($callable[$j]==0 || (!($gt[$j] eq ".") && $gq[$j]<20)) { #is this not callable or (not homref and GQ<20)?
				$notcallablecnt++;
			} elsif ($gt[$j] eq $gt1) { #otherwise, it's callable, so does it agree with the first callable callset's GT
				$agree++;
				if ($gqsum==0 || !($agreedatasets =~ /$dataset[$j]/)) {$gqsum=$gqsum+$gq[$j];} #only add to gqsum once for each dataset
				if (!($agreedatasets =~ /$dataset[$j]/)) {
					$agreedatasets="$agreedatasets,$dataset[$j]";
					$agreedatasetscnt++;
				}
				if (!($agreecallsets =~ /$callset[$j]/)) {
					$agreecallsets="$agreecallsets,$callset[$j]";
					$agreecallsetscnt++;
				}
				if (!($agreeplatforms =~ /$platform[$j]/)) {
					$agreeplatforms="$agreeplatforms,$platform[$j]";
					$agreeplatformscnt++;
				}
			} else { #otherwise, it disagrees with the first GT
				$disagree++;
			}
		}
	}
	
	#Next, use filtering information in addition to callable information
	my $gt1nofilt=""; #genotype of first callset without low cov or filter
	my $agreenofilt=0; #number of callsets that agree with gt1nofilt and not filtered
	my $disagreenofilt=0; #number of callsets that disagree with gt1nofilt and not filtered 
	my $notcallablefiltcnt=0; #number of call sets with low cov or genotype with . or filter
	my $agreecallsetsnofilt=""; #list of datasets that agree with gt1nofilt and not filtered
	my $agreecallsetsnofiltcnt=0; #number of datasets that agree with gt1nofilt and not filtered
	my $agreedatasetsnofilt=""; #list of datasets that agree with gt1nofilt and not filtered
	my $agreedatasetsnofiltcnt=0; #number of datasets that agree with gt1nofilt and not filtered
	my $agreeplatformsnofilt=""; #list of platforms that agree with gt1nofilt and not filtered
	my $agreeplatformsnofiltcnt=0; #number of platforms that agree with gt1nofilt and not filtered
	my $filtout=""; #output text to FILTER field
	my $arbitratetxt=""; #output to INFO if arbitrated
	if ($agree>0 && $disagree==0) {
		#all callsets agree, so check to make sure they're not all filtered
		#$gqsum=0;
		for (my $j=0;$j<$callsets;$j++) {
			if ($gt[$j] eq $gt1 && $callable[$j]>0 && $filt[$j]==0 && ($gt[$j] eq "." || $gq[$j]>=20) && !(($gt[$j] eq ".") && ($datasetfilttxt =~ /DS_${dataset[$j]}_filt/))) { #is this callset's GT the same as the first callable callset's GT, is this one callable, is it not filtered, is GQ>=20, and (is it not (filtered in another callset for this dataset and homozygous reference))
				#TODO: Should probably reverse the logic here to match above and below, and also probably add the Ion indel filter from below
				$agreenofilt++;
				if ($agreenofilt==1) {
					$agreecallsetsnofilt=$callset[$j];
					$agreecallsetsnofiltcnt++;
					$agreedatasetsnofilt=$dataset[$j];
					$agreedatasetsnofiltcnt++;
					$agreeplatformsnofilt=$platform[$j];
					$agreeplatformsnofiltcnt++;
				} else {
					if (!($agreecallsetsnofilt =~ /$callset[$j]/)) {
						$agreecallsetsnofilt="$agreecallsetsnofilt,$callset[$j]";
						$agreecallsetsnofiltcnt++;
					}
					#if (!($agreedatasetsnofilt =~ /$dataset[$j]/)) {
					#	$agreedatasetsnofilt="$agreedatasetsnofilt,$dataset[$j]";
					#	$agreedatasetsnofiltcnt++;
					#}
					if (!($agreeplatformsnofilt =~ /$platform[$j]/)) {
						$agreeplatformsnofilt="$agreeplatformsnofilt,$platform[$j]";
						$agreeplatformsnofiltcnt++;
					}
					#$gqsum=$gqsum+$gq[$j];
				}
			}
		}
		if ($gt1 eq ".") {
			$gt1="0/0";
			$filtout=".";
		}
		if ($gt1 eq "0/0" && $agreenofilt>0 && ((length($ref)>10 || length($alt)>10) || ($agreeplatformsnofilt eq "Ion" && (length($ref)>1 || length($alt)>1)))) { 
			$filtout="questionableindel"; #even if an unfiltered callset calls homozygous reference, make it uncertain if another callset has a filtered indel >10bp (or any indel if Ion misses it) because we found some callsets miss larger indels
		} elsif ($gt1 eq "0/0" && $agreenofilt>0) {
			$filtout="."; #High-confidence homozygous reference if not a long indel
		} elsif ($agreenofilt>0 && $gqsum<70) {
			$filtout="GQlessthan70"; #if all callable, unfiltered datasets with this call sum to a GQ < 70, then make uncertain because there is insufficient support
		} elsif ($agreenofilt>0) {
			$filtout="PASS"; #High-confidence variant if at least 1 unfiltered, callable callset
		} else {
			$filtout="allfilteredbutagree"; #otherwise, no unfiltered, callable callsets, so uncertain position
		}
	} elsif ($agree>0) {
		#If all genotypes that were callable did not agree, then see if filtering information helps to arbitrate and determine which callset to trust
		
		$gqsum=0;
		for (my $j=0;$j<$callsets;$j++) {
			if ($gt1nofilt eq "") { #if haven't found first unfiltered call set yet
				if ($callable[$j]==0 || $filt[$j]>0 || (!($gt[$j] eq ".") && $gq[$j]<20) || (($gt[$j] eq ".") && ($datasetfilttxt =~ /DS_${dataset[$j]}_filt/)) || (($gt[$j] eq ".") && ($platform[$j] eq "Ion") && (length($ref)>1 || length($alt)>1))) { #is this callset untrustworthy here
					$notcallablefiltcnt++;
				} else { #Otherwise, we trust it and get its genotype
					$gt1nofilt=$gt[$j];
					$agreenofilt++;
					$agreecallsetsnofilt=$callset[$j];
					$agreecallsetsnofiltcnt++;
					$agreedatasetsnofilt=$dataset[$j];
					$agreedatasetsnofiltcnt++;
					$agreeplatformsnofilt=$platform[$j];
					$agreeplatformsnofiltcnt++;
					$gqsum=$gqsum+$gq[$j];
				}
			} else { #compare to first unfiltered genotype if it is callable and unfiltered
				if ($callable[$j]==0 || $filt[$j]>0 || (!($gt[$j] eq ".") && $gq[$j]<20) || (($gt[$j] eq ".") && ($datasetfilttxt =~ /DS_${dataset[$j]}_filt/)) || (($gt[$j] eq ".") && ($platform[$j] eq "Ion") && (length($ref)>1 || length($alt)>1))) {  #is this callset untrustworthy here
					$notcallablefiltcnt++;
				} elsif ($gt[$j] eq $gt1nofilt) { #Otherwise, we trust it and compare to first trusted callset's genotype
					$agreenofilt++;
					if ($gqsum==0 || !($agreedatasetsnofilt =~ /$dataset[$j]/)) {$gqsum=$gqsum+$gq[$j];} #only add to gqsum once for each dataset
					if (!($agreecallsetsnofilt =~ /$callset[$j]/)) {
						$agreecallsetsnofilt="$agreecallsetsnofilt,$callset[$j]";
						$agreecallsetsnofiltcnt++;
					}
					if (!($agreedatasetsnofilt =~ /$dataset[$j]/)) {
						$agreedatasetsnofilt="$agreedatasetsnofilt,$dataset[$j]";
						$agreedatasetsnofiltcnt++;
					}
					if (!($agreeplatformsnofilt =~ /$platform[$j]/)) {
						$agreeplatformsnofilt="$agreeplatformsnofilt,$platform[$j]";
						$agreeplatformsnofiltcnt++;
					}
					
				} else { #otherwise, this callset disagrees with the first callset
					$disagreenofilt++;
				}
			}
		}
		if ($gt1nofilt eq "." ) {
			$gt1="0/0";
		} elsif (!($gt1nofilt eq "")) {
			$gt1=$gt1nofilt;
		}
		if (!($gt1nofilt eq "")) {
			$agreecallsets=$agreecallsetsnofilt;
			$agreecallsetscnt=$agreecallsetsnofiltcnt;
			$agreedatasets=$agreedatasetsnofilt;
			$agreedatasetscnt=$agreedatasetsnofiltcnt;
			$agreeplatforms=$agreeplatformsnofilt;
			$agreeplatformscnt=$agreeplatformsnofiltcnt;
		}

		if ($gt1 eq "0/0" && $agreenofilt>0 && ((length($ref)>10 || length($alt)>10) || ($agreeplatformsnofilt eq "Ion" && (length($ref)>1 || length($alt)>1)))) { 
			$filtout="questionableindel"; #even if an unfiltered callset calls homozygous reference, make it uncertain if another callset has a filtered indel >10bp (or any indel if Ion misses it) because we found some callsets miss larger indels
		} elsif ($gt1 eq "0/0" && $agreenofilt>0) {
			$filtout="."; #High-confidence homozygous reference if not a long indel
		} elsif ($agreenofilt>0 && $gqsum<70) {
			$filtout="GQlessthan70"; #if all callable, unfiltered datasets with this call sum to a GQ < 70, then make uncertain because there is insufficient support
		} elsif ($agreenofilt>0) {
			$filtout="PASS"; #High-confidence variant if at least 1 unfiltered, callable callset
		} else {
			$filtout="allfilteredbutagree"; #otherwise, no unfiltered, callable callsets, so uncertain position
		}


		if ($gt1 eq "0/0" && $agreenofilt>0 && $disagreenofilt==0 && ((length($ref)>10 || length($alt)>10) || ($agreeplatformsnofilt eq "Ion" && (length($ref)>1 || length($alt)>1)))) { 
			$filtout="questionableindel"; #even if an unfiltered callset calls homozygous reference, make it uncertain if another callset has a filtered indel >10bp (or any indel if Ion misses it) because we found some callsets miss larger indels
		} elsif ($gt1 eq "0/0" && $agreenofilt>0 && $disagreenofilt==0) {
			$filtout="."; #High-confidence homozygous reference if not a long indel
			$arbitratetxt=";arbitrated=TRUE";
		} elsif ($agreenofilt>0 && $disagreenofilt==0 && $gqsum<70) {
			$filtout="GQlessthan70"; #if all callable, unfiltered datasets with this call sum to a GQ < 70, then make uncertain because there is insufficient support
		} elsif ($agreenofilt>0 && $disagreenofilt==0 && $gqsum>=70) {
			$filtout="PASS"; #High-confidence variant if at least 1 unfiltered, callable callset
			$arbitratetxt=";arbitrated=TRUE";
		} elsif ($agreenofilt==0) {
			$filtout="allfilteredanddisagree"; #no unfiltered, callable callsets, so uncertain position
		} else {
			$filtout="discordantunfiltered"; #otherwise, unfiltered callable callsets had discordant genotype calls
		}
		#if ($pos==299593) { print "$gt1nofilt;$gt1;$agreenofilt;$agree;$disagreenofilt\n";exit;}
	} else {
		#Otherwise, all callsets were not callable or had low GQ, so get the first genotype of a callset that is variant
		for (my $j=0;$j<$callsets;$j++) {
			if ($gt1 eq "") {
				if (!($gt[$j] eq ".")) {
					$gt1=$gt[$j];
					$agree++;
					$agreedatasets=$dataset[$j];
					$agreedatasetscnt++;
					$agreeplatforms=$platform[$j];
					$agreeplatformscnt++;
					$gqsum=$gqsum+$gq[$j];
				}
			} else {
				if ($gt[$j] eq $gt1) {
					$agree++;
					if ($gqsum==0 || !($agreedatasets =~ /$dataset[$j]/)) {$gqsum=$gqsum+$gq[$j];} #only add to gqsum once for each dataset
					if (!($agreecallsets =~ /$callset[$j]/)) {
						$agreecallsets="$agreecallsets,$callset[$j]";
						$agreecallsetscnt++;
					}
					if (!($agreedatasets =~ /$dataset[$j]/)) {
						$agreedatasets="$agreedatasets,$dataset[$j]";
						$agreedatasetscnt++;
					}
					if (!($agreeplatforms =~ /$platform[$j]/)) {
						$agreeplatforms="$agreeplatforms,$platform[$j]";
						$agreeplatformscnt++;
					}
				} else {
					$disagree++;
				}
			}
		}
		if ($disagree>0) {
			$filtout="allfilteredanddisagree";
		} else {
			$filtout="allfilteredbutagree";
		}
	}
	   
	#We have now determined if the site is high-confidence, and if not, why not, so sum DP and AD for output annotations
	#sum DP for datasets that agree and count datasets and platforms that disagree and aren't filtered
	$agreedatasetsnofilt=""; #list of datasets that agree with gt1nofilt and not filtered; reset so that it can be used below to determine if AD's have been added from this dataset yet
	$agreedatasetsnofiltcnt=0; #number of datasets that agree with gt1nofilt and not filtered
	my $disagreecallsetsnofilt="";
	my $disagreedatasetsnofilt="";
	my $disagreeplatformsnofilt="";
	my $disagreecallsetsnofiltcnt=0;
	my $disagreedatasetsnofiltcnt=0;
	my $disagreeplatformsnofiltcnt=0;
	
	for (my $j=0;$j<$callsets;$j++) {
		if ($gt[$j] eq $gt1 || ($gt[$j] eq "." && $gt1 eq "0/0")) { #Does this callset agree with the output genotype
			if ($dp[$j] =~ /\d/) { $DPSum+=$dp[$j]; } #add its DP to DPSum
			if ($gqsum==0 || !($agreedatasets =~ /$dataset[$j]/)) {$gqsum=$gqsum+$gq[$j];} #only add to gqsum once for each dataset
			if ($ADallSum1+$ADallSum2+$ADallSum3==0 || !($agreedatasets =~ /$dataset[$j]/)) {$ADallSum1=$ADallSum1+$ad1[$j]; $ADallSum2=$ADallSum2+$ad2[$j]; $ADallSum3=$ADallSum3+$ad3[$j];} #only add to AD sums once for each dataset
			if (($callable[$j]>0 && $filt[$j]==0) && ($ADSum1+$ADSum2+$ADSum3==0 || !($agreedatasetsnofilt =~ /$dataset[$j]/))) {
				$ADSum1=$ADSum1+$ad1[$j]; $ADSum2=$ADSum2+$ad2[$j]; $ADSum3=$ADSum3+$ad3[$j];
				if (!($agreedatasetsnofilt =~ /$dataset[$j]/)) { 
						$agreedatasetsnofilt="$agreedatasetsnofilt,$dataset[$j]";
						$agreedatasetsnofiltcnt++;
				}
			} 
			
			if (!($agreecallsets =~ /$callset[$j]/)) {
				$agreecallsets="$agreecallsets,$callset[$j]";
				$agreecallsetscnt++;
			}
			if (!($agreedatasets =~ /$dataset[$j]/)) {
				$agreedatasets="$agreedatasets,$dataset[$j]";
				$agreedatasetscnt++;
			}
			if (!($agreeplatforms =~ /$platform[$j]/)) {
				$agreeplatforms="$agreeplatforms,$platform[$j]";
				$agreeplatformscnt++;
			}
		} else { #otherwise, this callset's GT is different from the output genotype
			if ($agreeplatformsnofiltcnt==1 && $filtout eq "PASS" && $gt1 eq "0/1" && $gt[$j] eq "1/1" && $gq[$j]>70) { $filtout="discordanthet"; } #fix problem with some CG and freebayes sites incorrectly calling hets when they should be homvar
			if (!($disagreecallsetsnofilt =~ /$callset[$j]/)) {
				$disagreecallsetsnofilt="$disagreecallsetsnofilt,$callset[$j]";
				$disagreecallsetsnofiltcnt++;
			}
			if (!($disagreedatasetsnofilt =~ /$dataset[$j]/)) {
				$disagreedatasetsnofilt="$disagreedatasetsnofilt,$dataset[$j]";
				$disagreedatasetsnofiltcnt++;
			}
			if (!($disagreeplatformsnofilt =~ /$platform[$j]/)) {
				$disagreeplatformsnofilt="$disagreeplatformsnofilt,$platform[$j]";
				$disagreeplatformsnofiltcnt++;
			}
			
		}
	}
	
	if ($gt1 eq "0/1" && $filtout eq "PASS" && ($ADSum1>4*$ADSum2 || $ADSum2>4*$ADSum1)) {$filtout="alleleimbalance";} #filter sites that have allele balance < 0.2 or > 0.8
	if ($agreeplatformscnt==1 && ($agreeplatforms =~ /CG/) && ($filtout eq "PASS")) { $filtout="cgonly";}
	$gqsum = int($gqsum+0.5); #round gqsum
	my $ADSumout = "$ADSum1,$ADSum2"; my $ADallSumout = "$ADallSum1,$ADallSum2";
	if ($alt =~ /,/) {$ADSumout = "$ADSum1,$ADSum2,$ADSum3"; $ADallSumout = "$ADallSum1,$ADallSum2,$ADallSum3";}
	
	#filter sites that are within 50bp of another passing call but none of the callsets that support the 2 calls match
	if (!($gt1 eq "0/0") && $filtout eq "PASS" && $pos - $prevpasspos < 50 && $pos - $prevpasspos >= 0 ) {
		my @currentcallsets = split(",", $agreecallsets);
		my $matchcallset=0;
		for (my $j=0;$j<=$#currentcallsets;$j++) {
			if ($agreecallsetsprev =~ /$currentcallsets[$j]/) { $matchcallset=1; }
		}
		if ($matchcallset==0) { $filtout="overlappingcall"; } 
		
	} 
	if (!($gt1 eq "0/0") && $filtout eq "PASS") {
		$prevpasspos=$pos;
		$agreecallsetsprev=$agreecallsets;
	}
	
	#output if platform appears to have a bias causing reads with evidence for a false variant
	my $disagreetxt="";
	if ($gt1 eq "0/0" && $filtout eq ".") {
		$DPSum = ".";
		$ADSumout = ".";
		$ADallSumout = ".";
		$gqsum = ".";
		if ($disagreeplatformsnofilt =~ /,(..*)/) {
			$disagreetxt=";platformbias=$1"; #output platforms that have evidence for a variant at a high-confidence hom ref site, since these are potentially caused by bias
		}
	} elsif (!($gt1 eq "0/0") && $filtout eq "PASS") {
		if ($disagreedatasetsnofilt =~ /,(..*)/) {
			$disagreetxt=";datasetsmissingcall=$1"; #output platforms that have an incorrect call at a high-confidence variant site
		}
	} elsif ($filtout eq "discordantunfiltered" ) { #output if there is a callset that uniquely calls a PASSing genotype when 2+ PASSing callsets support a different genotype
		if ($agreecallsetsnofiltcnt>1 && $disagreecallsetsnofiltcnt==1 && $disagreecallsetsnofilt =~ /,(..*)/) {
			$disagreetxt=";callsetwithotheruniqgenopassing=$1";
		} elsif ($disagreecallsetsnofiltcnt>1 && $agreecallsetsnofiltcnt==1 && $agreecallsetsnofilt =~ /,(..*)/) {
			$disagreetxt=";callsetwiththisuniqgenopassing=$1";
		}
	}
	my $qualout=5; #default QUAL to 5 if filtered
	if ($filtout eq "PASS") {
		$qualout=50; #50 if passing variant
	} elsif ($filtout eq ".") {
		$qualout=0; #0 if passing homref site
	}
	if (!($gtphased eq "") && $gt1 eq $gt[0]) { $gt1=$gtphased; } #keep phasing info from first callset

	print OUTPUT "$chrom\t$pos\t$id\t$ref\t$alt\t$qualout\t$filtout\tplatforms=$agreeplatformscnt;platformnames=$agreeplatforms;datasets=$agreedatasetscnt;datasetnames=$agreedatasets;callsets=$agreecallsetscnt;callsetnames=$agreecallsets${disagreetxt}${callabletxt}${filttxt}${arbitratetxt}${difficultregion}\tGT:PS:DP:ADALL:AD:GQ\t$gt1:$ps:$DPSum:$ADallSumout:$ADSumout:$gqsum\n";
	if ($filtout eq "PASS") {
		print OUTPUTARBITRATED "$chrom\t$pos\t$id\t$ref\t$alt\t$qualout\t$filtout\tplatforms=$agreeplatformscnt;platformnames=$agreeplatforms;datasets=$agreedatasetscnt;datasetnames=$agreedatasets;callsets=$agreecallsetscnt;callsetnames=$agreecallsets${disagreetxt}${callabletxt}${filttxt}${arbitratetxt}${difficultregion}\tGT:PS:DP:ADALL:AD:GQ\t$gt1:$ps:$DPSum:$ADallSumout:$ADSumout:$gqsum\n";
	}
	if ($filtout eq "PASS" && $agreeplatformsnofiltcnt>1) {
		print OUTPUT2PLATFORMS "$chrom\t$pos\t$id\t$ref\t$alt\t$qualout\t$filtout\tplatforms=$agreeplatformscnt;platformnames=$agreeplatforms;datasets=$agreedatasetscnt;datasetnames=$agreedatasets;callsets=$agreecallsetscnt;callsetnames=$agreecallsets${disagreetxt}${callabletxt}${filttxt}${arbitratetxt}${difficultregion}\tGT:PS:DP:ADALL:AD:GQ\t$gt1:$ps:$DPSum:$ADallSumout:$ADSumout:$gqsum\n";
	}
	if ($agreeplatformsnofilt eq "CG" && $gt1 eq "1/1") {
		print TESTOUT "$DPSum\t$gqsum\n";
	}
	
}



close CALLSETTABLEFILE;
close VCFALL;
close OUTPUT;
    close OUTPUTARBITRATED;
    close OUTPUT2PLATFORMS;
close TESTOUT;

