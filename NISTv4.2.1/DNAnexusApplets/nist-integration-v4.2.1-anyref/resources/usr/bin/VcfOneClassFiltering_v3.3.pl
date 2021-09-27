#!/usr/bin/perl -w
#
# VcfClassifyUsingFilters_v3.pl - arbitrate uncertain genotype calls using filters and callable regions, determined differently for each technology
# v3.1 - separate training for SNPs and non-SNPs
# v3.1.1 - allow different annotations for SNPs and non-SNPs

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use FileHandle;
use English;
use File::Basename;


my $line;


#create header
my $header = "";
# Check if the right number of annotations are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 4) {
    print "$#ARGV arguments\n usage: VcfOneClassFiltering_v3.3.pl TrainingVcf SensitiveVcf AnnotationsToFilter outputfilestart callsetname\n";
    print "example: perl VcfOneClassFiltering_v3.3.pl 20_vcfBeta-GS000025639-ASM_Filter_norm_sort_intersect2platforms.vcf 20_vcfBeta-GS000025639-ASM_Filter_norm_sort.vcf CG_Annotations.txt  20_vcfBeta-GS000025639-ASM_Filter_norm_sort_oneclassfilter CGnormal\n";
    print "Outputs are a bed file including regions of 50bp on either side of outlier or previously filtered vcf rows, as well as a new vcf with the FILTER field set to outlier for filtered rows\n"; 
    print "callsetname is the name of this callset from the callsettable input to VcfClassifyUsingFilters.pl\n"; 
    print "ParametersToFilter is a tab-delimited file with 3 columns and the following header describing the columns:\n"; 
    print "Annotation\t'Tails (-1 for left, 1 for right, 0 for both)'\t'0 for INFO or 1 for FORMAT'\t'SNPs'\t'indels''\t'BED file names\n"; exit;
}



my $trainvcf = "$ARGV[0]";
my $sensvcf = "$ARGV[1]";
my $annotations = "$ARGV[2]";
my $outfilestart = "$ARGV[3]";
my $callsetname = "$ARGV[4]";


unless ( open(VCFTRAIN, "${trainvcf}")) {
	print "\nCannot open the file: ${trainvcf}! \n\n";
	exit;
}
unless ( open(VCFSENS, "${sensvcf}")) {
	print "\nCannot open the file: ${sensvcf}! \n\n";
	exit;
}
unless ( open(ANNOTATIONS, "${annotations}")) {
	print "\nCannot open the file: ${annotations}! \n\n";
	exit;
}

# Output files are opened
unless ( open(OUTPUT, ">${outfilestart}_filtered.bed")) {
    print "\nCannot open the file: ${outfilestart}_filtered.bed to write to! \n\n";
    exit;
}
unless ( open(OUTVCF, ">${outfilestart}_filtered.vcf")) {
    print "\nCannot open the file: ${outfilestart}_filtered.vcf to write to! \n\n";
    exit;
}


#Read in information about each callset
my @annots=(("") x (20));
my @tails=(("") x (20));
my @infoformat=(("") x (20));
my @annsnp=((1) x (20));
my @annindel=((1) x (20));
my @fields;
$line=<ANNOTATIONS>; #skip header line
my $annotno=0;

while($line=<ANNOTATIONS>) {
	#print "$line\n";
	# Split up the line into an array
	@fields = split( "\t", $line);
	$annots[$annotno]=$fields[0];
	$tails[$annotno]=$fields[1];
	$infoformat[$annotno]=$fields[2];
	$annsnp[$annotno]=$fields[3];
	if ($fields[4] =~ /(.*)\n/) {$annindel[$annotno]=$1;}
	else { $annindel[$annotno]=$fields[4];}
	$annotno++;
}

#skip header lines in the training vcf 
    my $x=0;
	my $vgraph=0;
    while ($line = <VCFTRAIN>) {
        if ($line =~ /ID=BT,/) { $vgraph=1;} #check if BT format field is in the header - if so, then the vcf is from vgraph; if not, then use all rows 
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
        if ($x==1) {last;}
    }



my $ln_vcf_file = `awk \'END { print NR }\' $trainvcf `;

my $rand_row_number = 40000/int($ln_vcf_file);

print "Proportion of rows used to get 40000: $rand_row_number\n";

#Initialize 2D array to store annotations from SNPs and indels	
my @annotmatrixsnp;
@annotmatrixsnp=([(1000000.0) x 40000]);
for (my $j=1;$j<$annotno;$j++) {push @annotmatrixsnp, [(1000000.0) x 40000];}
#print "annotno:$annotno length:$#annotmatrixsnp\n";

my @annotmatrixindel;
@annotmatrixindel=([(1000000.0) x 40000]);
for (my $j=1;$j<$annotno;$j++) {push @annotmatrixindel, [(1000000.0) x 40000];}

#print $vgraph; exit;
my $rownosnp=0;
my $rownoindel=0;
my $skiprows=0;
while($line=<VCFTRAIN>) {
             
	# Split up the line into an array
	@fields = split( "\t", $line);
	
	#If length of REF or any ALT is greater than one, then treat as not a snp (i.e., an indel)
	my $nonsnp=0;
	my $ref = $fields[3];
	my $alt = $fields[4];
	my @alts = split( ",", $alt);
	if (length($ref)>1) { $nonsnp=1;} 
	if ($nonsnp==0) {
		for (my $j=0;$j<=$#alts;$j++) {
			if (length($alts[$j])>1) { $nonsnp=1; last; }
		}
	}
	
	#randomly select 40000 rows for SNPs
	if (($nonsnp==0 && ($rownosnp>=40000 || rand()>$rand_row_number)) || ($nonsnp==1 && ($rownoindel>=40000))) {
		$skiprows++;
		next;	
	}

	my $info = $fields[7];
	my $format = $fields[8];
	my @formats = split( ":", $format);
	my @chars = split( ":", $fields[9]);

	#First, check if # of technologies is >1 if the input is from the vgraph multimerge
	if ($vgraph==1) {
		my $techno=-1;
		for (my $k=0; $k<=$#chars; $k++) {
			if ($formats[$k] eq "BT") {
				my @charfld1=split( ",", $chars[$k]); #only use first value if field is multivalued (e.g., multiallelic GATK sites)
				if ($charfld1[0] =~ /(.*)\n/) {$techno=$1;}
				else { $techno=$charfld1[0];}
			}
		}
		if ($techno==0 || $techno==1 || $techno==-1) { $skiprows++; next;} #skip if # of techs from vgraph is 0 or 1 or if BT is missing (since filtered rows have no BT)
	}
	
	#Add annotations from this line to the matrix of annotations
	my $annotval=-1;
	for (my $j=0;$j<$annotno;$j++) {
		if ($infoformat[$j]==0) { #field is in INFO
			my @infoflds=split( ";", $info);
			for (my $k=0;$k<=$#infoflds;$k++) {
				my @infofld=split( "=", $infoflds[$k]);
				if ($#infofld==1 && $infofld[0] eq $annots[$j]) { 
					my @infofld1=split( ",", $infofld[1]); # use minimum of first 3 values if field is multivalued (e.g., multiallelic GATK sites)
					$annotval=$infofld1[0];
					if($#infofld1>0 && $annotval>$infofld1[1]) { $annotval=$infofld1[1];}
					if($#infofld1>1 && $annotval>$infofld1[2]) { $annotval=$infofld1[2];}
				}
			}
		} else { #field is in FORMAT
		
			for (my $k=0; $k<=$#chars; $k++) {
				if ($formats[$k] eq $annots[$j]) {
					my @charfld1=split( ",", $chars[$k]); #use minimum of first 3 values if field is multivalued (e.g., CGA_CEHQ or multiallelic GATK sites)
					if ($charfld1[0] =~ /(.*)\n/) {$annotval=$1;}
					else { 
						$annotval=$charfld1[0];
						if($#charfld1>0) {
							if ($charfld1[1] =~ /(.*)\n/) {if($annotval>$1) {$annotval=$1;}}
							else { 
								if ($annotval>$charfld1[1]) {$annotval=$charfld1[1];}
							}
							if($#charfld1>1) {
								if ($charfld1[2] =~ /(.*)\n/) {if($annotval>$1) {$annotval=$1;}}
								else { 
									if ($annotval>$charfld1[2]) {$annotval=$charfld1[2];}
								}
							}

						 }
					}
					#print "insidetest:$annotmatrixsnp[$j][$rowno];$rowno\n";
				}
			}
		}
		#Put annotation value in SNP or indel matrix
		if ($nonsnp==0) {
			$annotmatrixsnp[$j][$rownosnp]=$annotval;
		} else {
			$annotmatrixindel[$j][$rownoindel]=$annotval;
		}
	}
	if ($nonsnp==0) {
		$rownosnp++;
	} else {
		$rownoindel++;
	}
#	if ($rownoindel==10) {print "insidetest:$annotmatrixsnp[1][$rownosnp]\n"; exit;}
	
	#exit;
	
}
close VCFTRAIN;

#Step 2: For each annotation, sort the values from smallest to largest and find the appropriate cut-off(s)
my @annotmatrixsort;
print "$skiprows $rownosnp $rownoindel\n";

my @cutoffleftsnp = ((-1000000.0) x $annotno);
my @cutoffrightsnp = ((1000000.0) x $annotno);
my @cutoffleftindel = ((-1000000.0) x $annotno);
my @cutoffrightindel = ((1000000.0) x $annotno);

for (my $j=0;$j<$annotno;$j++) {
	#First, find cut-offs for SNPs
	my @annotmatrix1 = ((1000000.0) x $rownosnp); 
	my $rowno2=0;
	while ($rowno2<$rownosnp) { #first make new 1D array with column of interest
		$annotmatrix1[$rowno2]=$annotmatrixsnp[$j][$rowno2];
		$rowno2++;
	}
	@annotmatrixsort = sort { $a <=> $b } @annotmatrix1; #sort
	
	#Finding the number of rows that had a value for this annotation
	$rowno2=0;
	while ($rowno2<$rownosnp) {
		if ($annotmatrixsort[$rowno2] > 999999) { last; }
		$rowno2++;
	}
	#Find cut-off values for this annotation (Exclude extreme values less than or greater than 5%/(number of annotations)/(number of tails excluded))
	if ($tails[$j]==0) { #filter both tails
		$cutoffleftsnp[$j] = $annotmatrixsort[int(0.5+$rowno2*(0.05/$annotno/2))]; 
		$cutoffrightsnp[$j] = $annotmatrixsort[int(0.5+$rowno2*(1.0-0.05/$annotno/2))]; 
	} elsif ($tails[$j]==-1) { #filter left tail only
		$cutoffleftsnp[$j] = $annotmatrixsort[int(0.5+$rowno2*(0.05/$annotno))]; 
	} elsif ($tails[$j]==1) { #filter right tail only
		$cutoffrightsnp[$j] = $annotmatrixsort[int(0.5+$rowno2*(1.0-0.05/$annotno))]; 
	}		
	print "Number of records with values and left and right cutoffs for $annots[$j] for SNPs: $rowno2, $cutoffleftsnp[$j], $cutoffrightsnp[$j]\n";
	
	#Now, find cut-offs for indels
	@annotmatrix1 = ((1000000.0) x $rownoindel); 
	$rowno2=0;
	while ($rowno2<$rownoindel) { #first make new 1D array with column of interest
		$annotmatrix1[$rowno2]=$annotmatrixindel[$j][$rowno2];
		$rowno2++;
	}
	@annotmatrixsort = sort { $a <=> $b } @annotmatrix1; #sort
	$rowno2=0;
	while ($rowno2<$rownoindel) {
		#print OUTPUT "$j,$rowno2,$annotmatrixsort[$rowno2],$annotmatrixsort[1]\n";
		if ($annotmatrixsort[$rowno2] > 999999) { last; }
		$rowno2++;
	}
	#Only find cut-offs if there are >100 values to use
	if ($tails[$j]==0 && $rowno2>100) { #filter both tails
		$cutoffleftindel[$j] = $annotmatrixsort[int(0.499+$rowno2*(0.05/$annotno/2))-1]; 
		$cutoffrightindel[$j] = $annotmatrixsort[int(0.499+$rowno2*(1.0-0.05/$annotno/2))-1]; 
	} elsif ($tails[$j]==-1 && $rowno2>100) { #filter left tail only
		$cutoffleftindel[$j] = $annotmatrixsort[int(0.499+$rowno2*(0.05/$annotno))-1]; 
	} elsif ($tails[$j]==1 && $rowno2>100) { #filter right tail only
		$cutoffrightindel[$j] = $annotmatrixsort[int(0.499+$rowno2*(1.0-(0.05/$annotno)))-1]; 
		#print $rowno2*(1.0-(0.05/$annotno));
		#print "\n";
	}		
	print "Number of records with values and left and right cutoffs for $annots[$j] for indels: $rowno2, $cutoffleftindel[$j], $cutoffrightindel[$j]\n";
}



##finally, create a bed file using the vcf records that fall outside the cut-offs or are already filtered

#skip header lines in vcf file containing all sites
    $x=0;
    my $addfilt=1;
    while ($line = <VCFSENS>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
        if ($addfilt==1 && $line =~ /^\#\#INFO/) { #add new filter in header right before INFO header lines
        	print OUTVCF "##FILTER=<ID=outlier,Description=\"Filtered because one or more annotations were determined to be outliers using VcfOneClassFiltering_v3.1.pl\">\n";
        	$addfilt=0;
        }
        print OUTVCF $line;
        if ($x==1) {last;}
    }

my $bedchrom="";
my $bedstart=-1;
my $bedend=-1;
my $filtcnt=0;
my $totalcnt=0;
my @filtindivcntsnp=((0) x $annotno);
my @filtindivcntindel=((0) x $annotno);

while($line=<VCFSENS>) {
             
	# Split up the line into an array
	@fields = split( "\t", $line);
	
	
	my $nonsnp=0;
	my $nogeno=0; #don't output records with no genotype or partial genotypes
	my $ref = $fields[3];
	my $alt = $fields[4];
	my @alts = split( ",", $alt);
	
	#is this a SNP or not
	if (length($ref)>1) { $nonsnp=1;} 
	if ($nonsnp==0) {
		for (my $j=0;$j<=$#alts;$j++) {
			if (length($alts[$j])>1) { $nonsnp=1; last; }
		}
	}
	

	my $filt = $fields[6];
	my $info = $fields[7];
	my $format = $fields[8];
	my @formats = split( ":", $format);
	my @chars = split( ":", $fields[9]);

	my $filtline=0; #flag to set to 1 if line should be filtered
	if (!($filt eq "PASS" || $filt eq ".")) {
		$filtline=1;
	} else {
		for (my $j=0;$j<$annotno;$j++) {
			if ($infoformat[$j]==0) { #field is in INFO
				my @infoflds=split( ";", $info);
				for (my $k=0;$k<=$#infoflds;$k++) {
					my @infofld=split( "=", $infoflds[$k]);
					if ($#infofld==1 && $infofld[0] eq $annots[$j]) { 
						my @infofld1=split( ",", $infofld[1]); #only use first value if field is multivalued (e.g., multiallelic GATK sites)
						my $annotval=$infofld1[0];
						if($#infofld1>0 && $annotval>$infofld1[1]) { $annotval=$infofld1[1];}
						if($#infofld1>1 && $annotval>$infofld1[2]) { $annotval=$infofld1[2];}
						if ($nonsnp==0 && $annsnp[$j]==1 && ($annotval<$cutoffleftsnp[$j] || $annotval>$cutoffrightsnp[$j])) { $filtline=1; $filtindivcntsnp[$j]++;}
						if ($nonsnp==1 && $annindel[$j]==1 && ($annotval<$cutoffleftindel[$j] || $annotval>$cutoffrightindel[$j])) { $filtline=1; $filtindivcntindel[$j]++;}
					}
				}
			} else { #field is in FORMAT
	
				for (my $k=0; $k<=$#chars; $k++) {
					if ($formats[$k] eq $annots[$j] ) {
						my @charfld1=split( ",", $chars[$k]); #use minimum of first 3 values if field is multivalued (e.g., CGA_CEHQ or multiallelic GATK sites)
						my $charval=0;
						if ($charfld1[0] =~ /(.*)\n/) {$charval=$1;}
						else { 
							$charval=$charfld1[0];
							if($#charfld1>0) {
								if ($charfld1[1] =~ /(.*)\n/) {if($charval>$1) {$charval=$1;}}
								else { 
									if ($charval>$charfld1[1]) {$charval=$charfld1[1];}
								}
								if($#charfld1>1) {
									if ($charfld1[2] =~ /(.*)\n/) {if($charval>$1) {$charval=$1;}}
									else { 
										if ($charval>$charfld1[2]) {$charval=$charfld1[2];}
									}
								}

							 }
						}
						#print "insidetest:$annotmatrixsnp[$j][$rowno];$rowno\n";
						if ($nonsnp==0 && $annsnp[$j]==1 && (!($charval eq ".") && ($charval<$cutoffleftsnp[$j] || $charval>$cutoffrightsnp[$j]))) { $filtline=1; $filtindivcntsnp[$j]++;}
						if ($nonsnp==1 && $annindel[$j]==1 && (!($charval eq ".") && ($charval<$cutoffleftindel[$j] || $charval>$cutoffrightindel[$j]))) { $filtline=1; $filtindivcntindel[$j]++;}
					}
				}
			}
		}
		#also check if GQ<20
		for (my $k=0; $k<=$#chars; $k++) {
			if ($formats[$k] eq "GQ" ) {
				my $charval=0;
				my @charfld1=split( ",", $chars[$k]); #only use first value if field is multivalued (e.g., multiallelic GATK sites)
				if ($charfld1[0] =~ /(.*)\n/) {$charval=$1;}
				else { $charval=$charfld1[0];}
				#print "insidetest:$annotmatrixsnp[$j][$rowno];$rowno\n";
				if (!($charval eq ".") && ($charval<20)) { $filtline=1;}
			}
			if ($formats[$k] eq "GT" ) {
				if ($chars[$k] =~ /\./) {$nogeno=1;}
			}
		}
	}
	
	if ($filtline==1) {
		if ($bedend==-1) {
			$bedchrom = $fields[0];
			$bedstart = $fields[1]-1-50; #subtract 1 to create 0-based bed and add 50bp padding
			$bedend = $fields[1]+length($fields[3])+50; #add number of reference bases to uncertain region + 50bp padding
			
		} elsif ($fields[1]-1-50>$bedend || !($fields[0] eq $bedchrom)) { #print previous record if far from current row 
			print OUTPUT "$bedchrom\t$bedstart\t$bedend\tCS_${callsetname}_filt\n";
			$bedchrom = $fields[0];
			$bedstart = $fields[1]-1-50; #subtract 1 to create 0-based bed and add 50bp padding
			$bedend = $fields[1]+length($fields[3])+50; #add number of reference bases to uncertain region + 50bp padding
		} else { #extend previous record
			if ($bedend < $fields[1]+length($fields[3])+50) {$bedend=$fields[1]+length($fields[3])+50;} #add number of reference bases to uncertain region + 50bp padding
		}
		$filtcnt++;
		$filt="outlier";
	}
	
	if ($nogeno==0) {
		print OUTVCF "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$filt\t$fields[7]\t$fields[8]\t$fields[9]\n";
		$totalcnt++;
	}
}
print OUTPUT "$bedchrom\t$bedstart\t$bedend\tCS_${callsetname}_filt\n"; #print last record

print "Number of vcf records filtered: $filtcnt out of $totalcnt.\n";
for (my $j=0;$j<$annotno;$j++) {
	print "Number of snps filtered due to outlier $annots[$j]: $filtindivcntsnp[$j]\n";
	print "Number of non-snps filtered due to outlier $annots[$j]: $filtindivcntindel[$j]\n";

}

close VCFSENS;
close OUTPUT;
close OUTVCF;
    close ANNOTATIONS;

