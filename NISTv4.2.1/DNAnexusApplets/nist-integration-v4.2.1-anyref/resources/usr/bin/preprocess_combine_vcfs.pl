#To load following perl modules
use strict;
use warnings;
use Getopt::Long;

#Input parameters
my (@CALLSETTABLE, $FILEPATH);
my ($verbose);

#Declare variables
GetOptions(
	'cstable=s{1,}' => \@CALLSETTABLE,
	'p|filepath:s' => \$FILEPATH,
	"h|help|?" => \&usage);

#To check input files
my $CALLSETTABLE = $CALLSETTABLE[0];
if(defined($FILEPATH)){$FILEPATH = $FILEPATH} else {print usage();die "\n\nPlease specify the path of the vcf and bed files.\n\n"}

#To define usage of the preprocess_combine_vcfs
print "Usage = perl preprocess_combine_vcfs.pl -cstable $CALLSETTABLE ";
print "-p $FILEPATH \n\n";

sub usage {
	print "\nusage: perl preprocess_combine_vcfs.pl -cstable <txt> -p <FILEPATH> \n";
	print "\t-v \t\t\tThe input Tab-delimited text file containing columns Platform, Dataset, Callset, vcfAll.vcf.gz, callableBed.bed, and annotationsFile \n";
	print "\t-p, --filepath\t\t\t\tPath where vcf and bed files are located, ending in \\.\n";
	exit 1;
}

open (CALLSETTABLE, "<", $CALLSETTABLE) or die "Can't open the input table file. Please check if the file exists\n\n";

my $m = 1; #vcf number in order in callsettable
my $vcfs=""; #output variable that is printed that will contain names of processed vcfs for integration
my @line_head = (); #Callset table header line, which contains bed file names that may or may not be excluded from the callable regions of each callset

#loop through each line in callset table
while (my $line = <CALLSETTABLE>) {
	my $vcfno="$m";
	if ($m<10) {
		$vcfno="0$m";
	}
	chomp($line); #removes new line char
	
	my @line_data = split("\t", $line);
    
    if($line =~ m/^\#/) { #get header line and go to next line
        @line_head = split("\t", $line);
        next;
    }
    
    my $vcfstart=""; #begininning of vcf file name
    if($line_data[3] =~ /(.*)\.vcf.gz/) { $vcfstart=$1; } else { print "vcfAll does not end in .vcf.gz in cstable\n";exit 1; }
    
    #unzip vcf, keep first 10 columns, removes hom ref lines, removes chr prefix, breaks into "primitive" SNPs and indels, and sorts
    `gunzip -c ${FILEPATH}$line_data[3] | cut -f1-10 | grep -v '0/0' | grep -v '0|0' | vcfallelicprimitives -k -g | perl vcfsorter.pl genome.dict - > ${FILEPATH}${vcfstart}_nohomref.vcf`;
    
    #open formatted vcf from above
    open (VCFIN, "<", "${FILEPATH}${vcfstart}_nohomref.vcf") or die "Can't open the vcf ${FILEPATH}${vcfstart}_nohomref.vcf. Please check if the file exists\n\n";
	#open output vcf with additional formatting
	open (OUTINDIV, ">", "${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf") or die $!;
	#open output bed that will be used to remove regions with duplicate vcf lines (mostly due to a bug in freebayes parallel that outputs multiple vcf lines at edge of parallelized regions
	open (OUTDUPBED, ">", "${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_dup.bed") or die $!;
    
    my $lineprev=<VCFIN>;
    chomp($lineprev);
    my @line_data1prev = split("\t", $lineprev); #read in first line of vcf and put in @line_data1prev
    my $duppos=0;
    my $anydup=0;
    while (my $line1 = <VCFIN>) {
		chomp($line1);
		my @line_data1 = split("\t", $line1);
	
		if($lineprev =~ m/^\#CHROM/) { #change sample name so that they are ordered in the sequence of callsettable when using vcfcombine
			print OUTINDIV "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVCF${vcfno}_$line_data[2]\n";
			#print  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVCF${m}_$line_data[2]\n";
		} elsif ($lineprev =~ m/^\#/) {
			print OUTINDIV "$lineprev\n";
		} elsif ( $line_data1prev[1]==$line_data1[1] ) {
			$duppos=1; #duplicate position, so add line to bed file that excludes +-50bp around these vcf lines to avoid vcfcombine bug that ignores all remaining variants
			$anydup=1;
			my $bed1=$line_data1[1]-50;
			my $bed2=$line_data1[1]+50;
			print OUTDUPBED "$line_data1[0]\t$bed1\t$bed2\n";
			print "Duplicated vcf lines: $lineprev\n$line1\n$line_data1[0]\t$bed1\t$bed2\n";
		} elsif ($duppos==1) { #don't print the previous line if it was a duplicate
			$duppos=0;
		} else { #otherwise, print the previous line
			print OUTINDIV "$lineprev\n";
		}
		$lineprev=$line1;
		@line_data1prev = split("\t", $line1);
	}
	print OUTINDIV "$lineprev\n";
	close VCFIN;
	close OUTINDIV;
	close OUTDUPBED;
	
	#bgzip and index vcf output above
	`bgzip -c ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf > ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf.gz`;
	`tabix -p vcf ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf.gz`;
  	$vcfs="$vcfs ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf.gz";
	
	if (!($line_data[4] eq "none")) { #skip if no callable regions
		#Find lines with CALLABLE in column 4 (if callableloci bed file) or only 3 columns in bed file (assume all are callable in this case)
		if (`awk '{ print NF; exit }' ${FILEPATH}$line_data[4]` > 3) {
			if ($anydup==1) {
				`grep CALLABLE ${FILEPATH}$line_data[4] | cut -f1-3 | subtractBed -a stdin -b ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_dup.bed > ${FILEPATH}${vcfstart}_callable_tmp_1.bed`;
			} else {
				`grep CALLABLE ${FILEPATH}$line_data[4] > ${FILEPATH}${vcfstart}_callable_tmp_1.bed`;
			}
		} else {
			if ($anydup==1) {
				`cut -f1-3 ${FILEPATH}$line_data[4] | subtractBed -a stdin -b ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_dup.bed > ${FILEPATH}${vcfstart}_callable_tmp_1.bed`;
			} else {
				`mv ${FILEPATH}$line_data[4] ${FILEPATH}${vcfstart}_callable_tmp_1.bed`;
			}
		}
	
		#in v3.3, now subtract general difficult bed files from each callset's callable bed only when specified in callsettable
		my $c=6; #start with 7th column, which is the first bed file to exclude
		my $d=1;
		while ($c<=$#line_head) {
			if ($line_data[$c]==1) {
				`subtractBed -a ${FILEPATH}${vcfstart}_callable_tmp_${d}.bed -b ${line_head[$c]} > ${FILEPATH}${vcfstart}_callable_tmp_${c}.bed`;
				$d=$c;
			}
			$c++;
		}
		`awk 'BEGIN {FS = OFS = "\t"} {print \$1 "\t" \$2 "\t" \$3 "\t" "CS_$line_data[2]_callable"}' ${FILEPATH}${vcfstart}_callable_tmp_${d}.bed > ${FILEPATH}${vcfstart}_callable_processed.bed`;
		`rm ${FILEPATH}${vcfstart}_callable_tmp*.bed `;
	}
	$m=$m+1;

}
print "$vcfs\n";
#`vcfs="$vcfs"`;  #set environment variable that will contain names of processed vcfs for integration
close CALLSETTABLE;

#Add column to general difficult bed files to add annotations to union vcf
my $c=6; #start with 7th column, which is the first bed file to exclude
while ($c<=$#line_head) {
	my $bedname=$line_head[$c];
	if ($line_head[$c] =~ /(.*).bed.gz/) { 
		$bedname=$1; 
		`gunzip -c ${line_head[$c]} | awk 'BEGIN {FS = OFS = "\t"} {print \$1 "\t" \$2 "\t" \$3 "\t" "$bedname"}' > ${bedname}_diffbedannotcol.bed`;
	} elsif ($line_head[$c] =~ /(.*).bed/) { 
		$bedname=$1; 
		`awk 'BEGIN {FS = OFS = "\t"} {print \$1 "\t" \$2 "\t" \$3 "\t" "$bedname"}' ${line_head[$c]} > ${bedname}_diffbedannotcol.bed`;
	} 
	$c++;
}


 
