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
if(defined($FILEPATH)){$FILEPATH = $FILEPATH} else {print usage(); die "\n\nPlease specify the path of the vcf and bed files.\n\n"}

#To define usage of the preprocess_combine_vcfs
print "Usage = perl RunOneClassFilter.pl -cstable $CALLSETTABLE ";
print "-p $FILEPATH \n\n";

sub usage {
	print "\nusage: perl RunOneClassFilter.pl -cstable <txt> -p <FILEPATH> \n";
	print "\t-v \t\t\tThe input Tab-delimited text file containing columns Platform, Dataset, Callset, vcfAll.vcf.gz, callableBed.bed, and annotationsFile \n";
	print "\t-p, --filepath\t\t\t\tPath where vcf and bed files are located, ending in \\.\n";
	exit 1;
}

open (CALLSETTABLE, "<", $CALLSETTABLE) or die "Can't open the input table file. Please check if the file exists\n\n";

my $m = 1;
my $vcfs=""; #set environment variable that will contain names of processed vcfs for integration
while (my $line = <CALLSETTABLE>) {
	my $vcfno="$m";
	if ($m<10) {
		$vcfno="0$m";
	}
	chomp($line);
	my @line_data = ();
	@line_data = split("\t", $line);
    
    if($line =~ m/^\#/) {
        next;
    }
    my $vcfstart="";
    if($line_data[3] =~ /(.*)\.vcf.gz/) { $vcfstart=$1; } else { print "vcfAll does not end in .vcf.gz in cstable\n";exit 1; }
	if ($line_data[5] eq "none") { #no filtering if "none" in annotations table column
		`cp ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename_oneclassfilter_filtered.vcf`
	} else {
    	print "perl VcfOneClassFiltering_v3.3.pl ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf.gz_indivintersect2platforms.vcf ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf ${FILEPATH}$line_data[5] ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename_oneclassfilter $line_data[2]\n";
    	print `perl VcfOneClassFiltering_v3.3.pl ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf.gz_indivintersect2platforms.vcf ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename.vcf ${FILEPATH}$line_data[5] ${FILEPATH}VCF${vcfno}_${vcfstart}_nohomref_samplename_oneclassfilter $line_data[2]`;
	}
	$m=$m+1;

}
close CALLSETTABLE;
 
