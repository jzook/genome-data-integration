#To load following perl modules
use strict;
use warnings;
use Getopt::Long;

#Input parameters
my (@INPUT_VCF, $OUTPUT_FILE);
my ($verbose);

#Declare variables
GetOptions(
	'vcf=s{1,}' => \@INPUT_VCF,
	'v' => \$verbose,
	'o|output:s' => \$OUTPUT_FILE,
	"h|help|?" => \&usage);

#To check input files
my $INPUT_VCF = $INPUT_VCF[0];
if(defined($OUTPUT_FILE)){$OUTPUT_FILE = $OUTPUT_FILE} else {print usage();die "\n\nPlease specify the name of output file.\n\n"}

#To define usage of the cg_vcf
print "Usage = perl cg_vcf.pl -vcf $INPUT_VCF ";
print "-o $OUTPUT_FILE \n\n";

sub usage {
	print "\nusage: perl cg_vcf.pl -vcf <vcf> -o <TAB> \n";
	print "\t-v \t\t\tThe input VCF file.\n";
	print "\t-o, --output\t\t\t\tOutput file name. Output file format is a TAB format\n";
	exit 1;
}

open (VCF_FILE, "<", $INPUT_VCF) or die "Can't open the VCF file. Please check if the file exists\n\n";
open (OUT_INFO, ">", $OUTPUT_FILE) or die $!;

my $n = 0;

while (my $line = <VCF_FILE>) {
	chomp($line);
	my @line_data = ();
	@line_data = split("\t", $line);
    
    if($line =~ m/^\#/) {
        print OUT_INFO $line . "\n";
    }
    
    else {

		if((0+@line_data)==10) {
			my @line_geno = ();
			@line_geno = split(":", $line_data[9]);
			#exclude no calls and break-end SVs
			if(($line_data[4] =~ m/NOCALL/) || ($line_geno[0] =~ m/(\.)/) || ($line_data[5] =~ m/\[/) || ($line_data[5] =~ m/\]/)) {
				$n = $n+1;
			}
        
			else {
				print OUT_INFO $line . "\n";
			}			
		}
	}
}

