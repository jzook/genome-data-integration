#Annotate SNPs specific to one dataset after merging VCF
#v2 - run on datasets in csv automatically

#read in user parameters
my $csv = shift @ARGV;  #csv file
my $bammin = shift @ARGV;  #min bam file number to run
my $bammax = shift @ARGV;  #max bam file number to run

my $res=1;


$res=0;
# Open the csv file
unless ( open(INFILE, "${csv}") ) {
    print "\nCould not open file: $csv!\n\n";
    exit;
}

my $i=1;
while(defined($line = <INFILE>)) {

    my @fields = split( ",", $line);
my $parallel = "&"; 
#if ($i % 3) { $parallel = ""; } #run 3 in parallel at a time
if ($fields[3] <= $bammax && $fields[3] >= $bammin && $res == 0) {
    my $samCmd = "perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_conf2_mbq10_raw_merged_SNPs_120427_new_${fields[1]}annotate ${fields[2]} ${fields[3]} ${fields[4]}";
    print "$samCmd\n";
    $res = system($samCmd);
}   

}

exit;
