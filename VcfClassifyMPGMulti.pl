#!/usr/bin/perl -w
#
# VcfClassifyMPGMulti.pl - Classify highly confident genotypes using MostProbableGenotype annotation from GATK
# Highly confident Homozygous requires net MPGHom>12, net MPGHom/DP>0.16, and no individual MPGHom<-2
#	Note that MPGHom/DP>0.16 allows up to 5% mismatching bases with BQ=25
# Highly confident Heterozygous requires net MPGHet>20, net MPGHet/DP>0.68, and no individual MPGHet<-2
#	Note that MPGHet/DP>0.68 allows allele balance between 40 and 60% with all BQ=25
#
#
# Version 1.0.0 - Mar 1, 2012
#

use strict;
use warnings;


my $line;

# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV < 2) {
 print "usage: VcfClassifyMPGMulti.pl inputvcffilenamestart datasetNames(>=2 separated by spaces without call.vcf)\n";
 print "Example: perl /Volumes/SSD960/workspace/data/vcfannot/VcfClassifyMPGMulti.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ HSWG ILLWG XIll SOLWG 454WG CG ILLCLIA IllPCRFreeNoBQSR XPSolWGLS\n";
 exit;
}

# Assign the inputs
my $vcfcnt = 1;

	my $infilestart = "$ARGV[0]";
my @infiles;
my @infilehandlesDist;
my @outfilehandlesHomRefUncert;
my @outfilehandlesHetRefUncert;
my @outfilehandlesHomVarUncert;
while ($vcfcnt <= $#ARGV) {
	$infiles[$vcfcnt-1] = $ARGV[$vcfcnt];
	local *FILEIN;
	open (FILEIN, "${infilestart}${infiles[$vcfcnt-1]}call.vcf") || die;
	push(@infilehandlesDist, *FILEIN);
	local *FILE;
	open (FILE, ">${infilestart}${infiles[$vcfcnt-1]}call_HomRefUncert.vcf") || die;
	push(@outfilehandlesHomRefUncert, *FILE);
	local *FILE2;
	open (FILE2, ">${infilestart}${infiles[$vcfcnt-1]}call_HetRefUncert.vcf") || die;
	push(@outfilehandlesHetRefUncert, *FILE2);
	local *FILE3;
	open (FILE3, ">${infilestart}${infiles[$vcfcnt-1]}call_HomVarUncert.vcf") || die;
	push(@outfilehandlesHomVarUncert, *FILE3);
	$vcfcnt++;
}
$vcfcnt -= 1;


# Output file is opened
unless ( open(OUTPUT, ">${infilestart}allcall.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTHOMREF, ">${infilestart}allcall_HomRef.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_HomRef.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTHET, ">${infilestart}allcall_HetRef.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_HetRef.vcf to write to! \n\n";
exit;
}
unless ( open(OUTHOMVAR, ">${infilestart}allcall_HomVar.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_HomVar.vcf to write to! \n\n";
    exit;
}
    #skip header lines in individual files and write first one to output
    my $fhn = 0;
    foreach my $fh (@infilehandlesDist) {
	    my $x=0;
        while ($line = <$fh>) {
            if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
		if ($fhn==0) {print OUTPUT $line;}
		if ($fhn==0) {print OUTHOMREF $line;}
		if ($fhn==0) {print OUTHET $line;}
		if ($fhn==0) {print OUTHOMVAR $line;}
			my $fho = $outfilehandlesHomRefUncert[$fhn];
            	print $fho $line;
                        $fho = $outfilehandlesHetRefUncert[$fhn];
            	print $fho $line;
                        $fho = $outfilehandlesHomVarUncert[$fhn];
            	print $fho $line;
        if ($x==1) {last;}
        }
		$fhn++;
		
    }

my @fields;
my @lines;
my $info;



my $endline=1;
while($endline==1) {


    
            my $MPGHomRefTot=0; my $MPGHetTot=0; my $MPGHomVarTot=0;
            my $NoHomRef=0; my $NoHet=0; my $NoHomVar=0;
            my $MQ0=0; my $DPTot=0;

    my $vcfno = 0; #number of vcf file
    foreach my $fh (@infilehandlesDist) {
	    $line = <$fh>;
	    #print $line; next;
		if (defined($line)) {
		} else { #quit loop if reach end of file
			$endline = 0;
			last;
		}
		$lines[$vcfno] = $line;
		$vcfno++;
	
		# Split up the line into an array 
		@fields = split( "\t", $line);
	
		$info = $fields[7]; #extract INFO field

	
	#print "$info\n";
			my $tmp;
		if ($info =~ /MPGLikHomRef=(.*?);/) { 
			$tmp=$1;
			if ($1 =~ /Inf/) {next;}
			$MPGHomRefTot += $tmp;
			#print OUTPUT "MHR=$tmp,"; #test
		} 
		if ($info =~ /MPGLikHomVar=(.*?);/) { 
			my $tmp2=$1;
			if ($1 =~ /Inf/) {$MPGHomRefTot -= $tmp; next;}
			$MPGHomVarTot += $tmp2;
		} 
		if ($info =~ /MPGLikHetRef=(.*?);/) { 
			$MPGHetTot += $1;
		} 
		#if ($info =~ /MPGHomRef=((-|{0-9}|\.)*?);/) { 
		if ($info =~ /MPGHomRef=(.*?);/) { 
		    if  ($1 < -2) {$NoHomRef+=1;}
		} 
		if ($info =~ /MPGHetRef=(.*?);/) { 
		    if  ($1 < -2) {$NoHet+=1;}
		} 
		if ($info =~ /MPGHomVar=(.*?);/) { 
		    if  ($1 < -2) {$NoHomVar+=1;}
		} 
		if ($info =~ /DP=(.*?);/) { 
			$DPTot += $1;
		} 
        if ($info =~ /MQ0Fraction=(.*?);/) {
		    if ($1>0.2) {$MQ0 += 1;}
        }

    }
    #last;
    my $MPGHetRat;
    if ($MPGHomRefTot>$MPGHomVarTot) {
    	$MPGHetRat=$MPGHomRefTot-$MPGHetTot;
    } else {
    	$MPGHetRat=$MPGHomVarTot-$MPGHetTot;
    }
    my $geno=0;
    if ($MQ0/$vcfno<0.5 && $MPGHetTot-$MPGHomRefTot > 12 && ($MPGHetTot-$MPGHomRefTot)/$DPTot > 0.16 && $NoHomRef/$vcfno<0.25) {
    	$geno = 1;
	    my $vcfno = 0; #number of vcf file
		foreach my $fh (@outfilehandlesHomRefUncert) {
		
			# Split up the line into an array 
			@fields = split( "\t", $lines[$vcfno]);
			print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7];geno=$geno\t$fields[8]\t$fields[9]";
	
			$vcfno++;
		}
		print OUTHOMREF "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\tgeno=$geno\t$fields[8]\t$fields[9]";

    } elsif ($MQ0/$vcfno<0.5 && $MPGHetTot-$MPGHomVarTot > 12 && ($MPGHetTot-$MPGHomVarTot)/$DPTot > 0.16 && $NoHomVar/$vcfno<0.25) {
    	$geno = 3;
	    my $vcfno = 0; #number of vcf file
		foreach my $fh (@outfilehandlesHomVarUncert) {
		
			# Split up the line into an array 
			@fields = split( "\t", $lines[$vcfno]);
			print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7];geno=$geno\t$fields[8]\t$fields[9]";
	
			$vcfno++;
		}
		print OUTHOMVAR "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\tgeno=$geno\t$fields[8]\t$fields[9]";

    } elsif ($MQ0/$vcfno<0.3 && $MPGHetRat > 20 && $MPGHetRat/$DPTot > 0.68 && $NoHet/$vcfno<0.15) {
    	$geno = 2;
	    my $vcfno = 0; #number of vcf file
		foreach my $fh (@outfilehandlesHetRefUncert) {
		
			# Split up the line into an array 
			@fields = split( "\t", $lines[$vcfno]);
			print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7];geno=$geno\t$fields[8]\t$fields[9]";
	
			$vcfno++;
		}
		print OUTHET "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\tgeno=$geno\t$fields[8]\t$fields[9]";

    } else {
		$geno = 0;
		my $vcfno = 0; #number of vcf file
		foreach my $fh (@outfilehandlesHomRefUncert) {

                        # Split up the line into an array
	    @fields = split( "\t", $lines[$vcfno]);
	    print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7];geno=$geno\t$fields[8]\t$fields[9]";

	    $vcfno++;
	}
	
         $vcfno = 0; #number of vcf file
        foreach my $fh (@outfilehandlesHomVarUncert) {

                        # Split up the line into an array
            @fields = split( "\t", $lines[$vcfno]);
            print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7];geno=$geno\t$fields[8]\t$fields[9]";

            $vcfno++;
        }

         $vcfno = 0; #number of vcf file
        foreach my $fh (@outfilehandlesHetRefUncert) {

                        # Split up the line into an array
            @fields = split( "\t", $lines[$vcfno]);
            print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7];geno=$geno\t$fields[8]\t$fields[9]";

            $vcfno++;
        }
    }


		print OUTPUT "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\tgeno=$geno;MPGLikHomRefTot=$MPGHomRefTot;MPGLikHetRefTot=$MPGHetTot;MPGLikHomVarTot=$MPGHomVarTot;NoHomRef=$NoHomRef;NoHetRef=$NoHet;NoHomVar=$NoHomRef;DPTot=$DPTot;LowMQ0Fraction=$MQ0\t$fields[8]\t$fields[9]";

#last;
}
 

close OUTPUT;
close OUTHOMREF;
close OUTHET;
close OUTHOMVAR;
    foreach my $fh (@infilehandlesDist) {
        close $fh;
    }
    foreach my $fh (@outfilehandlesHomRefUncert) {
        close $fh;
    }
    foreach my $fh (@outfilehandlesHetRefUncert) {
        close $fh;
    }
    foreach my $fh (@outfilehandlesHomVarUncert) {
        close $fh;
    }

exit;
