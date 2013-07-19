#!/usr/bin/perl -w
#
# VcfCombineUGHaplo.pl - Combine genotype information from GATK UnifiedGenotyper and HaplotypeCaller
# Prior to running this:
#   Find union of indels and MNPs using GATK UnifiedGenotyper and HaplotypeCaller for each dataset
#   Trigger assembly at union of indel and MNP sites with HaplotypeCaller and force UnifiedGenotyper to call at union of SNP, indel, and MNP sites
# If HaplotypeCaller and UnifiedGenotyper agree at a site for a dataset, then use UnifiedGenotyper info for integration. If HaplotypeCaller and UnifiedGenotyper disagree at a site for a dataset, then: (1) for SNPs, if no nearby (<=20bp) indels/SNPs differ between HaplotypeCaller and UnifiedGenotyper, then use the UnifiedGenotyper call, or (2) for SNPs, if nearby (<=20bp) indels/SNPs differ between HaplotypeCaller and UnifiedGenotyper, then use the HaplotypeCaller call, or (3) for indels and MNPs, use the HaplotypeCaller call (or filter UnifiedGenotyper call)
# After running this script, run xxx.pl to:
#   (1) If datasets disagree and all datasets are using the UnifiedGenotyper call, then use normal arbitration process for SNPs or indels. (2) If datasets disagree and some datasets use the HaplotypeCaller call, then (a) if HaplotypeCaller calls all agree then use that call or (b) if HaplotypeCaller calls disagree then use arbitration process with HaplotypeCaller covariates, maybe emphasizing read length and homopolymer.
#
#
# Version 0.1 - Jan 10, 2013
# Version 0.3 - May 6, 2013 - changed so that positions that are reference in UG and don't have a nearby Haplotypecaller indel call have UGHap=UG instead of both, which caused correct UGdups calls to be ignored
#

use strict;
use warnings;


my $line;

# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != 2) {
 print "usage: VcfCombineUGHaplo.pl UnifiedGenotyper.vcf HaplotypeCaller.vcf output.vcf\n";
    print "example: perl /Applications/bioinfo/perl/VcfCombineUGHaplo.pl OMNI_CG_1KG_SOL_ILL_454_HS_ILLCLIA_XIll_IllPCRFree_IonEx_XPSolWGLS_HapCallAll_conf2_mbq10_raw_merged_allvar_121205_HSWGcallUG_1.vcf OMNI_CG_1KG_SOL_ILL_454_HS_ILLCLIA_XIll_IllPCRFree_IonEx_XPSolWGLS_HapCall_conf2_mbq10_raw_merged_IndelMNPvartrigger_121201_HSWGcall_1.vcf OMNI_CG_1KG_SOL_ILL_454_HS_ILLCLIA_XIll_IllPCRFree_IonEx_XPSolWGLS_HapCall_conf2_mbq10_raw_merged_IndelMNPvartrigger_130111_HSWGcall_UGHapMerge_1.vcf\n";
    exit;
}

# Assign the inputs
my $infileUG = "$ARGV[0]";
my $infileHaplo = "$ARGV[1]";
my $outfile = "$ARGV[2]";


# open files
unless ( open(INUG, "$infileUG")) {
    print "\nCannot open the file: $infileUG! \n\n";
    exit;
}
unless ( open(INHAP, "$infileHaplo")) {
    print "\nCannot open the file: $infileHaplo! \n\n";
    exit;
}
unless ( open(OUTPUT, ">$outfile")) {
    print "\nCannot open the file: $outfile to write to! \n\n";
    exit;
}

my $x=0;
while ($line = <INUG>) {
    if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, with last one being the line starting with #
    print OUTPUT $line;
    if ($x==1) {last;}
}


$x=0;
while ($line = <INHAP>) {
    if ($line =~ /^\#\#/) { } else {last;} #skip header lines, with last one being the line starting with #
}


my @fields;
my $lineUG;
my $chromposUG=0;
my $chromUG=1;
my $posUG=0;
my $idUG;
my $refUG;
my $altUG;
my $qualUG;
my $filterUG;
my $infoUG;
my $formatUG;
my $charUG;
my $gtUG;
my $varTypeUG;

my $lineHap;
my $chromposHap=0;
my $chromHap=1;
my $posHap=0;
my $idHap;
my $refHap;
my $altHap;
my $qualHap;
my $filterHap;
my $infoHap;
my $formatHap;
my $charHap;
my $gtHap;
my $varTypeHap="SNP";
my $varTypeHapprev="SNP";



my $endlineUG=1; #0 when reach end of file
my $endlineHap=1; #0 when reach end of file
my $UGread=1;  #1 when ready to read next VCF line
my $Hapread=1; #1 when ready to read next VCF line
my $chromprev=1; #previous chromosome name
my $posHapprev=-1; #previous Hap position
my $posUGprev=-1; #previous UG position
my $qualHapprev=0;

while($endlineUG==1 || $endlineHap==1) {


    if ($UGread==1) { readUG(); }
    if ($Hapread==1) { readHap(); }

    if ($posHapprev==$posHap || $posUGprev==$posUG) {last};

    if ($chromposUG==$chromposHap && $altUG eq $altHap && $gtUG eq $gtHap) { #same variant and genotype in both vcfs, so print UG
		print OUTPUT "$chromUG\t$posUG\t$idUG\t$refUG\t$altUG\t$qualUG\t$filterUG\t$infoUG;UGHap=both\t$formatUG\t$charUG";
        $UGread=1;
        $Hapread=1;
    } elsif ($chromposUG==$chromposHap && $altUG eq $altHap && !($gtUG eq $gtHap)) { #same variant but different genotype in vcfs, so print Hap
		print OUTPUT "$chromHap\t$posHap\t$idHap\t$refHap\t$altHap\t$qualHap\t$filterHap\t$infoHap;UGHap=HapVarUGdiff\t$formatHap\t$charHap";
        $UGread=1;
        $Hapread=1;
    } elsif ($chromposUG==$chromposHap && !($altUG eq $altHap)) { #different variants at same location in vcfs, so print Hap
		print OUTPUT "$chromHap\t$posHap\t$idHap\t$refHap\t$altHap\t$qualHap\t$filterHap\t$infoHap;UGHap=HapVarUGdiff\t$formatHap\t$charHap";
        $UGread=1;
        $Hapread=1;
    } elsif ($endlineUG==0 || (($chromposUG>$chromposHap))) { #no variant in UG at this position in Hap, so print Hap and increment Hap
		print OUTPUT "$chromHap\t$posHap\t$idHap\t$refHap\t$altHap\t$qualHap\t$filterHap\t$infoHap;UGHap=HapnoUG\t$formatHap\t$charHap";
        $UGread=0;
        $Hapread=1;
    } elsif (($endlineHap==0 || (($chromposUG<$chromposHap))) && ($gtUG eq "0/0" || $gtUG eq "./.\n") && ($posHap-$posUG>20 || $qualHap<=40 || $varTypeHap eq "SNP") && ($posUG-$posHapprev>20 || $qualHapprev<=40 || $varTypeHapprev eq "SNP")) { #no variant in Hap at this position in UG and UG genotype is HomRef or unknown and the closest variant called by the HaplotypeCaller >20 (or the variant quality <=40 or it's a SNP), so print UG and increment UG
		print OUTPUT "$chromUG\t$posUG\t$idUG\t$refUG\t$altUG\t$qualUG\t$filterUG\t$infoUG;UGHap=UG\t$formatUG\t$charUG";
        $UGread=1;
        $Hapread=0;
    } elsif (($endlineHap==0 || (($chromposUG<$chromposHap))) && ($gtUG eq "0/0" || $gtUG eq "./.\n") ) { #no variant in Hap at this position in UG and UG genotype is HomRef or unknown and Haplotypecaller has a confident indel call nearby, so print UG with UGHap=both and increment UG
		print OUTPUT "$chromUG\t$posUG\t$idUG\t$refUG\t$altUG\t$qualUG\t$filterUG\t$infoUG;UGHap=both\t$formatUG\t$charUG";
        $UGread=1;
        $Hapread=0;
#    } elsif (($endlineUG==0 || (($chromposUG<$chromposHap))) && ($posHap-$posUG>20 || $qualHap<=40) && ($posUG-$posHapprev>20 || $qualHapprev<=40) && $varTypeUG eq "SNP" ) { #no variant in Hap at this position in UG, UG genotype is not HomRef and the distance between this location and the closest variant called by the HaplotypeCaller >20 (or the variant quality <=40) and variant is a SNP, so print UG and increment UG
        #    } elsif (($endlineUG==0 || (($chromposUG<$chromposHap))) && ($posHap-$posUG>20 || $qualHap<=40) && ($posUG-$posHapprev>20 || $qualHapprev<=40) ) { #no variant in Hap at this position in UG, UG genotype is not HomRef and the distance between this location and the closest variant called by the HaplotypeCaller >20 (or the variant quality <=40), so print UG and increment UG
    } elsif (($endlineHap==0 || (($chromposUG<$chromposHap))) && ($posHap-$posUG>20 || $qualHap<=40 || $varTypeHap eq "SNP") && ($posUG-$posHapprev>20 || $qualHapprev<=40 || $varTypeHapprev eq "SNP") ) { #no variant in Hap at this position in UG, UG genotype is not HomRef and the distance between this location and the closest variant called by the HaplotypeCaller >20 (or the variant quality <=40 or it's a SNP), so print UG and increment UG
		print OUTPUT "$chromUG\t$posUG\t$idUG\t$refUG\t$altUG\t$qualUG\t$filterUG\t$infoUG;UGHap=UG\t$formatUG\t$charUG";
        $UGread=1;
        $Hapread=0;
    } else { ##no variant in Hap at this position in UG, UG genotype is not HomRef and the distance between this location and the closest variant called by the HaplotypeCaller >20 or the Haplotype variant is a SNP, so print filtered UG and increment UG
		print OUTPUT "$chromUG\t$posUG\t$idUG\t$refUG\t$altUG\t$qualUG\tnotInHaplotype\t$infoUG;UGHap=HapRefUGVar\t$formatUG\t$charUG";
        $UGread=1;
        $Hapread=0;
    }

}


close OUTPUT;
close INUG;
close INHAP;

exit;

sub readUG {
    $chromprev=$chromUG;
    $posUGprev=$posUG;
    $lineUG = <INUG>;
    #print $line; next;
    if (defined($lineUG)) {
    } else { #quit loop if reach end of file
        $endlineUG = 0;
        return;
    }

    # Split up the line into an array
    @fields = split( "\t", $lineUG);

    my $chrno=$fields[0];
    if ($chrno eq "M" || $chrno eq "MT") {
        $chrno=23;
    } elsif ($chrno eq "X") {$chrno=24;
    } elsif ($chrno eq "Y") {$chrno=25;
    } elsif ($chrno =~ /^\d*$/) {
    } else { $endlineUG = 0; return;}
    $chromposUG = $chrno*1000000000+$fields[1];

    $chromUG = $fields[0];
     $posUG = $fields[1];
     $idUG = $fields[2];
     $refUG = $fields[3];
     $altUG = $fields[4];
     $qualUG = $fields[5];
     $filterUG = $fields[6];
     $infoUG = $fields[7];
     $formatUG = $fields[8];
     $charUG = $fields[9];
    
    if (length($refUG)==1 && length($altUG)==1) { $varTypeUG = "SNP"; } else { $varTypeUG = "notSNP"; }
    
    my @formats = split( ":", $formatUG);
    my @chars = split( ":", $charUG);

    for (my $i=0; $i<=$#formats; $i++) {
        if ($formats[$i] eq "GT") { $gtUG=$chars[$i]; }
    }
    
    return;

}

sub readHap {
    $chromprev=$chromHap;
    $posHapprev=$posHap;
    $qualHapprev=$qualHap;
    $varTypeHapprev=$varTypeHap;
    $lineHap = <INHAP>;
    #print $line; next;
    if (defined($lineHap)) {
    } else { #quit loop if reach end of file
        $endlineHap = 0;
        return;
    }
    
    # Split up the line into an array
    @fields = split( "\t", $lineHap);
    
    my $chrno=$fields[0];
    if ($chrno eq "M" || $chrno eq "MT") {
        $chrno=23;
    } elsif ($chrno eq "X") {$chrno=24;
    } elsif ($chrno eq "Y") {$chrno=25;
    } elsif ($chrno =~ /^\d*$/) {
    } else { $endlineHap = 0; return;}
    $chromposHap = $chrno*1000000000+$fields[1];

    $chromHap = $fields[0];
    $posHap = $fields[1];
    $idHap = $fields[2];
    $refHap = $fields[3];
    $altHap = $fields[4];
    $qualHap = $fields[5];
    $filterHap = $fields[6];
    $infoHap = $fields[7];
    $formatHap = $fields[8];
    $charHap = $fields[9];

    my @formats = split( ":", $formatHap);
    my @chars = split( ":", $charHap);
    
    if (length($refHap)==1 && length($altHap)==1) { $varTypeHap = "SNP"; } else { $varTypeHap = "notSNP"; }

    for (my $i=0; $i<=$#formats; $i++) {
        if ($formats[$i] eq "GT") { $gtHap=$chars[$i]; }
    }

    return;
    
}

