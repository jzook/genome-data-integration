#!/usr/bin/perl -w
#
# VcfHighConfUGHaploMulti.pl - find high confidence SNPs and indels
# Prior to running this:
#   Find union of indels and MNPs using GATK UnifiedGenotyper and HaplotypeCaller for each dataset
#   Trigger assembly at union of indel and MNP sites with HaplotypeCaller and force UnifiedGenotyper to call at union of SNP, indel, and MNP sites
#   Run VcfCombineUGHaplo.pl for each dataset - Combine genotype information from GATK UnifiedGenotyper and HaplotypeCaller
# Run this to find high confidence SNPs and indels:
#   (1) If datasets disagree and all datasets are using the UnifiedGenotyper call, then use normal arbitration process for SNPs or indels. (2) If datasets disagree and some datasets use the HaplotypeCaller call, then (a) if HaplotypeCaller calls all agree then use that call or (b) if HaplotypeCaller calls disagree then use arbitration process with HaplotypeCaller covariates, maybe emphasizing read length and homopolymer.
#
#
# Version 1.0.0 - Jan 10, 2013
# Version 1.1 - Feb 17, 2013 - if there are multiple lines in the vcf at one position, then find the most probable and discard the others
# Ver 1.2 - Apr 8, 2013 - Don't use dups call if haplotypecaller makes a call since dup call would override HomRef haplotypecaller calls.
#			8/26/13 - require that a haplotypecaller reference base be used in >20% of the datasets since otherwise we sometimes have a 454 or Ion homopolymer indel override a real snp in other platforms
#			8/28/13 - add export of single alt variants to csv for machine learning
#			11/8/13 - add one to $HapRefCall if $info[$kk] =~ /UGHap=HapRefUGVar/ to avoid calling incorrect variants near complex variants 
#			2/9/14	- look for HapCall==both or ==Both
use strict;
use warnings;
use List::Util qw[min max];

my $xx=0; #debug counter
my $line;
my @ALTPL=(3,6,10,15); #number of PLs for each number of ALT alleles
    #HSWG ILLWG XIll 454WG ILLCLIA IllPCRFree XPSolWGLS IonEx HSWEx ILLWEx CG PlatGen ILL250
    my @exp = (3,4,6,10,11,12); #datasets to export as csv

#F(j/k) = (k*(k+1)/2)+j -> formula for determining PL order from genotype
my @GTs=("0/0","0/1","1/1","0/2","1/2","2/2","0/3","1/3","2/3","3/3","0/4","1/4","2/4","3/4","4/4");

# Check if the right number of paramters are given at start up
#
# The following form is needed at start up:

if ($#ARGV < 3) {
 print "usage: VcfHighConfUGHaploMulti.pl inputvcfUGHapMergefilenamestart chromosomeName datasetNames(>=2 separated by spaces without call.vcf)\n";
    print "example: perl /Applications/bioinfo/perl/VcfHighConfUGHaploMulti_1.1.pl AllFDAdatasets_130204_ 5 HSWG ILLWG XIll 454WG ILLCLIA IllPCRFree XPSolWGLS IonEx HSWEx ILLWEx CG\n";
    exit;
}

# Assign the inputs
my $vcfcnt = 2;

my $infilestart = "$ARGV[0]";
my $chrom = "$ARGV[1]";
my @infiles;
my @infilehandlesDist;
my @infilehandlesDups;
my @outfilehandlesHomRefUncert;
my @outfilehandlesHetRefUncert;
my @outfilehandlesHomVarUncert;
my @outfilehandlesHomUncert;
while ($vcfcnt <= $#ARGV) {
	$infiles[$vcfcnt-2] = $ARGV[$vcfcnt];
	local *FILEIN;
	open (FILEIN, "${infilestart}${infiles[$vcfcnt-2]}call_UGHapMerge_${chrom}.vcf") || die;
	push(@infilehandlesDist, *FILEIN);
	local *FILEIN2;
	open (FILEIN2, "${infilestart}${infiles[$vcfcnt-2]}call_UGdups_${chrom}.vcf") || die;
	push(@infilehandlesDups, *FILEIN2);
#	local *FILE;
#	open (FILE, ">${infilestart}${infiles[$vcfcnt-2]}call_UGHapMerge_HomRefUncert_${chrom}.vcf") || die;
#	push(@outfilehandlesHomRefUncert, *FILE);
	local *FILE2;
	open (FILE2, ">${infilestart}${infiles[$vcfcnt-2]}call_UGHapMerge_HetRefUncert_${chrom}.vcf") || die;
	push(@outfilehandlesHetRefUncert, *FILE2);
#	local *FILE3;
#	open (FILE3, ">${infilestart}${infiles[$vcfcnt-2]}call_UGHapMerge_HomVarUncert_${chrom}.vcf") || die;
#	push(@outfilehandlesHomVarUncert, *FILE3);
	local *FILE3;
	open (FILE3, ">${infilestart}${infiles[$vcfcnt-2]}call_UGHapMerge_HomUncert_${chrom}.vcf") || die;
	push(@outfilehandlesHomUncert, *FILE3);
	$vcfcnt++;
}
$vcfcnt -= 2;


# Output file is opened
unless ( open(OUTPUT, ">${infilestart}allcall_UGHapMerge_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTHOMREF, ">${infilestart}allcall_UGHapMerge_HomRef_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_HomRef_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTHET, ">${infilestart}allcall_UGHapMerge_HetRef_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_HetRef_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(OUTHOMVAR, ">${infilestart}allcall_UGHapMerge_HomVar_${chrom}.vcf")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_HomVar_${chrom}.vcf to write to! \n\n";
    exit;
}
unless ( open(EXPCSV, ">${infilestart}allcall_UGHapMerge_${chrom}.csv")) {
    print "\nCannot open the file: ${infilestart}allcall_UGHapMerge_${chrom}.csv to write to! \n\n";
    exit;
}

#print header for CSV
if ($chrom eq "1a") {
print EXPCSV "chrompos,ref,alt,genotypeClass";
for (my $i=0; $i<=$#exp; $i++) {
	print EXPCSV ",${infiles[$exp[$i]]}_PL0,${infiles[$exp[$i]]}_PL1,${infiles[$exp[$i]]}_PL2,${infiles[$exp[$i]]}_PLbyDP,${infiles[$exp[$i]]}_DP,${infiles[$exp[$i]]}_FS,${infiles[$exp[$i]]}_MQ,${infiles[$exp[$i]]}_NBQ,${infiles[$exp[$i]]}_OND,${infiles[$exp[$i]]}_ReadPosEndDist,${infiles[$exp[$i]]}_HaplotypeScore,${infiles[$exp[$i]]}_BaseQRankSum";    		
}
print EXPCSV "\n";
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
        #print $fho $line;
        $fho = $outfilehandlesHetRefUncert[$fhn];
        print $fho $line;
        #$fho = $outfilehandlesHomVarUncert[$fhn];
        #print $fho $line;
        $fho = $outfilehandlesHomUncert[$fhn];
        print $fho $line;
        if ($x==1) {last;}
    }
    $fhn++;
    
}
foreach my $fh (@infilehandlesDups) {
    my $x=0;
    while ($line = <$fh>) {
        if ($line =~ /^\#\#/) { } else {$x=1;} #skip header lines, but including line starting with #
        if ($x==1) {last;}
    }
    
}

my @fields;

my $info;

my @line=((0) x ($#infilehandlesDist+1));
my @chrompos=((0) x ($#infilehandlesDist+1));
my @chromposprev=((0) x ($#infilehandlesDist+1));
my @chrom=((0) x ($#infilehandlesDist+1));
my @pos=((0) x ($#infilehandlesDist+1));
my @id=((0) x ($#infilehandlesDist+1));
my @ref=((0) x ($#infilehandlesDist+1));
my @alt=((0) x ($#infilehandlesDist+1));
my @qual=((0) x ($#infilehandlesDist+1));
my @filter=((0) x ($#infilehandlesDist+1));
my @info=((0) x ($#infilehandlesDist+1));
my @format=((0) x ($#infilehandlesDist+1));
my @char=((0) x ($#infilehandlesDist+1));
my @gt=((0) x ($#infilehandlesDist+1));
my @PLtext=((0) x ($#infilehandlesDist+1));
my @varType=((0) x ($#infilehandlesDist+1));
    my @hapCall=((0) x ($#infilehandlesDist+1));

my @linedups=((0) x ($#infilehandlesDist+1));
my @chromposdups=((0) x ($#infilehandlesDist+1));
my @chromposdupsprev=((0) x ($#infilehandlesDist+1));
my @chromdups=((0) x ($#infilehandlesDist+1));
my @posdups=((0) x ($#infilehandlesDist+1));
my @iddups=((0) x ($#infilehandlesDist+1));
my @refdups=((0) x ($#infilehandlesDist+1));
my @altdups=((0) x ($#infilehandlesDist+1));
my @qualdups=((0) x ($#infilehandlesDist+1));
my @filterdups=((0) x ($#infilehandlesDist+1));
my @infodups=((0) x ($#infilehandlesDist+1));
my @formatdups=((0) x ($#infilehandlesDist+1));
my @chardups=((0) x ($#infilehandlesDist+1));
my @gtdups=((0) x ($#infilehandlesDist+1));
my @PLtextdups=((0) x ($#infilehandlesDist+1));
my @varTypedups=((0) x ($#infilehandlesDist+1));



my $minchrompos=0; #previous minimum chrompos
my $newline=1; #printed out a new line, so ready for next one

my @endline = ((1) x ($#infilehandlesDist+1));
my @endlinedups = ((1) x ($#infilehandlesDist+1));
my $endlineall=1;
while($endlineall==1) {
    
    
    
    my @PLTot=((0) x 15);    #allows up to 4 ALT alleles
    my @NoPLTot=((0) x 15);
    my @YesPLTot=((0) x 15);
    my $MQ0=0; my $DPTot=0; my $PLno=0;
    
    my $dataset=-1;
    foreach my $fh (@infilehandlesDist) {
        $dataset++;
		if ($endline[$dataset]==1 && $chrompos[$dataset]==$minchrompos) {
            $chromposprev[$dataset]=$chrompos[$dataset];
            $line[$dataset] = <$fh>;
            #print $line[$dataset];
            if (defined($line[$dataset])) {
            } else { #quit loop if reach end of file
                $endline[$dataset] = 0;
                $chrompos[$dataset] = 1000000000000;
                next;
            }
            
            # Split up the line into an array
            @fields = split( "\t", $line[$dataset]);
            my $chrno=$fields[0];
            if ($chrno eq "M" || $chrno eq "MT") {
                $chrno=23;
            } elsif ($chrno eq "X") {$chrno=24;
            } elsif ($chrno eq "Y") {$chrno=25;
            } elsif ($chrno =~ /^\d*/) {
            } else {
                $endline[$dataset] = 0; $chrompos[$dataset] = 1000000000000;next;
            }
            $chrompos[$dataset] = $chrno*1000000000+$fields[1];
            #print "chrpos=$chrompos[$dataset]\n";
            $chrom[$dataset] = $fields[0];
            $pos[$dataset] = $fields[1];
            $id[$dataset] = $fields[2];
            $ref[$dataset] = $fields[3];
            $alt[$dataset] = $fields[4];
            $qual[$dataset] = $fields[5];
            if ($qual[$dataset] eq ".") { $qual[$dataset]=0;}
            $filter[$dataset] = $fields[6];
            $info[$dataset] = $fields[7];
            $format[$dataset] = $fields[8];
            $char[$dataset] = $fields[9];

            my @formats = split( ":", $format[$dataset]);
            my @chars = split( ":", $char[$dataset]);
            $PLtext[$dataset]="";
            for (my $i=0; $i<=$#formats; $i++) {
                if ($formats[$i] eq "GT") { $gt[$dataset]=$chars[$i]; }
                if ($formats[$i] eq "PL") {
                    if ($chars[$i] =~ /(.*)\n/) {$PLtext[$dataset]=$1;}
                    else { $PLtext[$dataset]=$chars[$i];}
                }
            }
            if (($qual[$dataset]==0 || $qual[$dataset]>40) && ($info[$dataset] =~ /UGHap=HapnoUG/ ||  $info[$dataset] =~ /UGHap=HapVarUGdiff/ || $info[$dataset] =~ /UGHap=both/ || $info[$i] =~ /UGHap=Both/)) { $hapCall[$dataset]=1;
            } else { $hapCall[$dataset]=0;}
            
            #if ($chrompos[$dataset]==1020179655) {
            #	print "$line[$dataset]";
            #	print "$qual[$dataset],$hapCall[$dataset],$info[$dataset]\n";
            #	if ($info[$dataset] =~ /UGHap=HapVarUGdiff/) { print "HapVarUGdiff\n\n"; }
            #}
            

		}
    }

    $dataset=-1;
    foreach my $fh (@infilehandlesDups) {
        $dataset++;
		if ($endlinedups[$dataset]==1 && $chromposdups[$dataset]==$minchrompos) {
            $chromposdupsprev[$dataset]=$chromposdups[$dataset];
            $linedups[$dataset] = <$fh>;
            #print $linedups[$dataset];
            if (defined($linedups[$dataset])) {
            } else { #quit loop if reach end of file
                $endlinedups[$dataset] = 0;
                $chromposdups[$dataset] = 1000000000000;
                next;
            }
            
            # Split up the line into an array
            @fields = split( "\t", $linedups[$dataset]);
            my $chrno=$fields[0];
            if ($chrno eq "M" || $chrno eq "MT") {
                $chrno=23;
            } elsif ($chrno eq "X") {$chrno=24;
            } elsif ($chrno eq "Y") {$chrno=25;
            } elsif ($chrno =~ /^\d*/) {
            } else { $endlinedups[$dataset] = 0; $chromposdups[$dataset] = 1000000000000; next;}
            $chromposdups[$dataset] = $chrno*1000000000+$fields[1];
            #print "chrpos=$chromposdups[$dataset]\n";
            $chromdups[$dataset] = $fields[0];
            $posdups[$dataset] = $fields[1];
            $iddups[$dataset] = $fields[2];
            $refdups[$dataset] = $fields[3];
            $altdups[$dataset] = $fields[4];
            $qualdups[$dataset] = $fields[5];
            if ($qualdups[$dataset] eq ".") { $qualdups[$dataset]=0;}
            $filterdups[$dataset] = $fields[6];
            $infodups[$dataset] = $fields[7];
            $formatdups[$dataset] = $fields[8];
            $chardups[$dataset] = $fields[9];
            
            my @formats = split( ":", $formatdups[$dataset]);
            my @chars = split( ":", $chardups[$dataset]);
            $PLtextdups[$dataset]="";
            for (my $i=0; $i<=$#formats; $i++) {
                if ($formats[$i] eq "GT") { $gtdups[$dataset]=$chars[$i]; }
                if ($formats[$i] eq "PL") {
                    if ($chars[$i] =~ /(.*)\n/) {$PLtextdups[$dataset]=$1;}
                    else { $PLtextdups[$dataset]=$chars[$i];}
                }
            }
		}
    }
#    if ($minchrompos>2121992900) {
#        print "$endline[0],$endline[1],$minchrompos,@chrompos\n";
#        print "$endlinedups[0],$endlinedups[1],$minchrompos,@chromposdups\n\n";
#    }
    
    if (max(@endline)+max(@endlinedups)==0) { last; }
    $minchrompos = min(min(@chrompos),min(@chromposdups));
    
    #if ($minchrompos==1003765095) { print "$minchrompos\n";}

    ###TODO: take into account multiple Ref/Alt combinations
    #Currently, take predominant Ref and integrate different Alts for this ref and ignore other refs
    my @altallVar=(("") x 6);
    my $varType="SNP";
    my $altno=0; my $altVarLong=""; my $refVarLong=""; my $multiplealts=0; my $refmax="";
    my $altno1=0; my $altVarLong1="";
    my $altno2=0; my $altVarLong2="";
    my $altno3=0; my $altVarLong3="";
    my $ref1=""; my $ref1cnt=0; my $ref2=""; my $ref2cnt=0; my $ref3=""; my $ref3cnt=0; my $ref1hapcnt=0; my $ref2hapcnt=0; my $ref3hapcnt=0;
    my $hapCallRef=0;
    for (my $kk=0; $kk<=$#infilehandlesDist; $kk++) {
        if ($chrompos[$kk]==$minchrompos && $ref1 eq "") {
            my @alts = split( ",", $alt[$kk]);
            $ref1=$ref[$kk]; $altno1=$#alts+1; $altVarLong1=$alt[$kk];
        }
        if ($chrompos[$kk]==$minchrompos && $hapCall[$kk]==1 && $gt[$kk] eq "./.\n") {
            my @alts = split( ",", $alt[$kk]);
            if (length($ref[$kk])>1 && length($alts[0])>1) { $hapCallRef++; }
        }
        if ($info[$kk] =~ /UGHap=HapRefUGVar/) {
            $hapCallRef++;
            next;
        }
        if ($chrompos[$kk]==$minchrompos && !($qual[$kk] eq ".") && $qual[$kk]>10) {
            my @alts = split( ",", $alt[$kk]);
            
            if ($ref[$kk] eq $ref1) {
                $ref1cnt++;
                if ($hapCall[$kk]==1) {$ref1hapcnt++;}
                if ($#alts+1>$altno1) {$altno1=$#alts+1; $altVarLong1=$alt[$kk]; }
            } elsif ($ref2 eq "" || $ref[$kk] eq $ref2) {
                $ref2=$ref[$kk];
                $ref2cnt++;
                if ($hapCall[$kk]==1) {$ref2hapcnt++;}
                if ($#alts+1>$altno2) {$altno2=$#alts+1; $altVarLong2=$alt[$kk]; }
            } elsif ($ref3 eq "" || $ref[$kk] eq $ref3) {
                $ref3=$ref[$kk];
                $ref3cnt++;
                if ($hapCall[$kk]==1) {$ref3hapcnt++;}
                if ($#alts+1>$altno3) {$altno3=$#alts+1; $altVarLong3=$alt[$kk]; }
            }
        }
        if (($hapCall[$kk]==0 || $chrompos[$kk]!=$minchrompos) && $chromposdups[$kk]==$minchrompos && !($qualdups[$kk] eq ".") && $qualdups[$kk]>20) {
            my @alts = split( ",", $altdups[$kk]);
            
            if ($refdups[$kk] eq $ref1) {
                $ref1cnt++;
                if ($#alts+1>$altno1) {$altno1=$#alts+1; $altVarLong1=$altdups[$kk]; }
            } elsif ($ref2 eq "" || $refdups[$kk] eq $ref2) {
                $ref2=$refdups[$kk];
                $ref2cnt++;
                if ($#alts+1>$altno2) {$altno2=$#alts+1; $altVarLong2=$altdups[$kk]; }
            } elsif ($ref3 eq "" || $refdups[$kk] eq $ref3) {
                $ref3=$refdups[$kk];
                $ref3cnt++;
                if ($#alts+1>$altno3) {$altno3=$#alts+1; $altVarLong3=$altdups[$kk]; }
            }
        }
    }
    
    if ($minchrompos==1020179655) {
        print "$ref1hapcnt,$ref2hapcnt,$ref3hapcnt,$hapCallRef\n";
        print "$ref1cnt,$ref2cnt,$ref3cnt,$ref1,$ref2,$ref2\n";
        print "$hapCall[0],$hapCall[1],$hapCall[2],$hapCall[3],$hapCall[4],$hapCall[5],$hapCall[6],$hapCall[7],$hapCall[8],$hapCall[9],$hapCall[10],$hapCall[11],$hapCall[12],$hapCall[13]\n";
        print "$gt[0],$gt[1],$gt[2],$gt[3],$gt[4],$gt[5],$gt[6],$gt[7],$gt[8],$gt[9],$gt[10],$gt[11],$gt[12],$gt[13]\n";
        print "$ref[0],$ref[1],$ref[2],$ref[3],$ref[4],$ref[5],$ref[6],$ref[7],$ref[8],$ref[9],$ref[10],$ref[11],$ref[12],$ref[13]\n";
        print "$alt[0];$alt[1];$alt[2];$alt[3];$alt[4];$alt[5];$alt[6];$alt[7];$alt[8];$alt[9];$alt[10];$alt[11];$alt[12];$alt[13]\n";
        print "$refdups[0],$refdups[1],$refdups[2],$refdups[3],$refdups[4],$refdups[5],$refdups[6],$refdups[7],$refdups[8],$refdups[9],$refdups[10],$refdups[11],$refdups[12],$refdups[13]\n";
        print "Chrompos:@chrompos\n";
    }
    my $HapCallVar=0;
    if (($ref3hapcnt+$ref2hapcnt+$ref1hapcnt<$hapCallRef) && ($hapCallRef>2 || $ref1cnt+$ref2cnt+$ref3cnt<4)) {
#        if ($minchrompos==1226460240) {
#            print "1next\n";
#        }
        #print "$minchrompos,";
        next; #skip this record if the HaplotypeCaller is HomRef and most records have no genotype call (except dups)
    } elsif ($ref3hapcnt>$ref2hapcnt && $ref3hapcnt>$ref1hapcnt && $ref3hapcnt>0.5*($ref1cnt+$ref2cnt+$ref3cnt)) {
        $refmax=$ref3; $altVarLong=$altVarLong3; $altno=$altno3; $HapCallVar=$ref3hapcnt;
    } elsif ($ref2hapcnt>$ref1hapcnt && $ref2hapcnt>0.5*($ref1cnt+$ref2cnt+$ref3cnt)) {
        $refmax=$ref2; $altVarLong=$altVarLong2; $altno=$altno2; $HapCallVar=$ref2hapcnt;
    } elsif ($ref1hapcnt>0.5*($ref1cnt+$ref2cnt+$ref3cnt)) {
        $refmax=$ref1; $altVarLong=$altVarLong1; $altno=$altno1; $HapCallVar=$ref1hapcnt;
    } elsif ($ref3cnt>$ref2cnt && $ref3cnt>$ref1cnt) {
        $refmax=$ref3; $altVarLong=$altVarLong3; $altno=$altno3;
    } elsif ($ref2cnt>$ref1cnt) {
        $refmax=$ref2; $altVarLong=$altVarLong2; $altno=$altno2;
    } else {#includes if no datasets have a qual>20, then use first datasets' ref/alt
        $refmax=$ref1; $altVarLong=$altVarLong1; $altno=$altno1;
    }
    
#    if ($minchrompos==1020179655) {
#        print "1\n";
#    }
    
    #F(j/k) = (k*(k+1)/2)+j -> formula for determining PL order from genotype
    #    my @GTs=("0/0","0/1","1/1","0/2","1/2","2/2","0/3","1/3","2/3","3/3","0/4","1/4","2/4","3/4","4/4");
    #my @ALTPL=(3,6,10,15); #number of PLs for each number of ALT alleles
    
    my @allalts = split( ",", $altVarLong);
    #if there are  different alt alleles for refmax, then combine them
    for (my $kk=0; $kk<=$#infilehandlesDist; $kk++) {
        if ($chrompos[$kk]==$minchrompos) {
            if ($ref[$kk] eq $refmax) {
                my @altsVar = split( ",", $alt[$kk]);
                for (my $ii=0; $ii<=$#altsVar; $ii++) {
                    my $match=0;
                    for (my $jj=0; $jj<=$#allalts; $jj++) {
                        if ($altsVar[$ii] eq $allalts[$jj]) { $match=1;}
                    }
                    if ($match==0 && $altno<4) {
                        $altVarLong="$altVarLong,$altsVar[$ii]";
                        @allalts = split( ",", $altVarLong);
                        $altno++;
                    }
                }
            }
        }
        if (($hapCall[$kk]==0 || $chrompos[$kk]!=$minchrompos) && $chromposdups[$kk]==$minchrompos) {
            if ($refdups[$kk] eq $refmax) {
                my @altsVar = split( ",", $altdups[$kk]);
                for (my $ii=0; $ii<=$#altsVar; $ii++) {
                    my $match=0;
                    for (my $jj=0; $jj<=$#allalts; $jj++) {
                        if ($altsVar[$ii] eq $allalts[$jj]) { $match=1;}
                    }
                    if ($match==0 && $altno<4) {
                        $altVarLong="$altVarLong,$altsVar[$ii]";
                        @allalts = split( ",", $altVarLong);
                        $altno++;
                    }
                }
            }
        }
    }
    if (length($refmax)>1) {$varType="INDEL";}
    my $j=0;
    @allalts = split( ",", $altVarLong);
    foreach my $alt (@allalts) {
        if (length($alt)>1) {$varType="INDEL";}
        $j++;
    }
    if ($altno>4) {$altno=4;} # only allow 4 ALT alleles
    $PLno=$ALTPL[$altno-1];
    
    if ($minchrompos==1020179655) {
        print "allalts:";
        foreach my $alt (@allalts) {
            print "$alt;";
        }
        print "\n";
    }

#    if ($minchrompos==1020179655) {
#        print "2\n";
#    }
    #store best call from dups and non-dups
    my @linebest=((0) x ($#infilehandlesDist+1));
    my @chromposbest=((0) x ($#infilehandlesDist+1));
    my @chrombest=((0) x ($#infilehandlesDist+1));
    my @posbest=((0) x ($#infilehandlesDist+1));
    my @idbest=((0) x ($#infilehandlesDist+1));
    my @refbest=((0) x ($#infilehandlesDist+1));
    my @altbest=((0) x ($#infilehandlesDist+1));
    my @qualbest=((0) x ($#infilehandlesDist+1));
    my @filterbest=((0) x ($#infilehandlesDist+1));
    my @infobest=((0) x ($#infilehandlesDist+1));
    my @formatbest=((0) x ($#infilehandlesDist+1));
    my @charbest=((0) x ($#infilehandlesDist+1));
    my @gtbest=((0) x ($#infilehandlesDist+1));
    my @PLtextbest=((0) x ($#infilehandlesDist+1));
    my @minorRef=((0) x ($#infilehandlesDist+1));
    #Correct PLs if alt alleles don't match - currently doesn't correct all combinations, but fixes most
    for (my $kk=0; $kk<=$#infilehandlesDist; $kk++) {
        if ($minchrompos==1020179655) {
            print "$chrompos[$kk],$chromposdups[$kk],$refdups[$kk],$refmax,$qual[$kk],$qualdups[$kk],$hapCall[$kk]\n";
        }
        if ($chrompos[$kk]!=$minchrompos || ($chromposdups[$kk]==$minchrompos && $refdups[$kk] eq $refmax && $qual[$kk]<=$qualdups[$kk] && $hapCall[$kk]==0)) {next;}
        if (!($ref[$kk] eq $refmax)) { $minorRef[$kk]=1; next;}
        
        my @alts = split( ",", $alt[$kk]);
        my @PLs = split( ",", $PLtext[$kk]);
        if (($PLtext[$kk] eq "")) {
            my @PLcorr= ((0) x $ALTPL[$altno-1]);
            $PLtext[$kk]=join(',', @PLcorr);
            next;
        }
        my @PLcorr= ((max(@PLs)) x $ALTPL[$altno-1]); #set values to max if for alt bases it doesn't have
        $PLcorr[0]=$PLs[0];
        if ($alts[0] eq $allalts[0]) {
            $PLcorr[1]=$PLs[1];
            $PLcorr[2]=$PLs[2];
            if ($#alts>0 && $alts[1] eq $allalts[1]) {
                $PLcorr[3]=$PLs[3];
                $PLcorr[4]=$PLs[4];
                $PLcorr[5]=$PLs[5];
                if ($#alts>1 && $alts[2] eq $allalts[2]) {
                    $PLcorr[6]=$PLs[6];
                    $PLcorr[7]=$PLs[7];
                    $PLcorr[8]=$PLs[8];
                    $PLcorr[9]=$PLs[9];
                    if ($#alts>2 && $alts[3] eq $allalts[3]) {
                        $PLcorr[10]=$PLs[10];
                        $PLcorr[11]=$PLs[11];
                        $PLcorr[12]=$PLs[12];
                        $PLcorr[13]=$PLs[13];
                        $PLcorr[14]=$PLs[14];
                    }
                }
            } elsif ($#alts>0 && $altno>2 && $alts[1] eq $allalts[2]) {
                $PLcorr[6]=$PLs[3];
                $PLcorr[7]=$PLs[4];
                $PLcorr[9]=$PLs[5];
            }
        } elsif ($altno>1 && $alts[0] eq $allalts[1]) {
            $PLcorr[3]=$PLs[1];
            $PLcorr[5]=$PLs[2];
            if ($#alts>0 && $altno>2 && $alts[1] eq $allalts[2]) {
                $PLcorr[6]=$PLs[3];
                $PLcorr[8]=$PLs[4];
                $PLcorr[9]=$PLs[5];
            }
        } elsif ($altno>2 && $alts[0] eq $allalts[2]) {
            $PLcorr[6]=$PLs[1];
            $PLcorr[9]=$PLs[2];
        } elsif ($altno>3 && $alts[0] eq $allalts[3]) {
            $PLcorr[10]=$PLs[1];
            $PLcorr[14]=$PLs[2];
        }
        
        if ($#alts>0 && $altno>1 && $alts[1] eq $allalts[1]) {
            $PLcorr[3]=$PLs[3];
            $PLcorr[5]=$PLs[5];
        } elsif ($#alts>0 && $altno>2 && $alts[1] eq $allalts[2]) {
            $PLcorr[6]=$PLs[3];
            $PLcorr[9]=$PLs[5];
        } elsif ($#alts>0 && $altno>3 && $alts[1] eq $allalts[3]) {
            $PLcorr[10]=$PLs[3];
            $PLcorr[14]=$PLs[5];
        }
        
        if ($#alts>1 && $altno>2 && $alts[2] eq $allalts[2]) {
            $PLcorr[6]=$PLs[6];
            $PLcorr[9]=$PLs[9];
        } elsif ($#alts>1 && $altno>3 && $alts[2] eq $allalts[3]) {
            $PLcorr[10]=$PLs[6];
            $PLcorr[14]=$PLs[9];
        }
        
        if ($#alts>2 && $altno>3 && $alts[3] eq $allalts[3]) {
            $PLcorr[10]=$PLs[10];
            $PLcorr[14]=$PLs[14];
        }
        $PLtextbest[$kk]=join(',', @PLcorr);
        
        $linebest[$kk] = $line[$kk];
        $chromposbest[$kk] = $chrompos[$kk];
        $chrombest[$kk] = $chrom[$kk];
        $posbest[$kk] = $pos[$kk];
        $idbest[$kk] = $id[$kk];
        $refbest[$kk] = $ref[$kk];
        $altbest[$kk] = $alt[$kk];
        $qualbest[$kk] = $qual[$kk];
        $filterbest[$kk] = $filter[$kk];
        $infobest[$kk] = $info[$kk];
        $formatbest[$kk] = $format[$kk];
        $charbest[$kk] = $char[$kk];
        $gtbest[$kk]=$gt[$kk];
    }

#    if ($minchrompos==1020179655) {
#        print "3\n";
#    }
    #Correct PLs if alt alleles don't match - currently doesn't correct all combinations, but fixes most
    for (my $kk=0; $kk<=$#infilehandlesDist; $kk++) {
        if (($chromposdups[$kk]!=$minchrompos || ($chrompos[$kk]==$minchrompos && ($hapCall[$kk]==1 || ($ref[$kk] eq $refmax && $qual[$kk]>$qualdups[$kk]))))) {next;}
        if (!($refdups[$kk] eq $refmax)) { $minorRef[$kk]=1; next;}
        
        my @alts = split( ",", $altdups[$kk]);
        my @PLs = split( ",", $PLtextdups[$kk]);
        if (($PLtextdups[$kk] eq "")) {
            my @PLcorr= ((0) x $ALTPL[$altno-1]);
            $PLtext[$kk]=join(',', @PLcorr);
            next;
        }
        my @PLcorr= ((max(@PLs)) x $ALTPL[$altno-1]); #set values to max if for alt bases it doesn't have
        $PLcorr[0]=$PLs[0];
        if ($alts[0] eq $allalts[0]) {
            $PLcorr[1]=$PLs[1];
            $PLcorr[2]=$PLs[2];
            if ($#alts>0 && $alts[1] eq $allalts[1]) {
                $PLcorr[3]=$PLs[3];
                $PLcorr[4]=$PLs[4];
                $PLcorr[5]=$PLs[5];
                if ($#alts>1 && $alts[2] eq $allalts[2]) {
                    $PLcorr[6]=$PLs[6];
                    $PLcorr[7]=$PLs[7];
                    $PLcorr[8]=$PLs[8];
                    $PLcorr[9]=$PLs[9];
                    if ($#alts>2 && $alts[3] eq $allalts[3]) {
                        $PLcorr[10]=$PLs[10];
                        $PLcorr[11]=$PLs[11];
                        $PLcorr[12]=$PLs[12];
                        $PLcorr[13]=$PLs[13];
                        $PLcorr[14]=$PLs[14];
                    }
                }
            } elsif ($#alts>0 && $altno>2 && $alts[1] eq $allalts[2]) {
                $PLcorr[6]=$PLs[3];
                $PLcorr[7]=$PLs[4];
                $PLcorr[9]=$PLs[5];
            }
        } elsif ($altno>1 && $alts[0] eq $allalts[1]) {
            $PLcorr[3]=$PLs[1];
            $PLcorr[5]=$PLs[2];
            if ($#alts>0 && $altno>2 && $alts[1] eq $allalts[2]) {
                $PLcorr[6]=$PLs[3];
                $PLcorr[8]=$PLs[4];
                $PLcorr[9]=$PLs[5];
            }
        } elsif ($altno>2 && $alts[0] eq $allalts[2]) {
            $PLcorr[6]=$PLs[1];
            $PLcorr[9]=$PLs[2];
        } elsif ($altno>3 && $alts[0] eq $allalts[3]) {
            $PLcorr[10]=$PLs[1];
            $PLcorr[14]=$PLs[2];
        }
        
        if ($#alts>0 && $altno>1 && $alts[1] eq $allalts[1]) {
            $PLcorr[3]=$PLs[3];
            $PLcorr[5]=$PLs[5];
        } elsif ($#alts>0 && $altno>2 && $alts[1] eq $allalts[2]) {
            $PLcorr[6]=$PLs[3];
            $PLcorr[9]=$PLs[5];
        } elsif ($#alts>0 && $altno>3 && $alts[1] eq $allalts[3]) {
            $PLcorr[10]=$PLs[3];
            $PLcorr[14]=$PLs[5];
        }
        
        if ($#alts>1 && $altno>2 && $alts[2] eq $allalts[2]) {
            $PLcorr[6]=$PLs[6];
            $PLcorr[9]=$PLs[9];
        } elsif ($#alts>1 && $altno>3 && $alts[2] eq $allalts[3]) {
            $PLcorr[10]=$PLs[6];
            $PLcorr[14]=$PLs[9];
        }
        
        if ($#alts>2 && $altno>3 && $alts[3] eq $allalts[3]) {
            $PLcorr[10]=$PLs[10];
            $PLcorr[14]=$PLs[14];
        }
        $PLtextbest[$kk]=join(',', @PLcorr);

        #use dup vcf instead of first vcf
        $linebest[$kk] = $linedups[$kk];
        $chromposbest[$kk] = $chromposdups[$kk];
        $chrombest[$kk] = $chromdups[$kk];
        $posbest[$kk] = $posdups[$kk];
        $idbest[$kk] = $iddups[$kk];
        $refbest[$kk] = $refdups[$kk];
        $altbest[$kk] = $altdups[$kk];
        $qualbest[$kk] = $qualdups[$kk];
        $filterbest[$kk] = $filterdups[$kk];
        $infobest[$kk] = $infodups[$kk];
        $formatbest[$kk] = $formatdups[$kk];
        $charbest[$kk] = $chardups[$kk];
        $gtbest[$kk]=$gtdups[$kk];

    }
    
#    if ($minchrompos==1020179655) {
#        print "4\n";
#    }
  
    my $HapNoVar=0; #number of datasets with haplotypecaller variant call within 20 bases but no variant at this position
    my $allPLannot=""; #PLs for each dataset at this position for allcall.vcf
    my @DP=((0) x ($#infilehandlesDist+1));
    my @PLminbest=((0) x ($#infilehandlesDist+1));
    
    for (my $i=0;$i<=$#infilehandlesDist;$i++) {
		if ($chrompos[$i]!=$minchrompos) {
            if ($chrompos[$i]<=$minchrompos+20) {
                if (!($qual[$i] eq ".") && $qual[$i]>40 && ($info[$i] =~ /UGHap=HapnoUG/ ||  $info[$i] =~ /UGHap=HapVarUGdiff/ || $info[$i] =~ /UGHap=Both/ || $info[$i] =~ /UGHap=both/)) { $HapNoVar++;}
            }
		} else {
        
            #if this dataset has a haplotypecaller variant call within 20 bases but no variant at this position, then note this
            if ($info[$i] =~ /UGHap=HapRefUGVar/) {
                $HapNoVar++;
                next;
            }

            if (!($PLtextbest[$i] eq "")) {
                $allPLannot="$allPLannot;PL${infiles[$i]}=$PLtextbest[$i]";
                my @PLs = split( ",", $PLtextbest[$i]);
                #$PLno=$#PLs;
                
                #print "$info\n";
                my $j=0; my $PL0=0; my $PLmin=1000000;
                foreach my $PL (@PLs) {
                    $PLTot[$j] += $PL;
                    if ($PL>20) { #confidently not this genotype in dataset
                        $NoPLTot[$j]+=1;
                    } elsif ($PL==0) { $PL0=$j; }
                    if ($PL>0 && $PL<$PLmin) { $PLmin=$PL; }
                    $j++;
                }
                if ($PLmin>20) { $YesPLTot[$PL0]++; }
                $PLminbest[$i]=$PLmin;
                
            
                if ($infobest[$i] =~ /DP=(.*?);/) {
                    $DPTot += $1;
                    $DP[$i]=$1;
                }
                if ($infobest[$i] =~ /MQ0Fraction=(.*?);/) {
                    if ($1>0.2) {$MQ0 += 1;}
                }
            }
        }
       
    }
    
    #last;
    my $geno=0; my $genoLowConf=0;
    
    #find genotype with lowest PL (highest likelihood) and next lowest
    my $PLmin1=1000000; my $PLmin2=1000000; my $NoPL=0; my $YesPL=0; my $PLTottxt="";
    for (my $j=0; $j<$PLno; $j++) {
        if ($PLTot[$j]<$PLmin1) {
            $geno=$j+1;
            $PLmin2=$PLmin1;
            $PLmin1=$PLTot[$j];
            $NoPL=$NoPLTot[$j];
            $YesPL=$YesPLTot[$j];
        } elsif ($PLTot[$j]<$PLmin2) {
            $PLmin2=$PLTot[$j];
        }
        #make PLTot into text format for output
        if ($j==0) {
            $PLTottxt=$PLTot[$j];
        } else {
            $PLTottxt="$PLTottxt,$PLTot[$j]";
        }
    }
    #find dataset with strongest genotype call to use its info in output
    my $maxdataset = 0; my $maxothPL=0;
    for (my $i=0;$i<=$#infilehandlesDist;$i++) {
		if ($chromposbest[$i]!=$minchrompos) {
            next;
		}
        if (!($PLtextbest[$i] eq "")) {
            my @PLs = split( ",", $PLtextbest[$i]);
            my $PLnoi=$#PLs;
            if ($PLnoi>=$geno-1 && $chromposbest[$i]==$minchrompos && $PLs[$geno-1]==0) {
                my $minPL=1000000;
                for (my $j=0; $j<=$PLnoi; $j++) {
                    if ($j==$geno-1) {
                        if ($PLs[$j]!=0) {$minPL=0; last;} #skip if a different genotype
                        next;
                    }
                    if ($PLs[$j]<$minPL) {
                        $minPL=$PLs[$j];
                    }
                }
                if ($minPL>$maxothPL) {
                    $maxdataset=$i;
                    $maxothPL=$minPL;
                }
            
            }
         }
            if ($minchrompos==2137104098) {
            	print "$maxdataset,$i,$infiles[$i]\n";
            	print "PLs:$PLtextbest[$i]\n";
            	print "$maxdataset,$maxothPL\n\n";
            }
            
    }
#    if ($minchrompos==1020179655) {
#        print "5\n";
#    }

    
    if ($HapNoVar/$#infilehandlesDist>0.5) {$geno=1;}
    my $PLmin=$PLmin2-$PLmin1; #likelihood of genotype
    if ($DPTot==0) {$DPTot=1;} #prevent divide by zero error
    my $PLbyDP=(int($PLmin/$DPTot*100+0.499)/100);
    $genoLowConf = $geno;
    if ($minchrompos>1020179655) {
        #exit;
    }

    if ($minchrompos==1226460240) {
        print "$chrombest[$maxdataset]\t$posbest[$maxdataset]\t$idbest[$maxdataset]\t$refmax\t$altVarLong\t$PLmin\tPASS\t$infobest[$maxdataset];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$geno-1]:$DP[$j]:$PLtextbest[$j]\n";
    }
    
    if (!defined($NoPL)) {
        print "$chrompos[$j],$chromposprev[$j],$minchrompos,$PLno\n";
        print "${infiles[$j]}\t$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$geno-1]:$DP[$j]:$PLtextbest[$j]\n";
        $xx++;
        if ($xx>10) {exit;}
    }
    if ((($geno==1 || ($HapNoVar<2 && ($geno==3 || $geno==6 || $geno==10 || $geno==15))) && $PLmin > 80 && $PLbyDP > 0.8 && $YesPL>1 && $NoPL/$#infilehandlesDist<0.25)) {
	    my $j = 0; #number of vcf file
        if ($chromposbest[$j]==1226460240) {
            print "Hom:@chrompos\n";
        }
		if ($geno==1) {
            print OUTHOMREF "$chrombest[$maxdataset]\t$posbest[$maxdataset]\t$idbest[$maxdataset]\t$refmax\t$altVarLong\t$PLmin\tPASS\t$infobest[$maxdataset];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]$allPLannot\tGT:DP:PL\t$GTs[$geno-1]:$DPTot:$PLTottxt\n";
#            foreach my $fh (@outfilehandlesHomRefUncert) {
#                if ($chromposbest[$j]==$chrompos[$maxdataset] && $refbest[$j] eq $refmax) {
#                    # Split up the line into an array
#                    @fields = split( "\t", $linebest[$j]);
#                    print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$geno-1]:$DP[$j]:$PLtextbest[$j]\n";
#                }
#                $j++;
#            }
		} else {
            print OUTHOMVAR "$chrom[$maxdataset]\t$pos[$maxdataset]\t$id[$maxdataset]\t$refmax\t$altVarLong\t$PLmin\tPASS\t$info[$maxdataset];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]$allPLannot\tGT:DP:PL\t$GTs[$geno-1]:$DPTot:$PLTottxt\n";
#            foreach my $fh (@outfilehandlesHomVarUncert) {
#                if ($chromposbest[$j]==$chrompos[$maxdataset] && $refbest[$j] eq $refmax) {
#                    # Split up the line into an array
#                    @fields = split( "\t", $linebest[$j]);
#                    print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$geno-1]:$DP[$j]:$PLtextbest[$j]\n";
#                }
#                $j++;
#            }
        }
        foreach my $fh (@outfilehandlesHomUncert) {
            if ($chromposbest[$j]==$chrompos[$maxdataset] && $refbest[$j] eq $refmax) {
                # Split up the line into an array
                @fields = split( "\t", $linebest[$j]);
                print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$geno-1]:$DP[$j]:$PLtextbest[$j]\n";
            }
            $j++;
        }
        
    } elsif (($geno!=1 && $geno!=3 && $geno!=6 && $geno!=10 && $geno!=15) && $MQ0/$#infilehandlesDist<0.3 && $PLmin > 100 && $PLmin/$DPTot > 3.4 && $YesPL>1 && $NoPL/$#infilehandlesDist<0.15 && $HapNoVar<2) {
        if ($chromposbest[$j]==1226460240) {
            print "Het:@chrompos\n";
        }
        my $j = 0; #number of vcf file
		foreach my $fh (@outfilehandlesHetRefUncert) {
            if ($chromposbest[$j]==$chrompos[$maxdataset] && $refbest[$j] eq $refmax) {
                # Split up the line into an array
                @fields = split( "\t", $linebest[$j]);
#                if ($chromposbest[$j]==1016359050) {
                if ($chrombest[$maxdataset] eq "0") {
                    print "@chrompos\n";
                    print "$chrompos[$maxdataset]:@chromposbest\n";
                    print "$refmax:@refbest; altvarlong:$altVarLong\n";
                    print "$chrompos[$j],$chromposprev[$j],$minchrompos,$PLno\n";
                    print "${infiles[$j]}\t$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$genoLowConf-1]:$DP[$j]:$PLtextbest[$j]\n";
                }
                print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$geno-1]:$DP[$j]:$PLtextbest[$j]\n";
            }
            
			$j++;
		}
                if (($chrombest[$maxdataset] eq "0")) {
                    print "@chrompos\n";
                    print "$chrompos[$maxdataset]:@chromposbest\n";
                    print "$refmax:@refbest; altvarlong:$altVarLong\n";
					print  "$chrombest[$maxdataset]\t$posbest[$maxdataset]\t$idbest[$maxdataset]\t$refmax\t$altVarLong\t$PLmin\tPASS\t$infobest[$maxdataset];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]$allPLannot\tGT:DP:PL\t$GTs[$geno-1]:$DPTot:$PLTottxt\n\n";
                }
		print OUTHET "$chrombest[$maxdataset]\t$posbest[$maxdataset]\t$idbest[$maxdataset]\t$refmax\t$altVarLong\t$PLmin\tPASS\t$infobest[$maxdataset];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]$allPLannot\tGT:DP:PL\t$GTs[$geno-1]:$DPTot:$PLTottxt\n";
        
    } else {
        if ($chromposbest[$j]==1226460240) {
            print "Uncert:@chrompos\n";
        }
        $genoLowConf = $geno;
        $geno = 0;
        my $j = 0; #number of vcf file
#        foreach my $fh (@outfilehandlesHomRefUncert) {
#            
#            if ($chromposbest[$j]==$chrompos[$maxdataset] && $refbest[$j] eq $refmax) {
#                # Split up the line into an array
#                @fields = split( "\t", $linebest[$j]);
#                print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;genoLowConf=$genoLowConf;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$genoLowConf-1]:$DP[$j]:$PLtextbest[$j]\n";
#            }
#            
#            $j++;
#        }
#        
#        $j = 0; #number of vcf file
#        foreach my $fh (@outfilehandlesHomVarUncert) {
#            
#            if ($chromposbest[$j]==$chrompos[$maxdataset] && $refbest[$j] eq $refmax) {
#                # Split up the line into an array
#                @fields = split( "\t", $linebest[$j]);
#                print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;genoLowConf=$genoLowConf;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$genoLowConf-1]:$DP[$j]:$PLtextbest[$j]\n";
#            }
#            
#            $j++;
#        }
#        
        $j = 0; #number of vcf file
        foreach my $fh (@outfilehandlesHetRefUncert) {
            
            if ($chromposbest[$j]==$chrompos[$maxdataset] && $refbest[$j] eq $refmax) {
                # Split up the line into an array
                @fields = split( "\t", $linebest[$j]);
                print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;genoLowConf=$genoLowConf;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$genoLowConf-1]:$DP[$j]:$PLtextbest[$j]\n";
            }
            $j++;
        }

        $j = 0; #number of vcf file
        foreach my $fh (@outfilehandlesHomUncert) {
            
            if ($chromposbest[$j]==$chrompos[$maxdataset] && $refbest[$j] eq $refmax) {
                # Split up the line into an array
                @fields = split( "\t", $linebest[$j]);
#                if ($chromposbest[$j]==1016359050) {
#                    print "@chrompos\n";
#                    print "$chrompos[$maxdataset]:@chromposbest\n";
#                    print "$refmax:@refbest\n";
#                    print "$chrompos[$j],$chromposprev[$j],$minchrompos,$PLno\n";
#                    print "${infiles[$j]}\t$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$genoLowConf-1]:$DP[$j]:$PLtextbest[$j]\n";
#                }
                print $fh "$fields[0]\t$fields[1]\t$fields[2]\t$refmax\t$altVarLong\t$fields[5]\tPASS\t$fields[7];geno=$geno;genoLowConf=$genoLowConf;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]\tGT:DP:PL\t$GTs[$genoLowConf-1]:$DP[$j]:$PLtextbest[$j]\n";
            }
            
            $j++;
        }
    }
    
    my $filtered="PASS";
    if ($geno==0) {$filtered="uncertain";}
    print OUTPUT "$chrombest[$maxdataset]\t$posbest[$maxdataset]\t$idbest[$maxdataset]\t$refmax\t$altVarLong\t$PLmin\tPASS\t$infobest[$maxdataset];geno=$geno;genoLowConf=$genoLowConf;PLtot=$PLmin;NoPLnum=$NoPL;DPTot=$DPTot;PLbyDP=$PLbyDP;HapNoVar=$HapNoVar;HapCallVar=$HapCallVar;LowMQ0Fraction=$MQ0;maxdataset=$infiles[$maxdataset]$allPLannot\tGT:DP:PL\t$GTs[$genoLowConf-1]:$DPTot:$PLTottxt\n";

	#output comma-delimited annotations for sites with data in all datasets
    my $dcnt=0;
    for (my $i=0; $i<=$#exp; $i++) { #first check that each dataset has info here
    	if ($chromposbest[$exp[$i]]==$chrompos[$maxdataset] && $refbest[$exp[$i]] eq $refmax && !($gtbest[$exp[$i]] eq "./.") && !($PLtextbest[$exp[$i]] eq "") && ($PLtextbest[$exp[$i]] =~ /.*,.*,/)) { $dcnt++; }
    }
    if ($altno==1 && $dcnt==$#exp+1) { #if all datasets have info, then export
			my $genoconf=0;
    		if ((($geno==1 || ($HapNoVar<2 && ($geno==3 || $geno==6 || $geno==10 || $geno==15))) && $PLmin > 80 && $PLbyDP > 0.8 && $YesPL>1 && $NoPL/$#infilehandlesDist<0.25)) {
    			if ($geno==1) { 
    				$genoconf=1;
    			} else {
    				$genoconf=3;
    			}
    		} elsif (($geno!=1 && $geno!=3 && $geno!=6 && $geno!=10 && $geno!=15) && $MQ0/$#infilehandlesDist<0.3 && $PLmin > 100 && $PLmin/$DPTot > 3.4 && $YesPL>1 && $NoPL/$#infilehandlesDist<0.15 && $HapNoVar<2) {
    			$genoconf=2;
    		}
		print EXPCSV "$chrompos[$maxdataset],$refmax,$altVarLong,$genoconf";    		

    	for (my $i=0; $i<=$#exp; $i++) {
            @fields = split( "\t", $linebest[$exp[$i]]);
			my $DP=0;
            if ($fields[7] =~ /DP=(.*?);/) { $DP=$1; }
			my $AB=1;
            if ($fields[7] =~ /ABH..=(.*?);/) { $AB=$1; }
			my $QD=0;
            if ($fields[7] =~ /QD=(.*?);/) { $QD=$1; }
			my $FS=0;
            if ($fields[7] =~ /FS=(.*?);/) { $FS=$1; }
			my $MQ=0;
            if ($fields[7] =~ /MQ=(.*?);/) { $MQ=$1; }
			my $NBQ=0;
            if ($fields[7] =~ /NBQ=(.*?);/) { $NBQ=$1; }
			my $ReadPosEndDist=0;
            if ($fields[7] =~ /ReadPosEndDist=(.*?);/) { $ReadPosEndDist=$1; }
			my $HaplotypeScore=0;
            if ($fields[7] =~ /HaplotypeScore=(.*?);/) { $HaplotypeScore=$1; }
			my $BaseQRankSum=0;
            if ($fields[7] =~ /BaseQRankSum=(.*?);/) { $BaseQRankSum=$1; }
			my $OND=0;
            if ($fields[7] =~ /OND=(.*?);/) { $OND=$1; }

            my $PLbyDP=0;
            if ($DP>0) { $PLbyDP=$PLminbest[$exp[$i]]/($DP+0.0); }
            my @PLs = split( ",", $PLtextbest[$exp[$i]]);
    		print EXPCSV ",$PLs[0],$PLs[1],$PLs[2],$PLbyDP,$DP,$FS,$MQ,$NBQ,$OND,$ReadPosEndDist,$HaplotypeScore,$BaseQRankSum";
		}
		print EXPCSV "\n";
	}    		
    
    #last;
}


close OUTPUT;
close OUTHOMREF;
close OUTHET;
close OUTHOMVAR;
close EXPCSV;
foreach my $fh (@infilehandlesDist) {
    close $fh;
}
#foreach my $fh (@outfilehandlesHomRefUncert) {
#    close $fh;
#}
foreach my $fh (@outfilehandlesHetRefUncert) {
    close $fh;
}
#foreach my $fh (@outfilehandlesHomVarUncert) {
#    close $fh;
#}
foreach my $fh (@outfilehandlesHomUncert) {
    close $fh;
}

exit;



sub readline {
    my $dataset = $_[0];
    $chromposprev[$dataset]=$chrompos[$dataset];
    $line[$dataset] = <$infilehandlesDist[$dataset]>;
    #print $line; next;
    if (defined($line[$dataset])) {
    } else { #quit loop if reach end of file
        $endline[$dataset] = 0;
        return;
    }
    
    # Split up the line into an array
    @fields = split( "\t", $line[$dataset]);
    my $chrno=$fields[0];
    if ($chrno eq "M" || $chrno eq "MT") {
        $chrno=23;
    } elsif ($chrno eq "X") {$chrno=24;
    } elsif ($chrno eq "Y") {$chrno=25;
    } elsif ($chrno =~ /^\d*$/) {
    } else { $endline[$dataset] = 0; return;}
    $chrompos[$dataset] = $chrno*1000000000+$fields[1];
    $chrom[$dataset] = $fields[0];
    $pos[$dataset] = $fields[1];
    $id[$dataset] = $fields[2];
    $ref[$dataset] = $fields[3];
    $alt[$dataset] = $fields[4];
    $qual[$dataset] = $fields[5];
    $filter[$dataset] = $fields[6];
    $info[$dataset] = $fields[7];
    $format[$dataset] = $fields[8];
    $char[$dataset] = $fields[9];

    my @formats = split( ":", $format[$dataset]);
    my @chars = split( ":", $char[$dataset]);
    
    for (my $i=0; $i<=$#formats; $i++) {
        if ($formats[$i] eq "GT") { $gt[$dataset]=$chars[$i]; }
        if ($formats[$i] eq "PL") { $PLtext[$dataset]=$chars[$i]; }
    }

    return;
    
}


