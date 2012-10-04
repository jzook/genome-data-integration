#!/usr/bin/perl -w
#
# sqlupdateVcfUsingVQSR_v1.1.pl - arbitrate uncertain genotype calls using VQSR tranches for SSEs, AB/QD, Alignment, and Mapping artifacts
#
#
# Version 1.0.0 - Jul 24, 2012
# Version 1.1.0 - Jul 26, 2012 - don't use records with "Infinity" in any of the MPGLik fields (mostly SOLWG and 454WG)
#				Also add TrancheSSEmin, TrancheABQDmin, etc to report the minimum tranche that at least 2 datasets pass, which should be useful for assessing likelihood of genotype being correct
# Version 1.1.0 - Jul 26, 2012 - don't use records with "Infinity" in any of the MPGLik fields (mostly SOLWG and 454WG)
#
#Genotype categories:
# genoAll<1000 - called from combined datasets without filtering (genoAll/10 gives the number of datasets with MPG ratio of at least 2 in support of this genotype)
# genoxxx<1000 - called from datasets after filtering by SSE tranche (genoxxx/10 gives the number of datasets with MPG ratio of at least 2 in support of this genotype)
# genoxxx in 1000s - called from datasets after filtering by align tranche (genoxxx/10 gives the number of datasets with MPG ratio of at least 2 in support of this genotype)
# genoxxx in 2000s - called from datasets after filtering by number of bases in mapped read (genoxxx/10 gives the number of datasets with MPG ratio of at least 2 in support of this genotype)
#
print "test";
use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use FileHandle;
use English;
use File::Basename;
#use Try::Tiny;
use DBI;

my $line;
my $line2;

# Check if the right number of parameters are given at start up
#
# The following form is needed at start up:

if ($#ARGV != -1) {
 print "$#ARGV arguments\n usage: sqlupdateVcfUsingVQSR_v1.1.pl \n";
 exit;
}

# Connect to datRefCI50ase file
my $dbh = DBI->connect("dbi:SQLite:dbname=/Volumes/SSD960/workspace/data/NA12878vcfCompare_120315/NA12878vcfs.db", "", "",
	{PrintError => 1, RaiseError =>1, AutoCommit => 0});


# Prepare UPDATE statements to update genotypes and number of high quality mapped datasets
  my $sthUpdAll = $dbh->prepare("UPDATE chrposall SET genoAll=?,genoMapGood=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=?,genoSeq=?,genoMapGoodSeq=?,genoSeqV=?,genoMapGoodSeqV=?,TrancheSSEmin2V=?,TrancheABQDmin2V=?,TrancheAlignmin2V=?,TrancheMapmin2V=?,genoAllTranche=?,genoMapGoodAllTranche=?,genoSeqIndiv=?,genoMapGoodSeqIndiv=?,genoSeqABQD=?,genoMapGoodSeqABQD=?,genoAllTrancheABQD=?,genoMapGoodAllTrancheABQD=?,genoSeqIndivABQD=?,genoMapGoodSeqIndivABQD=? WHERE chrompos=?");
#  my $sthUpdGenoAll = $dbh->prepare("UPDATE chrposall SET genoAll=?,genoMapGood=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=? WHERE chrompos=?");
#  my $sthUpdGenoSeq = $dbh->prepare("UPDATE chrposall SET genoSeq=?,genoMapGoodSeq=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=? WHERE chrompos=?");
#  my $sthUpdGenoAllTranche = $dbh->prepare("UPDATE chrposall SET genoAllTranche=?,genoMapGoodAllTranche=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=? WHERE chrompos=?");
#  my $sthUpdGenoSeqIndiv = $dbh->prepare("UPDATE chrposall SET genoSeqIndiv=?,genoMapGoodSeqIndiv=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=? WHERE chrompos=?");
#  my $sthUpdGenoSeqABQD = $dbh->prepare("UPDATE chrposall SET genoSeqABQD=?,genoMapGoodSeqABQD=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=? WHERE chrompos=?");
#  my $sthUpdGenoAllTrancheABQD = $dbh->prepare("UPDATE chrposall SET genoAllTrancheABQD=?,genoMapGoodAllTrancheABQD=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=? WHERE chrompos=?");
#  my $sthUpdGenoSeqIndivABQD = $dbh->prepare("UPDATE chrposall SET genoSeqIndivABQD=?,genoMapGoodSeqIndivABQD=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=? WHERE chrompos=?");

  my $sthSelChrPos = $dbh->prepare("SELECT chrompos, bam, DP, ReadMeanLen, MPGLikHetRef, MPGLikHomRef, MPGLikHomVar, TrancheSSEHet, TrancheABQDHet, TrancheAlignHet, TrancheMapHet, TrancheSSEHomV, TrancheABQDHomV, TrancheAlignHomV, TrancheMapHomV, TrancheSSEHomR, TrancheABQDHomR, TrancheAlignHomR, TrancheMapHomR FROM vcfannot ORDER BY chrompos");



#Create arrays containing CIs for each annotation for each genotype for each bam file
my $bams = 11; #number of bam files

my @geno = ((0.0) x $bams);
my @DP = ((0.0) x $bams);
my @ReadMeanLen = ((0.0) x $bams);
my @MPGLikHetRef = ((0.0) x $bams);
my @MPGLikHomRef = ((0.0) x $bams);
my @MPGLikHomVar = ((0.0) x $bams);
my @TrancheSSE = ((0.0) x $bams);
my @TrancheABQD = ((0.0) x $bams);
my @TrancheAlign = ((0.0) x $bams);
my @TrancheMap = ((0.0) x $bams);


my $Log_file="sqlupdateVcfUsingVQSR.log";
open (OUTLOG,">$Log_file") or die "Cannot open the log file submitted\n";
my @timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
my $i = 0;



my $chrposprev = 0;
  my $chrpos;
  my $bam;
  $sthSelChrPos->execute();
my $jj=-1; #ignore datasets>$bams since one dataset is in twice
while (my $row = $sthSelChrPos->fetchrow_arrayref) {
    $chrpos = @$row[0];
   $bam = @$row[1];
    if ($chrpos == $chrposprev || $chrposprev==0) {
	$chrposprev = $chrpos;
    	$jj++;
    } else {
	my $MPGLikHetRefSum=0;
	my $MPGLikHomRefSum=0;
	my $MPGLikHomVarSum=0;
	my $DPSum=0;
	my $HomR=0;
	my $HomV=0;
	my $HetR=0;
	my $HomRmap=0;
	my $HomVmap=0;
	my $HetRmap=0;
	
	my $igenoAll=0;
	my $igenoMapGood=0;
	my $iTrancheSSEmin2=0;
	my $iTrancheABQDmin2=0;
	my $iTrancheAlignmin2=0;
	my $iTrancheMapmin2=0;
	my $igenoSeq=0;
	my $igenoMapGoodSeq=0;
	my $igenoSeqV=0;
	my $igenoMapGoodSeqV=0;
	my $iTrancheSSEmin2V=0;
	my $iTrancheABQDmin2V=0;
	my $iTrancheAlignmin2V=0;
	my $iTrancheMapmin2V=0;
	my $igenoAllTranche=0;
	my $igenoMapGoodAllTranche=0;
	my $igenoSeqIndiv=0;
	my $igenoMapGoodSeqIndiv=0;
	my $igenoSeqABQD=0;
	my $igenoMapGoodSeqABQD=0;
	my $igenoAllTrancheABQD=0;
	my $igenoMapGoodAllTrancheABQD=0;
	my $igenoSeqIndivABQD=0;
	my $igenoMapGoodSeqIndivABQD=0;


my $TrancheSSEmin1=100.0; my $TrancheSSEmin2=100.0; my $TrancheABQDmin1=100.0; my $TrancheABQDmin2=100.0; my $TrancheAlignmin1=100.0; my $TrancheAlignmin2=100.0; my $TrancheMapmin1=100.0; my $TrancheMapmin2=100.0;
        #first, check if there is a consensus genotype without filtering
        for (my $kk=0; $kk<=$jj; $kk++) {
                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
        }
        if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
        	$igenoAll=1+10*$HomR; $igenoMapGood=$HomRmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
        } elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
        	$igenoAll=3+10*$HomV; $igenoMapGood=$HomVmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
        } elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
        	$igenoAll=2+10*$HetR; $igenoMapGood=$HetRmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
        } else { #genotype cannot be determined from combined data, so try filtering methods
	        #First, see if 5x more datasets agree than disagree
			if ($igenoSeqV==0 && $DPSum>0 && $HomR>1 && $HomV+$HetR>0 && $HomR/($HomV+$HetR+0.0)>=5 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
				$igenoSeqV=15001+10*$HomR; $igenoMapGoodSeqV=$HomRmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
			} elsif ($igenoSeqV==0 && $DPSum>0 && $HomV>1 && $HomR+$HetR>0 && $HomV/($HomR+$HetR+0.0)>=5 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
				$igenoSeqV=15003+10*$HomV; $igenoMapGoodSeqV=$HomVmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
			} elsif ($igenoSeqV==0 && $DPSum>0 && $HetR>1 && $HomV+$HomR>0 && $HetR/($HomV+$HomR+0.0)>=5 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
				$igenoSeqV=15002+10*$HetR; $igenoMapGoodSeqV=$HetRmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
	        }
	        
	        my @tranches = (99.5, 99.0, 95.0, 90.0);
	        my @readlens = (40,60,80,90);

                #first filtering method starts with highest tranche, filters SSEs then alignment then readlength, and then tries lower tranches, but tests if datasets agree after each filter
	        my $tt=0;
                my $nexttranche=1;
	        while ($tt<4 && $nexttranche==1) {
                	$nexttranche=0;
	                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                for (my $kk=0; $kk<=$jj; $kk++) {
	                        if ($TrancheSSE[$kk]<=$tranches[$tt]) {
	                                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                        }
	                }
	                if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
				        	$igenoSeq=1+10*$HomR; $igenoMapGoodSeq=$HomRmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
	                } elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
				        	$igenoSeq=3+10*$HomV; $igenoMapGoodSeq=$HomVmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
	                } elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
				        	$igenoSeq=2+10*$HetR; $igenoMapGoodSeq=$HetRmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
	                } else { #genotype cannot be determined from filtering SSEs, so also filter alignment artifacts
							#First, see if 5x more datasets agree than disagree
							if ($igenoSeqV==0 && $DPSum>0 && $HomR>1 && $HomV+$HetR>0 && $HomR/($HomV+$HetR+0.0)>=5 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
								$igenoSeqV=10001+10*$HomR; $igenoMapGoodSeqV=$HomRmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
							} elsif ($igenoSeqV==0 && $DPSum>0 && $HomV>1 && $HomR+$HetR>0 && $HomV/($HomR+$HetR+0.0)>=5 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
								$igenoSeqV=10003+10*$HomV; $igenoMapGoodSeqV=$HomVmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
							} elsif ($igenoSeqV==0 && $DPSum>0 && $HetR>1 && $HomV+$HomR>0 && $HetR/($HomV+$HomR+0.0)>=5 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
								$igenoSeqV=10002+10*$HetR; $igenoMapGoodSeqV=$HetRmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
							}
	                        $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                        $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                        for (my $kk=0; $kk<=$jj; $kk++) {
	                                if ($TrancheSSE[$kk]<=$tranches[$tt] && $TrancheAlign[$kk]<=$tranches[$tt]) {
	                                        $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                        $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                        $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                        $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                }
	                        }
	                        if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
									$igenoSeq=1001+10*$HomR; $igenoMapGoodSeq=$HomRmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
							} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
									$igenoSeq=1003+10*$HomV; $igenoMapGoodSeq=$HomVmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
							} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
									$igenoSeq=1002+10*$HetR; $igenoMapGoodSeq=$HetRmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
	                        } else { #genotype cannot be determined from filtering SSEs and alignment artifacts, so also filter short reads
									#First, see if 5x more datasets agree than disagree
									if ($igenoSeqV==0 && $DPSum>0 && $HomR>1 && $HomV+$HetR>0 && $HomR/($HomV+$HetR+0.0)>=5 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
										$igenoSeqV=11001+10*$HomR; $igenoMapGoodSeqV=$HomRmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
									} elsif ($igenoSeqV==0 && $DPSum>0 && $HomV>1 && $HomR+$HetR>0 && $HomV/($HomR+$HetR+0.0)>=5 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
										$igenoSeqV=11003+10*$HomV; $igenoMapGoodSeqV=$HomVmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
									} elsif ($igenoSeqV==0 && $DPSum>0 && $HetR>1 && $HomV+$HomR>0 && $HetR/($HomV+$HomR+0.0)>=5 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
										$igenoSeqV=11002+10*$HetR; $igenoMapGoodSeqV=$HetRmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
									}
	                                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                                for (my $kk=0; $kk<=$jj; $kk++) {
	                                        if ($TrancheSSE[$kk]<=$tranches[$tt] && $TrancheAlign[$kk]<=$tranches[$tt] && $ReadMeanLen[$kk]>$readlens[$tt]) {
	                                                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                        }
	                                }
                                                
	                                if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
											$igenoSeq=2001+10*$HomR; $igenoMapGoodSeq=$HomRmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
									} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
											$igenoSeq=2003+10*$HomV; $igenoMapGoodSeq=$HomVmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
									} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
											$igenoSeq=2002+10*$HetR; $igenoMapGoodSeq=$HetRmap;$iTrancheSSEmin2=$TrancheSSEmin2; $iTrancheABQDmin2=$TrancheABQDmin2; $iTrancheAlignmin2=$TrancheAlignmin2; $iTrancheMapmin2=$TrancheMapmin2;
	                                } else { #genotype cannot be determined from filtering this tranche, so try the next one
											#First, see if 5x more datasets agree than disagree
											if ($igenoSeqV==0 && $DPSum>0 && $HomR>1 && $HomV+$HetR>0 && $HomR/($HomV+$HetR+0.0)>=5 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
												$igenoSeqV=12001+10*$HomR; $igenoMapGoodSeqV=$HomRmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
											} elsif ($igenoSeqV==0 && $DPSum>0 && $HomV>1 && $HomR+$HetR>0 && $HomV/($HomR+$HetR+0.0)>=5 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
												$igenoSeqV=12003+10*$HomV; $igenoMapGoodSeqV=$HomVmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
											} elsif ($igenoSeqV==0 && $DPSum>0 && $HetR>1 && $HomV+$HomR>0 && $HetR/($HomV+$HomR+0.0)>=5 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
												$igenoSeqV=12002+10*$HetR; $igenoMapGoodSeqV=$HetRmap;$iTrancheSSEmin2V=$TrancheSSEmin2; $iTrancheABQDmin2V=$TrancheABQDmin2; $iTrancheAlignmin2V=$TrancheAlignmin2; $iTrancheMapmin2V=$TrancheMapmin2;
											}
	                                        $nexttranche=1;
	                                }
	                        }
	                }
	                $tt++;
                }

                #second filtering method  filters SSEs, alignment, and readlength simultaneously for each tranche, and then tries lower tranches, but tests if datasets agree after each filter
	        $tt=0;
                $nexttranche=1;
	        while ($tt<4 && $nexttranche==1) {
                	$nexttranche=0;
	                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                for (my $kk=0; $kk<=$jj; $kk++) {
                                if ($TrancheSSE[$kk]<=$tranches[$tt] && $TrancheAlign[$kk]<=$tranches[$tt] && $ReadMeanLen[$kk]>$readlens[$tt]) {
                                        $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
                                        $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
                                        $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
                                        $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
                                }
                        }
                         
                        if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
								$igenoAllTranche=1+10*$HomR; $igenoMapGoodAllTranche=$HomRmap;
						} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
								$igenoAllTranche=3+10*$HomV; $igenoMapGoodAllTranche=$HomVmap;
						} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
								$igenoAllTranche=2+10*$HetR; $igenoMapGoodAllTranche=$HetRmap;
                        } else { #genotype cannot be determined from filtering this tranche, so try the next one
                                $nexttranche=1;
                        }
	                $tt++;
                }

                #third filtering method starts with highest tranche, filters only SSEs then only alignment then only readlength, and then tries lower tranches, but tests if datasets agree after each filter
	        $tt=0;
                $nexttranche=1;
	        while ($tt<4 && $nexttranche==1) {
                	$nexttranche=0;
	                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                for (my $kk=0; $kk<=$jj; $kk++) {
	                        if ($TrancheSSE[$kk]<=$tranches[$tt]) {
	                                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                        }
	                }
	                if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
							$igenoSeqIndiv=1+10*$HomR; $igenoMapGoodSeqIndiv=$HomRmap;
					} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
							$igenoSeqIndiv=3+10*$HomV; $igenoMapGoodSeqIndiv=$HomVmap;
					} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
							$igenoSeqIndiv=2+10*$HetR; $igenoMapGoodSeqIndiv=$HetRmap;
	                } else { #genotype cannot be determined from filtering SSEs, so instead filter alignment artifacts
	                        $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                        $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                        for (my $kk=0; $kk<=$jj; $kk++) {
	                                if ($TrancheAlign[$kk]<=$tranches[$tt]) {
	                                        $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                        $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                        $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                        $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                }
	                        }
	                        if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
									$igenoSeqIndiv=1001+10*$HomR; $igenoMapGoodSeqIndiv=$HomRmap;
							} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
									$igenoSeqIndiv=1003+10*$HomV; $igenoMapGoodSeqIndiv=$HomVmap;
							} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
									$igenoSeqIndiv=1002+10*$HetR; $igenoMapGoodSeqIndiv=$HetRmap;
	                        } else { #genotype cannot be determined from filtering SSEs or alignment artifacts, so instead filter short reads
	                                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                                for (my $kk=0; $kk<=$jj; $kk++) {
	                                        if ($ReadMeanLen[$kk]>$readlens[$tt]) {
	                                                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                        }
	                                }
                                                
	                                if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
											$igenoSeqIndiv=2001+10*$HomR; $igenoMapGoodSeqIndiv=$HomRmap;
									} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
											$igenoSeqIndiv=2003+10*$HomV; $igenoMapGoodSeqIndiv=$HomVmap;
									} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
											$igenoSeqIndiv=2002+10*$HetR; $igenoMapGoodSeqIndiv=$HetRmap;
	                                } else { #genotype cannot be determined from filtering this tranche, so try the next one
	                                        $nexttranche=1;
	                                }
	                        }
	                }
	                $tt++;
                }


                #first filtering method starts with highest tranche, filters SSEs then alignment then readlength, and then tries lower tranches, but tests if datasets agree after each filter
	        $tt=0;
                $nexttranche=1;
	        while ($tt<4 && $nexttranche==1) {
                	$nexttranche=0;
	                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                for (my $kk=0; $kk<=$jj; $kk++) {
	                        if ($TrancheSSE[$kk]<=$tranches[$tt]) {
	                                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                        }
	                }
	                if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
				        	$igenoSeqABQD=1+10*$HomR; $igenoMapGoodSeqABQD=$HomRmap;
	                } elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
				        	$igenoSeqABQD=3+10*$HomV; $igenoMapGoodSeqABQD=$HomVmap;
	                } elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
				        	$igenoSeqABQD=2+10*$HetR; $igenoMapGoodSeqABQD=$HetRmap;
	                } else { #genotype cannot be determined from filtering SSEs, so also filter alignment artifacts
	                        $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                        $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                        for (my $kk=0; $kk<=$jj; $kk++) {
	                                if ($TrancheSSE[$kk]<=$tranches[$tt] && $TrancheAlign[$kk]<=$tranches[$tt]) {
	                                        $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                        $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                        $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                        $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                }
	                        }
	                        if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
									$igenoSeqABQD=1001+10*$HomR; $igenoMapGoodSeqABQD=$HomRmap;
							} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
									$igenoSeqABQD=1003+10*$HomV; $igenoMapGoodSeqABQD=$HomVmap;
							} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
									$igenoSeqABQD=1002+10*$HetR; $igenoMapGoodSeqABQD=$HetRmap;
	                        } else { #genotype cannot be determined from filtering SSEs and alignment artifacts, so also filter short reads
	                                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                                for (my $kk=0; $kk<=$jj; $kk++) {
	                                        if ($TrancheSSE[$kk]<=$tranches[$tt] && $TrancheAlign[$kk]<=$tranches[$tt] && $ReadMeanLen[$kk]>$readlens[$tt]) {
	                                                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                        }
	                                }
	                                if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
											$igenoSeqABQD=2001+10*$HomR; $igenoMapGoodSeqABQD=$HomRmap;
									} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
											$igenoSeqABQD=2003+10*$HomV; $igenoMapGoodSeqABQD=$HomVmap;
									} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
											$igenoSeqABQD=2002+10*$HetR; $igenoMapGoodSeqABQD=$HetRmap;
	                                } else { #genotype cannot be determined from filtering SSEs, alignment artifacts, and short reads, so also filter AB/QD
	                                        $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                                        $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                                        for (my $kk=0; $kk<=$jj; $kk++) {
	                                                if ($TrancheSSE[$kk]<=$tranches[$tt] && $TrancheAlign[$kk]<=$tranches[$tt] && $ReadMeanLen[$kk]>$readlens[$tt] && $TrancheABQD[$kk]<=$tranches[$tt]) {
	                                                        $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                                        $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                                        $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                                        $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                                }
	                                        }
                                                
	                                        if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
													$igenoSeqABQD=3001+10*$HomR; $igenoMapGoodSeqABQD=$HomRmap;
											} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
													$igenoSeqABQD=3003+10*$HomV; $igenoMapGoodSeqABQD=$HomVmap;
											} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
													$igenoSeqABQD=3002+10*$HetR; $igenoMapGoodSeqABQD=$HetRmap;
	                                        } else { #genotype cannot be determined from filtering this tranche, so try the next one
	                                                $nexttranche=1;
	                                        }
                                        }
	                        }
	                }
	                $tt++;
                }

                #second filtering method  filters SSEs, alignment, and readlength simultaneously for each tranche, and then tries lower tranches, but tests if datasets agree after each filter
	        $tt=0;
                $nexttranche=1;
	        while ($tt<4 && $nexttranche==1) {
                	$nexttranche=0;
	                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                for (my $kk=0; $kk<=$jj; $kk++) {
                                if ($TrancheSSE[$kk]<=$tranches[$tt] && $TrancheAlign[$kk]<=$tranches[$tt] && $ReadMeanLen[$kk]>$readlens[$tt] && $TrancheABQD[$kk]<=$tranches[$tt]) {
                                        $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
                                        $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
                                        $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
                                        $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
                                }
                        }
                        
                        if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
								$igenoAllTrancheABQD=1+10*$HomR; $igenoMapGoodAllTrancheABQD=$HomRmap;
						} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
								$igenoAllTrancheABQD=3+10*$HomV; $igenoMapGoodAllTrancheABQD=$HomVmap;
						} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
								$igenoAllTrancheABQD=2+10*$HetR; $igenoMapGoodAllTrancheABQD=$HetRmap;
                        } else { #genotype cannot be determined from filtering this tranche, so try the next one
                                $nexttranche=1;
                        }
	                $tt++;
                }

                #third filtering method starts with highest tranche, filters only SSEs then only alignment then only readlength, and then tries lower tranches, but tests if datasets agree after each filter
	        $tt=0;
                $nexttranche=1;
	        while ($tt<4 && $nexttranche==1) {
                	$nexttranche=0;
	                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                for (my $kk=0; $kk<=$jj; $kk++) {
	                        if ($TrancheSSE[$kk]<=$tranches[$tt]) {
	                                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                        }
	                }
	                if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
							$igenoSeqIndivABQD=1+10*$HomR; $igenoMapGoodSeqIndivABQD=$HomRmap;
					} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
							$igenoSeqIndivABQD=3+10*$HomV; $igenoMapGoodSeqIndivABQD=$HomVmap;
					} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
							$igenoSeqIndivABQD=2+10*$HetR; $igenoMapGoodSeqIndivABQD=$HetRmap;
	                } else { #genotype cannot be determined from filtering SSEs, so instead filter alignment artifacts
	                        $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                        $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                        for (my $kk=0; $kk<=$jj; $kk++) {
	                                if ($TrancheAlign[$kk]<=$tranches[$tt]) {
	                                        $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                        $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                        $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                        $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                }
	                        }
	                        if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
									$igenoSeqIndivABQD=1001+10*$HomR; $igenoMapGoodSeqIndivABQD=$HomRmap;
							} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
									$igenoSeqIndivABQD=1003+10*$HomV; $igenoMapGoodSeqIndivABQD=$HomVmap;
							} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
									$igenoSeqIndivABQD=1002+10*$HetR; $igenoMapGoodSeqIndivABQD=$HetRmap;
	                        } else { #genotype cannot be determined from filtering SSEs or alignment artifacts, so instead filter short reads
	                                $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                                $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0;
	                                for (my $kk=0; $kk<=$jj; $kk++) {
	                                        if ($ReadMeanLen[$kk]>$readlens[$tt]) {
	                                                $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                                $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                                $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                                $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                        }
	                                }
	                                if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
											$igenoSeqIndivABQD=2001+10*$HomR; $igenoMapGoodSeqIndivABQD=$HomRmap;
									} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
											$igenoSeqIndivABQD=2003+10*$HomV; $igenoMapGoodSeqIndivABQD=$HomVmap;
									} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
											$igenoSeqIndivABQD=2002+10*$HetR; $igenoMapGoodSeqIndivABQD=$HetRmap;
	                                } else { #genotype cannot be determined from filtering SSEs, alignment artifacts, or short reads, so instead filter AB/QD
	                                        $MPGLikHetRefSum=0; $MPGLikHomRefSum=0; $MPGLikHomVarSum=0; $DPSum=0;
	                                        $HomR=0; $HomV=0; $HetR=0; $HomRmap=0; $HomVmap=0; $HetRmap=0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0; $TrancheSSEmin1=100.0; $TrancheSSEmin2=100.0; $TrancheABQDmin1=100.0; $TrancheABQDmin2=100.0; $TrancheAlignmin1=100.0; $TrancheAlignmin2=100.0; $TrancheMapmin1=100.0; $TrancheMapmin2=100.0; 
	                                        for (my $kk=0; $kk<=$jj; $kk++) {
	                                                if ($TrancheABQD[$kk]<=$tranches[$tt]) {
	                                                        $MPGLikHetRefSum+=$MPGLikHetRef[$kk];
	                                                        $MPGLikHomRefSum+=$MPGLikHomRef[$kk];
	                                                        $MPGLikHomVarSum+=$MPGLikHomVar[$kk];
	                                                        $DPSum+=$DP[$kk];
	                                                        my $uncert=0; 
	                                                        if ($geno[$kk]==1 && ($MPGLikHetRef[$kk]-$MPGLikHomRef[$kk])>2) {
	                                                                $HomR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomRmap++; }
	                                                        } elsif ($geno[$kk]==3 && ($MPGLikHetRef[$kk]-$MPGLikHomVar[$kk])>2) {
	                                                                $HomV++;
	                                                                if ($TrancheMap[$kk]<=95) { $HomVmap++; }
	                                                        } elsif ($geno[$kk]==2 && ($MPGLikHomRef[$kk]-$MPGLikHetRef[$kk])>2 && ($MPGLikHomVar[$kk]-$MPGLikHetRef[$kk])>2) {
	                                                                $HetR++;
	                                                                if ($TrancheMap[$kk]<=95) { $HetRmap++; }
	                                                        } else {$uncert=1;}
	                                                        if ($uncert==0) {
	                                                        if ($TrancheSSE[$kk]<=$TrancheSSEmin1 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSEmin1; 
	                                                                $TrancheSSEmin1=$TrancheSSE[$kk];
	                                                        } elsif ($TrancheSSE[$kk]<$TrancheSSEmin2 && $TrancheSSE[$kk]>=0) {
	                                                                $TrancheSSEmin2=$TrancheSSE[$kk];
	                                                        } 
	                                                        if ($TrancheABQD[$kk]<=$TrancheABQDmin1 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQDmin1; 
	                                                                $TrancheABQDmin1=$TrancheABQD[$kk];
	                                                        } elsif ($TrancheABQD[$kk]<$TrancheABQDmin2 && $TrancheABQD[$kk]>=0) {
	                                                                $TrancheABQDmin2=$TrancheABQD[$kk];
	                                                        } 
	                                                        if ($TrancheAlign[$kk]<=$TrancheAlignmin1 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlignmin1; 
	                                                                $TrancheAlignmin1=$TrancheAlign[$kk];
	                                                        } elsif ($TrancheAlign[$kk]<$TrancheAlignmin2 && $TrancheAlign[$kk]>=0) {
	                                                                $TrancheAlignmin2=$TrancheAlign[$kk];
	                                                        } 
	                                                        if ($TrancheMap[$kk]<=$TrancheMapmin1 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMapmin1; 
	                                                                $TrancheMapmin1=$TrancheMap[$kk];
	                                                        } elsif ($TrancheMap[$kk]<$TrancheMapmin2 && $TrancheMap[$kk]>=0) {
	                                                                $TrancheMapmin2=$TrancheMap[$kk];
	                                                        } 
	                                                        } 
	                                                }
	                                        }
                                                
                                            if ($DPSum>0 && $HomR>1 && $HomV+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomRefSum>12 && ($MPGLikHetRefSum-$MPGLikHomRefSum)/($DPSum+0.0)>0.16) {
													$igenoSeqIndivABQD=3001+10*$HomR; $igenoMapGoodSeqIndivABQD=$HomRmap;
											} elsif ($DPSum>0 && $HomV>1 && $HomR+$HetR==0 && $MPGLikHetRefSum-$MPGLikHomVarSum>12 && ($MPGLikHetRefSum-$MPGLikHomVarSum)/($DPSum+0.0)>0.16) {
													$igenoSeqIndivABQD=3003+10*$HomV; $igenoMapGoodSeqIndivABQD=$HomVmap;
											} elsif ($DPSum>0 && $HetR>1 && $HomV+$HomR==0 && $MPGLikHomRefSum-$MPGLikHetRefSum>20 && ($MPGLikHomRefSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68  && $MPGLikHomVarSum-$MPGLikHetRefSum>20 && ($MPGLikHomVarSum-$MPGLikHetRefSum)/($DPSum+0.0)>0.68) {
													$igenoSeqIndivABQD=3002+10*$HetR; $igenoMapGoodSeqIndivABQD=$HetRmap;
	                                        } else { #genotype cannot be determined from filtering this tranche, so try the next one
	                                                $nexttranche=1;
	                                        }
                                        }
	                        }
	                }
	                $tt++;
                }

        }
#  my $sthUpdAll = $dbh->prepare("UPDATE chrposall SET genoAll=?,genoMapGood=?,TrancheSSEmin2=?,TrancheABQDmin2=?,TrancheAlignmin2=?,TrancheMapmin2=?,genoSeq=?,genoMapGoodSeq=?,genoSeqV=?,genoMapGoodSeqV=?,TrancheSSEmin2V=?,TrancheABQDmin2V=?,TrancheAlignmin2V=?,TrancheMapmin2V=?,genoAllTranche=?,genoMapGoodAllTranche=?,genoSeqIndiv=?,genoMapGoodSeqIndiv=?,genoSeqABQD=?,genoMapGoodSeqABQD=?,genoAllTrancheABQD=?,genoMapGoodAllTrancheABQD=?,genoSeqIndivABQD=?,genoMapGoodSeqIndivABQD=? WHERE chrompos=?");
    if ($igenoSeqV==0) {
    	$igenoSeqV=$igenoSeq;
    	$igenoMapGoodSeqV=$igenoMapGoodSeq;
    }
#    if ($igenoSeqV>15000) {    	$sthUpdAll->execute($igenoAll,$igenoMapGood,$iTrancheSSEmin2,$iTrancheABQDmin2,$iTrancheAlignmin2,$iTrancheMapmin2,$igenoSeq,$igenoMapGoodSeq,$igenoSeqV,$igenoMapGoodSeqV,$iTrancheSSEmin2V,$iTrancheABQDmin2V,$iTrancheAlignmin2V,$iTrancheMapmin2V,$igenoAllTranche,$igenoMapGoodAllTranche,$igenoSeqIndiv,$igenoMapGoodSeqIndiv,$igenoSeqABQD,$igenoMapGoodSeqABQD,$igenoAllTrancheABQD,$igenoMapGoodAllTrancheABQD,$igenoSeqIndivABQD,$igenoMapGoodSeqIndivABQD,$chrposprev);}
      	$sthUpdAll->execute($igenoAll,$igenoMapGood,$iTrancheSSEmin2,$iTrancheABQDmin2,$iTrancheAlignmin2,$iTrancheMapmin2,$igenoSeq,$igenoMapGoodSeq,$igenoSeqV,$igenoMapGoodSeqV,$iTrancheSSEmin2V,$iTrancheABQDmin2V,$iTrancheAlignmin2V,$iTrancheMapmin2V,$igenoAllTranche,$igenoMapGoodAllTranche,$igenoSeqIndiv,$igenoMapGoodSeqIndiv,$igenoSeqABQD,$igenoMapGoodSeqABQD,$igenoAllTrancheABQD,$igenoMapGoodAllTrancheABQD,$igenoSeqIndivABQD,$igenoMapGoodSeqIndivABQD,$chrposprev);

   	@geno = ((0.0) x $bams);
	@DP = ((0.0) x $bams);
	@ReadMeanLen = ((0.0) x $bams);
	@MPGLikHetRef = ((0.0) x $bams);
	@MPGLikHomRef = ((0.0) x $bams);
	@MPGLikHomVar = ((0.0) x $bams);
	@TrancheSSE = ((0.0) x $bams);
	@TrancheABQD = ((0.0) x $bams);
	@TrancheAlign = ((0.0) x $bams);
	@TrancheMap = ((0.0) x $bams);
	$chrposprev = $chrpos;
	$jj=0;
	$i++;
	  if (($i % 10000) == 10) {
	    $dbh->commit(); #commit inserts into database every 10000 records to save time
           # last; #DEBUG
	   # print $i . " $chrposprev\n";
	  }
    }
  # chrompos, bam, DP, ReadMeanLen, MPGLikHetRef, MPGLikHomRef, MPGLikHomVar, TrancheSSEHet, TrancheABQDHet, TrancheAlignHet, TrancheMapHet, TrancheSSEHomV, TrancheABQDHomV, TrancheAlignHomV, TrancheMapHomV, TrancheSSEHomR, TrancheABQDHomR, TrancheAlignHomR, TrancheMapHomR
	if (@$row[4] =~ /Inf/ || @$row[5] =~ /Inf/ || @$row[6] =~ /Inf/) {next;} #skip this record if any of the MPGLik fields have Infinity
	$DP[$jj] = @$row[2];
	$ReadMeanLen[$jj] = @$row[3];
	$MPGLikHetRef[$jj] = @$row[4];
	$MPGLikHomRef[$jj] = @$row[5];
	$MPGLikHomVar[$jj] = @$row[6];
	if ($MPGLikHetRef[$jj]<$MPGLikHomRef[$jj] && $MPGLikHetRef[$jj]<$MPGLikHomVar[$jj]) {
        	$TrancheSSE[$jj] = @$row[7];
        	$TrancheABQD[$jj] = @$row[8];
		$TrancheAlign[$jj] = @$row[9];
		$TrancheMap[$jj] = @$row[10];
                $geno[$jj] = 2;
        } elsif ($MPGLikHomVar[$jj]<$MPGLikHetRef[$jj] && $MPGLikHomVar[$jj]<$MPGLikHomRef[$jj]) {
        	$TrancheSSE[$jj] = @$row[11];
        	$TrancheABQD[$jj] = @$row[12];
		$TrancheAlign[$jj] = @$row[13];
		$TrancheMap[$jj] = @$row[14];
                $geno[$jj] = 3;
        } elsif ($MPGLikHomRef[$jj]<$MPGLikHetRef[$jj] && $MPGLikHomRef[$jj]<$MPGLikHomVar[$jj]) {
        	$TrancheSSE[$jj] = @$row[15];
        	$TrancheABQD[$jj] = @$row[16];
		$TrancheAlign[$jj] = @$row[17];
		$TrancheMap[$jj] = @$row[18];
                $geno[$jj] = 1;
	} else {
                $geno[$jj] = 0;
	}

#  last;
}

#run updates for last record


    $dbh->commit();
print OUTLOG "done:";
@timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
close OUTLOG;

$dbh->disconnect;


