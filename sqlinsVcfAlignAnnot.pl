#!/usr/bin/perl -w

#

# sqlinsVcfAlignAnnot.pl - insert entire records with new alignment-related annotations to vcfannot table

#

#

# Version 1.0.0 - Jul 23, 2012

#



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


# Check if the right number of paramters are given at start up

#

# The following form is needed at start up:



if ($#ARGV != 2) {

 print "$#ARGV arguments\n usage: sqlinsVcfAlignAnnot.pl inputvcffilestart(without xxcall.vcf) bamAbbrev bamNum \n";

 exit;

}

my $infile = $ARGV[0];

my $bamAbb = $ARGV[1];

my $bam = $ARGV[2];



# Input file is opened 
my $fn;
$fn="${infile}${bamAbb}call_HomRef_AB_recal90.vcf";
unless ( open(INHOMRAB, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_HomRef_Align_recal90.vcf";
unless ( open(INHOMRALIGN, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_HomRef_map_recal90.vcf";
unless ( open(INHOMRMAP, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_HomVar_AB_recal90.vcf";
unless ( open(INHOMVAB, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_HomVar_Align_recal90.vcf";
unless ( open(INHOMVALIGN, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_HomVar_map_recal90.vcf";
unless ( open(INHOMVMAP, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_Het_SSE_recal90.vcf";
unless ( open(INHETRSSE, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_Het_AB_recal90.vcf";
unless ( open(INHETRAB, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_Het_Align_recal90.vcf";
unless ( open(INHETRALIGN, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}
$fn="${infile}${bamAbb}call_Het_map_recal90.vcf";
unless ( open(INHETRMAP, $fn)) { print "\nCannot open the file: $fn to read! \n\n"; exit;}


# Connect to database file
my $dbh = DBI->connect("dbi:SQLite:dbname=/Volumes/SSD960/workspace/data/NA12878vcfCompare_120315/NA12878vcfs.db", "", "", 
	{PrintError => 1, RaiseError =>1, AutoCommit => 0});


# Prepare INSERT statement
  my $sthInsAnnot = $dbh->prepare("INSERT INTO vcfannot (chrompos, bam, DP, FS, MQ, MQ0, MQ0frac, AB, BaseQRankSum, QD, ReadPosRankSum, MapQRankSum, OND, ABCI5, ABCI50, ABCI95, ReadPosEndDist, ReadMeanPos, HaplotypeScore, ReadMeanLen, MPGLikHetRef, MPGLikHomRef, MPGLikHomVar, VQSLODSSEHet, VQSLODABQDHet, VQSLODAlignHet, VQSLODMapHet, VQSLODSSEHomV, VQSLODABQDHomV, VQSLODAlignHomV, VQSLODMapHomV, VQSLODSSEHomR, VQSLODABQDHomR, VQSLODAlignHomR, VQSLODMapHomR, TrancheSSEHet, TrancheABQDHet, TrancheAlignHet, TrancheMapHet, TrancheSSEHomV, TrancheABQDHomV, TrancheAlignHomV, TrancheMapHomV, TrancheSSEHomR, TrancheABQDHomR, TrancheAlignHomR, TrancheMapHomR) VALUES (?,${bam},?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");

my $Log_file="${infile}${bamAbb}call_sqlupdvcfannot.log";
open (OUTLOG,">$Log_file") or die "Cannot open the log file submitted\n";
my @timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
my $i = 0;

my $endline=0;
my $headline=0;
my $incHomR=0;
my $incHomV=0;
my $incHet=0;
my $lineABQDHomR; my $lineAlignHomR; my $lineMapHomR;
my $lineSSEHet; my $lineABQDHet; my $lineAlignHet; my $lineMapHet;
my $lineABQDHomV; my $lineAlignHomV; my $lineMapHomV;

while ($endline==0) {
	if ($incHomR==0) { #skip if didn't use this line yet
		while (defined($lineABQDHomR = <INHOMRAB>) && $headline==0) {
			if ($lineABQDHomR =~ /^\#/) { next; } else {last;} #skip header lines
		}
		while (defined($lineAlignHomR = <INHOMRALIGN>) && $headline==0) {
			if ($lineAlignHomR =~ /^\#/) { next; } else {last;} #skip header lines
		}
		while (defined($lineMapHomR = <INHOMRMAP>) && $headline==0) {
			if ($lineMapHomR =~ /^\#/) { next; } else {last;} #skip header lines
		}
	}
	if ($incHomV==0) { #skip if didn't use this line yet
		while (defined($lineABQDHomV = <INHOMVAB>) && $headline==0) {
			if ($lineABQDHomV =~ /^\#/) { next; } else {last;} #skip header lines
		}
		while (defined($lineAlignHomV = <INHOMVALIGN>) && $headline==0) {
			if ($lineAlignHomV =~ /^\#/) { next; } else {last;} #skip header lines
		}
		while (defined($lineMapHomV = <INHOMVMAP>) && $headline==0) {
			if ($lineMapHomV =~ /^\#/) { next; } else {last;} #skip header lines
		}
	}
	if ($incHet==0) { #skip if didn't use this line yet
		while (defined($lineSSEHet = <INHETRSSE>) && $headline==0) {
			if ($lineSSEHet =~ /^\#/) { next; } else {last;} #skip header lines
		}
		while (defined($lineABQDHet = <INHETRAB>) && $headline==0) {
			if ($lineABQDHet =~ /^\#/) { next; } else {last;} #skip header lines
		}
		while (defined($lineAlignHet = <INHETRALIGN>) && $headline==0) {
			if ($lineAlignHet =~ /^\#/) { next; } else {last;} #skip header lines
		}
		while (defined($lineMapHet = <INHETRMAP>) && $headline==0) {
			if ($lineMapHet =~ /^\#/) { next; } else {last;} #skip header lines
		}
	}
	$headline=1;
    my @fields;
#last;

    # Split up the line into an array 

	if (defined($lineABQDHomR)) {} else {last;}
	@fields = split( "\t", $lineABQDHomR);
	if ($fields[0] =~ /^[0-9]/) {} else {$fields[0]=23;}
	my $chromposABQDHomR = 1000000000*$fields[0] + $fields[1];
	@fields = split( "\t", $lineABQDHomV);
	if ($fields[0] =~ /^[0-9]/) {} else {$fields[0]=23;}
	my $chromposABQDHomV = 1000000000*$fields[0] + $fields[1];
	@fields = split( "\t", $lineABQDHet);
	if ($fields[0] =~ /^[0-9]/) {} else {$fields[0]=23;}
	my $chromposABQDHet = 1000000000*$fields[0] + $fields[1];

	#check which vcfs have a row at this position
	$incHomR=1; $incHomV=1; $incHet=1;
	if ($chromposABQDHomR<=$chromposABQDHomV && $chromposABQDHomR<=$chromposABQDHet) {
		$incHomR=0;
		if ($chromposABQDHomR>=23000000000) {last;}
	}
	if ($chromposABQDHomV<=$chromposABQDHomR && $chromposABQDHomV<=$chromposABQDHet) {
		$incHomV=0;
		if ($chromposABQDHomV>=23000000000) {last;}
	}
	if ($chromposABQDHet<=$chromposABQDHomV && $chromposABQDHet<=$chromposABQDHomR) {
		$incHet=0;
		if ($chromposABQDHet>=23000000000) {last;}
	}
	
	my $TrancheSSEHomR=-1; my $TrancheABQDHomR=-1; my $TrancheAlignHomR=-1; my $TrancheMapHomR=-1; 
	my $VQSLODSSEHomR=-1; my $VQSLODABQDHomR=-1; my $VQSLODAlignHomR=-1; my $VQSLODMapHomR=-1; 
	my $TrancheSSEHomV=-1; my $TrancheABQDHomV=-1; my $TrancheAlignHomV=-1; my $TrancheMapHomV=-1; 
	my $VQSLODSSEHomV=-1; my $VQSLODABQDHomV=-1; my $VQSLODAlignHomV=-1; my $VQSLODMapHomV=-1; 
	my $TrancheSSEHet=-1; my $TrancheABQDHet=-1; my $TrancheAlignHet=-1; my $TrancheMapHet=-1; 
	my $VQSLODSSEHet=-1; my $VQSLODABQDHet=-1; my $VQSLODAlignHet=-1; my $VQSLODMapHet=-1; 
	if ($incHomR==0) { #skip if this line isn't in this vcf
		@fields = split( "\t", $lineABQDHomR);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheABQDHomR = $1;
		} else {$TrancheABQDHomR=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODABQDHomR = $1;
		} 
		@fields = split( "\t", $lineAlignHomR);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheAlignHomR = $1;
		} else {$TrancheAlignHomR=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODAlignHomR = $1;
		} 
		@fields = split( "\t", $lineMapHomR);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheMapHomR = $1;
		} else {$TrancheMapHomR=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODMapHomR = $1;
		} 
	}
	if ($incHomV==0) { #skip if this line isn't in this vcf
		@fields = split( "\t", $lineABQDHomV);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheABQDHomV = $1;
		} else {$TrancheABQDHomV=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODABQDHomV = $1;
		} 
		@fields = split( "\t", $lineAlignHomV);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheAlignHomV = $1;
		} else {$TrancheAlignHomV=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODAlignHomV = $1;
		} 
		@fields = split( "\t", $lineMapHomV);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheMapHomV = $1;
		} else {$TrancheMapHomV=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODMapHomV = $1;
		} 
	}
	if ($incHet==0) { #skip if this line isn't in this vcf
		@fields = split( "\t", $lineSSEHet);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheSSEHet = $1;
		} else {$TrancheSSEHet=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODSSEHet = $1;
		} 
		@fields = split( "\t", $lineABQDHet);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheABQDHet = $1;
		} else {$TrancheABQDHet=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODABQDHet = $1;
		} 
		@fields = split( "\t", $lineAlignHet);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheAlignHet = $1;
		} else {$TrancheAlignHet=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODAlignHet = $1;
		} 
		@fields = split( "\t", $lineMapHet);
		if ($fields[6] =~ /to(.*)/) { #find TruthSensitivityTranche
			$TrancheMapHet = $1;
		} else {$TrancheMapHet=0;} #PASS
		if ($fields[7] =~ /VQSLOD=(.*?);/) { #find VQSLOD
			$VQSLODMapHet = $1;
		} 
	}

	
	#use whatever the last @fields was since they're all the same for these annotations
	my $chrompos = 1000000000*$fields[0] + $fields[1];
	my $info = $fields[7]; #extract INFO field
	my $DP=0;
	if ($info =~ /DP=(.*?);/) { 
		$DP=$1;
	} 
	my $FS=0;
	if ($info =~ /FS=(.*?);/) { 
		$FS=$1;
	} 
	my $MQ=0;
	if ($info =~ /MQ=(.*?);/) { 
		$MQ=$1;
	} 
	my $MQ0=0;
	if ($info =~ /MQ0=(.*?);/) { 
		$MQ0=$1;
	} 
	my $MQ0frac=0;
	if ($info =~ /MQ0frac=(.*?);/) { 
		$MQ0frac=$1;
	} 
	my $AB=0;
	if ($info =~ /AB=(.*?);/) { 
		$AB=$1;
	} 
	my $BaseQRankSum=0;
	if ($info =~ /BaseQRankSum=(.*?);/) { 
		$BaseQRankSum=$1;
	} 
	my $QD=0;
	if ($info =~ /QD=(.*?);/) { 
		$QD=$1;
	} 
	my $ReadPosRankSum=0;
	if ($info =~ /ReadPosRankSum=(.*?);/) { 
		$ReadPosRankSum=$1;
	} 
	my $MapQRankSum=0;
	if ($info =~ /MapQRankSum=(.*?);/) { 
		$MapQRankSum=$1;
	} 
	my $OND=0;
	if ($info =~ /OND=(.*?);/) { 
		$OND=$1;
	} 
	my $ABCI5=0;
	if ($info =~ /ABCI5=(.*?);/) { 
		$ABCI5=$1;
	} 
	my $ABCI50=0;
	if ($info =~ /ABCI50=(.*?);/) { 
		$ABCI50=$1;
	} 
	my $ABCI95=0;
	if ($info =~ /ABCI95=(.*?);/) { 
		$ABCI95=$1;
	} 
	my $ReadPosEndDist=0;
	if ($info =~ /ReadPosEndDist=(.*?);/) { 
		$ReadPosEndDist=$1;
	} 
	my $ReadMeanPos=0;
	if ($info =~ /ReadMeanPos=(.*?);/) { 
		$ReadMeanPos=$1;
	} 
	my $HaplotypeScore=0;
	if ($info =~ /HaplotypeScore=(.*?);/) { 
		$HaplotypeScore=$1;
	} 
	my $ReadMeanLen=0;
	if ($info =~ /ReadMeanLen=(.*?);/) { 
		$ReadMeanLen=$1;
	} 
	my $MPGLikHetRef=0;
	if ($info =~ /MPGLikHetRef=(.*?);/) { 
		$MPGLikHetRef=$1;
	} 
	my $MPGLikHomRef=0;
	if ($info =~ /MPGLikHomRef=(.*?);/) { 
		$MPGLikHomRef=$1;
	} 
	my $MPGLikHomVar=0;
	if ($info =~ /MPGLikHomVar=(.*?);/) { 
		$MPGLikHomVar=$1;
	} 

	#chrompos, bam, DP, FS, MQ, MQ0, MQ0frac, AB, BaseQRankSum, QD, ReadPosRankSum, MapQRankSum, OND, ABCI5, ABCI50, ABCI95, ReadPosEndDist, ReadMeanPos, HaplotypeScore, ReadMeanLen, MPGLikHetRef, MPGLikHomRef, MPGLikHomVar, VQSLODSSEHet, VQSLODABQDHet, VQSLODAlignHet, VQSLODMapHet, VQSLODSSEHomV, VQSLODABQDHomV, VQSLODAlignHomV, VQSLODMapHomV, VQSLODSSEHomR, VQSLODABQDHomR, VQSLODAlignHomR, VQSLODMapHomR, TrancheSSEHet, TrancheABQDHet, TrancheAlignHet, TrancheMapHet, TrancheSSEHomV, TrancheABQDHomV, TrancheAlignHomV, TrancheMapHomV, TrancheSSEHomR, TrancheABQDHomR, TrancheAlignHomR, TrancheMapHomR
    $sthInsAnnot->execute($chrompos,$DP,$FS,$MQ,$MQ0,$MQ0frac,$AB,$BaseQRankSum,$QD,$ReadPosRankSum,$MapQRankSum,$OND,$ABCI5,$ABCI50,$ABCI95,$ReadPosEndDist,$ReadMeanPos,$HaplotypeScore,$ReadMeanLen,$MPGLikHetRef,$MPGLikHomRef,$MPGLikHomVar,$VQSLODSSEHet,$VQSLODABQDHet,$VQSLODAlignHet,$VQSLODMapHet,$VQSLODSSEHomV,$VQSLODABQDHomV,$VQSLODAlignHomV,$VQSLODMapHomV,$VQSLODSSEHomR,$VQSLODABQDHomR,$VQSLODAlignHomR,$VQSLODMapHomR,$TrancheSSEHet,$TrancheABQDHet,$TrancheAlignHet,$TrancheMapHet,$TrancheSSEHomV,$TrancheABQDHomV,$TrancheAlignHomV,$TrancheMapHomV,$TrancheSSEHomR,$TrancheABQDHomR,$TrancheAlignHomR,$TrancheMapHomR);
  if (($i % 40000) == 13) {
    $dbh->commit(); #commit inserts into database every 40000 records to save time
#    print $i/40000 . " ";
#	last;
  }
  $i++;
#  last;
#  print $info;
#  if ($i==4) {last;}
}
    $dbh->commit();
print OUTLOG "done:";
@timeData = localtime(time);
print OUTLOG join(' ', @timeData) . "\n";
close OUTLOG; 
close INHOMRAB;
close INHOMRALIGN;
close INHOMRMAP;
close INHOMVAB;
close INHOMVALIGN;
close INHOMVMAP;
close INHETRSSE;
close INHETRAB;
close INHETRALIGN;
close INHETRMAP;
$dbh->disconnect;



