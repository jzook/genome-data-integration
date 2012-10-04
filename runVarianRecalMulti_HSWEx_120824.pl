#!/usr/bin/env perl
#run VariantRecalibrator for AB, SSEs, Alignment, and mapping biases
use POSIX;

#read in user parameters
my $dataset  = shift @ARGV; #reference fasta file

my $res = 0;



 if ($res == 0) {
 my $samCmd = "perl ~/bioinfo/runVariantRecalHetRef.pl HSWEx";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "perl ~/bioinfo/runVariantRecalHomRef.pl HSWEx";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "perl ~/bioinfo/runVariantRecalHomVar.pl HSWEx";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

