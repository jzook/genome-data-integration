my $res=0;
#bpsh 7 ~/bioinfo/bfast-0.6.5a/bfast/bfast fasta2brg -f ../references/100616_ERCC_Reference_noA_good.fasta -A 0 -t

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl 20120117_ceu_trio_b37_decoy/NA12878-NA12878.HiSeq.WGS.b37_decoy_Broad.recalibrated_splitMNP-nomnp_snps 1 1 1 2";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl 20120117_ceu_trio_b37_decoy/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf2.mbq10.snp.recal99 1 1 1 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl 20120117_ceu_trio_b37_decoy/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf202.snps.recal99 1 1 1 7";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
$res=0;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl 20120117_ceu_trio_b37_decoy/NA12878-CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.samtools.D200-nomnp_snps 1 1 1 5";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl CG/NA12878-CG_NA12878_all_nochr_sort_splitMNP-nomnp_snps 4 9 4 3";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl CG/NA12878-CG_NA12878_VQHIGH_nochr_sort_splitMNP-nomnp_snps 4 9 4 4";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl alignment/NA12878.all.SOLID.bfast.CEU.gatk.conf2010.snps.recal999 2 6 2 2";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl alignment/NA12878.all.LS454.ssaha2.CEU.samtools.filtD100_SNPs 3 7 3 5";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl Xprize/Illumina_SNPs_b37 1 5 1 6";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfHomHet.pl CEPH_NA12878/Variations/snps/ILLCLIA_allsnps 1 10 5 6";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
