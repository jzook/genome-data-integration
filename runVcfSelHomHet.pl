my $res=1;
#bpsh 7 ~/bioinfo/bfast-0.6.5a/bfast/bfast fasta2brg -f ../references/100616_ERCC_Reference_noA_good.fasta -A 0 -t

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl 20120117_ceu_trio_b37_decoy/NA12878-NA12878.HiSeq.WGS.b37_decoy_Broad.recalibrated_splitMNP-nomnp_snps &";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl 20120117_ceu_trio_b37_decoy/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf2.mbq10.snp.recal99";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl 20120117_ceu_trio_b37_decoy/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf202.snps.recal99";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl 20120117_ceu_trio_b37_decoy/NA12878-CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.samtools.D200-nomnp_snps &";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl CG/NA12878-CG_NA12878_all_nochr_sort_splitMNP-nomnp_snps &";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl CG/NA12878-CG_NA12878_VQHIGH_nochr_sort_splitMNP-nomnp_snps";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=1;

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl alignment/NA12878.all.SOLID.bfast.CEU.gatk.conf2010.snps.recal999";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl alignment/NA12878.all.LS454.ssaha2.CEU.samtools.filtD100_SNPs";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl Xprize/Illumina_SNPs_b37";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet_ILLCLIA.pl CEPH_NA12878/Variations/snps/ILLCLIA_allsnps";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/VcfSelHomHet.pl XPrize_Fosmid/solid_wgs/XPSolWGLS.gatk.conf202.mbq10.snp.recal99.vcf";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
