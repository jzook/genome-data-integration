my $res=0;

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_conf2_mbq10_raw_merged_SNPs_120212_HSWGannotate 1 1 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_conf2_mbq10_raw_merged_SNPs_120427_new_HSWGannotate 1 1 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_HSWGannotate 1 1 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_HSWExannotate 1 2 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_ILLWGannotate 1 3 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

#$res=0;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_ILLWExannotate 1 4 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

#$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_XIllannotate 1 5 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

#$res=0;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_SOLWGannotate 2 6 2";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_454WGannotate 3 7 3";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_454WExannotate 3 8 3";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

#$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_CGannotate 4 9 4";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_ILLCLIAannotate 1 10 5";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_conf2_mbq10_raw_merged_SNPs_120525_new_IllPCRFreeNoBQSRannotate 1 11 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
$res=1;
 if ($res == 0) {
 my $samCmd = "~/ActivePerl-5.12/bin/perl ~/bioinfo/sqlinsertVcfParams.pl OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_conf2_mbq10_raw_merged_SNPs_120427_IllPCRFreeNoBQSRannotate 1 11 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   


