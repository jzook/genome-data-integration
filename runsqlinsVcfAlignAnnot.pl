my $res=1;
#bpsh 7 /Volumes/SSD960/workspace/data/vcfannot/bfast-0.6.5a/bfast/bfast fasta2brg -f ../references/100616_ERCC_Reference_noA_good.fasta -A 0 -t
print "test";
 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ HSWG 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ ILLWG 3";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ XIll 5";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ SOLWG 6";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ 454WG 7";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ CG 9";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ ILLCLIA 10";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ IllPCRFreeNoBQSR 11";
  print "$samCmd\n";
 $res = system($samCmd);
 } 

 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ XPSolWGLS 12";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

$res=0;
 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ HSWEx 2";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ ILLWEx 4";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res == 0) {
 my $samCmd = "perl /Volumes/SSD960/workspace/data/vcfannot/sqlinsVcfAlignAnnot.pl /Volumes/SSD960/workspace/data/vcfannot/OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_ 454WEx 8";
  print "$samCmd\n";
 $res = system($samCmd);
 }   


$res=1;
if ($res == 0) {
    my $samCmd = "sqlite3 /Volumes/SSD960/workspace/data/NA12878vcfCompare_120315/NA12878vcfs.db 'CREATE INDEX chrposannot_idx on vcfannot (chrompos)'";
    print "$samCmd\n";
    $res = system($samCmd);
}

if ($res == 0) {
    my $samCmd = "sqlite3 /Volumes/SSD960/workspace/data/NA12878vcfCompare_120315/NA12878vcfs.db 'CREATE INDEX bamannot_idx on vcfannot (bam)'";
    print "$samCmd\n";
    $res = system($samCmd);
}

