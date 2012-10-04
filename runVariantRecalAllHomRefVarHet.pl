#!/usr/bin/env perl
#run VariantRecalibrator for AB, SSEs, Alignment, and mapping biases all together for HomVar, HomRef, and HetRef
use POSIX;

#read in user parameters
my $dataset  = shift @ARGV; #e.g., HSWG, CG, etc.

my $res = 1;



#$res=0;
 # HomVar: AB, QD, Alignment, and mapping annotations
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVarUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomVar.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an QD -an ABCI50 -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -an MQ -an MQ0Fraction -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 -nt 8 --percentBadVariants 0.03 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVarUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomVar.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an QD -an ABCI50 -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -an MQ -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 -nt 8 --percentBadVariants 0.03 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }

 }
#$res=0;
 # HomVar: AB, QD, Alignment, and mapping annotations
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T ApplyRecalibration -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVarUncert.vcf -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all.tranches_1.6 -o OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomVar_all_recal90.vcf -nt 8 --ts_filter_level 90.0";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

#$res=0;
 # HomRef: AB, QD, Alignment, and mapping annotations
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an QD -an ABCI50 -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -an MQ -an MQ0Fraction -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 -nt 8 --qualThreshold -1.0 --percentBadVariants 0.03 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an QD -an ABCI50 -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -an MQ -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 -nt 8 --qualThreshold -1.0 --percentBadVariants 0.03 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }
 }

 # HomRef: AB, QD, Alignment, and mapping annotations
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T ApplyRecalibration -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all.tranches_1.6 -o OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_all_recal90.vcf -nt 8 --ts_filter_level 90.0";
  print "$samCmd\n";
 $res = system($samCmd);
 }  

$res=0;
 # Het: SSEs, AB, QD, Alignment, and mapping annotations
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HetRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HetRef.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an FS -an BaseQRankSum -an QD -an ABCI50 -an HaplotypeScore -an ReadPosRankSum -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -an MQ -an MQRankSum -an MQ0Fraction -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 -nt 8 --percentBadVariants 0.03 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HetRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HetRef.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an FS -an BaseQRankSum -an QD -an ABCI50 -an HaplotypeScore -an ReadPosRankSum -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -an MQ -an MQRankSum -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 -nt 8 --percentBadVariants 0.03 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }
 }

 # Het: SSEs, AB, QD, Alignment, and mapping annotations
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T ApplyRecalibration -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HetRefUncert.vcf -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all.tranches_1.6 -o OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_Het_all_recal90.vcf -nt 8 --ts_filter_level 90.0";
  print "$samCmd\n";
 $res = system($samCmd);
 }  


 
