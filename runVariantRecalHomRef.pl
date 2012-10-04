#!/usr/bin/env perl
#run VariantRecalibrator for AB, SSEs, Alignment, and mapping biases for HomRef
use POSIX;

#read in user parameters
my $dataset  = shift @ARGV; #reference fasta file

my $res = 0;



 # SSEs
 if ($res == 2) { #TODO: add covariates for Hom SSEs
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an FS -an BaseQRankSum -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_SSE.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_SSE.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_SSE.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 # AB and QD
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an QD -an ABCI50 -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8 --percentBadVariants 0.03 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
  
 if ($res != 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an QD -an ABCI50 -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8 --percentBadVariants 0.08 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
 }

#$res=1;
 # Align
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) { #if it doesn't work, then increase percentBadVariants
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8  --percentBadVariants 0.08 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res != 0) { #if it doesn't work, then increase percentBadVariants
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8  --percentBadVariants 0.2 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res != 0) { #if it doesn't work, then increase percentBadVariants
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an ReadPosEndDist -an ReadMeanPos -an ReadMeanLen -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8  --percentBadVariants 0.08 --maxGaussians 2";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
}
#$res=1;
 # Mapping
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an MQ -an MQ0Fraction -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8 --percentBadVariants 0.08 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) { #if it doesn't work, then remove MQ0Fraction
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an MQ -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8 --percentBadVariants 0.08 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res != 0) { #if it doesn't work, then remove MQ0Fraction
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T VariantRecalibrator -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -resource:consensus1.6,known=false,training=true,truth=true,prior=20.0 OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_allcall_HomRef_no0qualor00GT.vcf -resource:dbSNP135,known=true,training=false,truth=false,prior=5.0 /data/results/justin/erccs/references/dbsnp_135.b37.vcf -an MQ -an DP -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.tranches_1.6 -rscriptFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 --qualThreshold -1.0 -nt 8 --percentBadVariants 0.15 --maxGaussians 1";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
}

#$res=0;
 # SSE
 if ($res == 2) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T ApplyRecalibration -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_SSE.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_SSE.tranches_1.6 -o OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_SSE_recal90.vcf -nt 8 --ts_filter_level 90.0";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 # AB and QD
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T ApplyRecalibration -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB.tranches_1.6 -o OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_AB_recal90.vcf -nt 8 --ts_filter_level 90.0";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 # Align
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T ApplyRecalibration -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align.tranches_1.6 -o OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_Align_recal90.vcf -nt 8 --ts_filter_level 90.0";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 # Mapping
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx22g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T ApplyRecalibration -input OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRefUncert.vcf -recalFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.recal_1.6 -tranchesFile OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map.tranches_1.6 -o OMNI_CG_1KG_SOL_ILL_454_HS_HSJaf_Hsam_ILLCLIA_XIll_IllPCRFree_EFosComp_IonEx_XFosIll_XFosSol_conf2_mbq10_raw_merged_SNPs_120620_${dataset}call_HomRef_map_recal90.vcf -nt 8 --ts_filter_level 90.0";
  print "$samCmd\n";
 $res = system($samCmd);
 }   




