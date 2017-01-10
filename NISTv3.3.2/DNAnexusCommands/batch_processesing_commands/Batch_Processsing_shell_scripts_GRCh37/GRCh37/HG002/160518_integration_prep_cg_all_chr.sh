#Generate vcf and bed from Complete Genomics vcfBeta
#Adjustments to options needed for each run:
#input vcf, -ivcf
#split by chromosome, -ichrom

#Run for HG002 on 5/18/16, took approximately 15 min for a chromsome
#chrMT failed integration-prep
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=1 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=2 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=3 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=4 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=5 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=6 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=7 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=8 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=9 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=10 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=11 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=12 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=13 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=14 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=15 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=16 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=17 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=18 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=19 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=20 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=21 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=22 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=MT --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=X --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=Y --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/

#Justin updated applet to output bed file with "callable" loci for use in comparison
#Ran on 5/23/16 to get "callable" loci bed
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=1 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=2 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=3 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=4 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=5 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=6 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=7 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=8 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=9 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=10 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=11 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=12 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=13 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=14 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=15 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=16 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=17 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=18 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=19 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=20 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=21 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=22 --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/160523_output_for_bed/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=MT --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=X --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG002/Complete_Genomics/vcfBeta-GS000037263-ASM.vcf.bz2 -ichrom=Y --destination=/HG002/Complete_Genomics/Integration_prepare_cg_output/
