#Generate vcf and bed from Complete Genomics vcfBeta
#Adjustments to options needed for each run:
#input vcf, -ivcf
#split by chromosome, -ichrom

#Run for HG005 on 7/22/16, took approximately 14 min for a chromsome
# MT Failed 
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=1 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=2 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=3 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=4 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=5 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=6 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=7 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=8 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=9 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=10 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=11 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=12 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=13 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=14 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=15 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=16 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=17 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=18 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=19 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=20 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=21 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=22 --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=MT --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=X --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG005/Complete_Genomics/vcfBeta-GS000037265-ASM.vcf.bz2 -ichrom=Y --destination=/HG005/Complete_Genomics/Integration_prepare_cg_output/

