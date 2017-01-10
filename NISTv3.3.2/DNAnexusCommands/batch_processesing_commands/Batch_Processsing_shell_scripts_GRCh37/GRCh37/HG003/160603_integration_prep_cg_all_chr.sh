#Generate vcf and bed from Complete Genomics vcfBeta
#Adjustments to options needed for each run:
#input vcf, -ivcf
#split by chromosome, -ichrom

#Run for HG003 on 6/3/16, took approximately 15 min for a chromsome
#MT failed
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=1 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=2 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=3 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=4 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=5 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=6 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=7 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=8 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=9 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=10 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=11 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=12 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=13 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=14 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=15 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=16 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=17 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=18 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=19 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=20 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=21 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=22 --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=MT --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=X --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/HG003/Complete_Genomics/vcfBeta-GS000037264-ASM.vcf.bz2 -ichrom=Y --destination=/HG003/Complete_Genomics/Integration_prepare_cg_output/

