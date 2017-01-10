#Generate vcf and bed from Complete Genomics vcfBeta
#vcfbeta file for NA12878 downloaded from ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/CompleteGenomics_RMDNA_11272014/ASM/vcfBeta-GS000025639-ASM.vcf.bz2
#Adjustments to options needed for each run:
#input vcf, -ivcf
#split by chromosome, -ichrom

#Run 3/29/16, took approximately 13 min for a chromsome
#Run again on 4/14/16 -- Justin made adjustments to applet -- took 14 min per chromosome
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=1 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=2 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=3 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=4 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=5 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=6 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=7 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=8 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=9 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=10 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=11 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=12 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=13 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=14 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=15 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=16 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=17 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=18 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=19 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
#dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=20 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=21 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=22 --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=MT --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=X --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/
dx run -y GIAB:/Workflow/integration-prepare-cg -ivcf_in=/NA12878/Complete_Genomics/vcfBeta-GS000025639-ASM.vcf.bz2 -ichrom=Y --destination=/NA12878/Complete_Genomics/Integration_prepare_cg_output/

#Ran again on 7/29/16 for integration v3.3, Justin has ammended tool to output callable and not callable bed.  Justin already ran chr20 on 7/8 so did not re-run
#ran for 14 min, MT and Y failed
