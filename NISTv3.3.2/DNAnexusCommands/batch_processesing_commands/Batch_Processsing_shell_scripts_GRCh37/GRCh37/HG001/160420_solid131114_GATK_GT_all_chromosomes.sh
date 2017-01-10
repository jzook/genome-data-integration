#To run GATK-tenotype-gvcfs to convert gvcf output to vcf
#Adjustments to options needed for each run:
#input "vcf", -ivcfs
#input .vcf.tbi, -ivcfs
#output vcf filename prefix, -iprefix

#Run 4/21/16 on solid 130920 run, took ~6-8 min per chromosome
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_1_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_1_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_1_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_2_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_2_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_2_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_3_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_3_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_3_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_4_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_4_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_4_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_5_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_5_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_5_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_6_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_6_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_6_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_7_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_7_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_7_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_8_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_8_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_8_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_9_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_9_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_9_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_10_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_10_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_10_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_11_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_11_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_11_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_12_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_12_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_12_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_13_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_13_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_13_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_14_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_14_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_14_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_15_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_15_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_15_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_16_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_16_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_16_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_17_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_17_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_17_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_18_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_18_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_18_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_19_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_19_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_19_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_20_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_20_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_20_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_21_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_21_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_21_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_22_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_22_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_22_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_MT_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_MT_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_MT_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_X_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_X_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_X_hg19_solid5500_SE75bp_GATKHC
dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5 -ivcfs=HG001_Y_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz -ivcfs=HG001_Y_hg19_solid5500_SE75bp_GATKHC_gvcf.vcf.gz.tbi -iprefix=HG001_Y_hg19_solid5500_SE75bp_GATKHC
