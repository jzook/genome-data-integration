#Adjustments to options needed for each run:
#integration applet
#input vcf, -ivcf
#callable loci for chromosome, -icallablelocibed
#exome assay targets, -itargetsbed
#split by chromosome, -ichrom
#output directory, -destination

#Run for HG004 on 7/8/16, took 1 min per chromosome
#No MT or Y run because callable loci for chrMT and Y failed
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_1_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=1 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_2_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=2 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_3_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=3 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_4_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=4 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_5_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=5 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_6_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=6 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_7_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=7 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_8_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=8 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_9_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=9 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_10_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=10 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_11_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=11 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_12_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=12 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_13_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=13 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_14_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=14 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_15_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=15 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_16_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=16 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_17_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=17 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_18_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=18 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_19_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=19 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_20_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=20 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_21_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=21 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_22_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=22 --destination=/HG004/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG004/IonExome/AmpliseqExome.20141120.NA24143.vcf -icallablelocibed=/HG004/IonExome/callableLoci_output/HG004_X_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG004/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=X --destination=/HG004/IonExome/Integration_prepare_ion_output/














