
#Adjustments to options needed for each run:
#integration applet
#input vcf, -ivcf
#callable loci for chromosome, -icallablelocibed
#exome assay targets, -itargetsbed
#split by chromosome, -ichrom
#output directory, -destination

#Run for HG002 on 5/18/16, took 1 min per chromosome
#Run again on 8/24/16, applet updated to output callable bed
#chrMT was not run because callable loci failed.  chrY failed integration-prep possible due to no variants.

dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_1_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=1 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_2_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=2 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_3_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=3 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_4_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=4 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_5_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=5 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_6_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=6 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_7_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=7 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_8_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=8 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_9_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=9 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_10_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=10 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_11_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=11 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_12_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=12 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_13_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=13 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_14_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=14 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_15_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=15 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_16_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=16 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_17_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=17 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_18_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=18 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_19_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=19 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_20_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=20 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_21_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=21 --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_22_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=22 --destination=/HG002/IonExome/Integration_prepare_ion_output/
#dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_MT_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=MT --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_X_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=X --destination=/HG002/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG002/IonExome/AmpliseqExome.20141120.NA24385.vcf -icallablelocibed=/HG002/IonExome/callableLoci_output/HG002_Y_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG002/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=Y --destination=/HG002/IonExome/Integration_prepare_ion_output/














