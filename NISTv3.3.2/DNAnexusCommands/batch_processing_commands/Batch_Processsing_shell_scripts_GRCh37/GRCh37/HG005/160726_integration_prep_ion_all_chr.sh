#Adjustments to options needed for each run:
#integration applet
#input vcf, -ivcf
#callable loci for chromosome, -icallablelocibed
#exome assay targets, -itargetsbed
#split by chromosome, -ichrom
#output directory, -destination

#Run for HG005 on 7/26/16, took 1 min per chromosome
#No MT run because failed callable loci
#Y failed callable loci
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_1_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=1 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_2_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=2 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_3_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=3 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_4_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=4 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_5_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=5 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_6_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=6 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_7_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=7 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_8_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=8 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_9_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=9 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_10_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=10 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_11_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=11 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_12_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=12 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_13_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=13 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_14_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=14 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_15_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=15 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_16_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=16 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_17_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=17 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_18_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=18 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_19_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=19 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_20_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=20 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_21_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=21 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_22_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=22 --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_X_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=X --destination=/HG005/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG005/IonExome/AmpliseqExome.20141120.NA24631.vcf -icallablelocibed=/HG005/IonExome/callableLoci_output/HG005_Y_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG005/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=Y --destination=/HG005/IonExome/Integration_prepare_ion_output/














