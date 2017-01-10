#Adjustments to options needed for each run:
#integration applet
#input vcf, -ivcf
#callable loci for chromosome, -icallablelocibed
#exome assay targets, -itargetsbed
#split by chromosome, -ichrom
#output directory, -destination

#Run for HG003 on 6/3/16, took 1 min per chromosome
#Re-run on 8/25/16 because applet updated on 7/8 to output callable bed
#No MT run because callable loci for chrMT failed
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_1_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=1 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_2_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=2 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_3_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=3 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_4_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=4 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_5_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=5 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_6_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=6 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_7_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=7 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_8_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=8 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_9_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=9 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_10_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=10 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_11_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=11 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_12_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=12 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_13_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=13 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_14_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=14 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_15_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=15 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_16_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=16 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_17_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=17 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_18_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=18 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_19_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=19 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_20_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=20 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_21_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=21 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_22_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=22 --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_X_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=X --destination=/HG003/IonExome/Integration_prepare_ion_output/
dx run -y GIAB:/Workflow/integration-prepare-ion -ivcf_in=/HG003/IonExome/AmpliseqExome.20141120.NA24149.vcf -icallablelocibed=/HG003/IonExome/callableLoci_output/HG003_Y_hs37d5_IonExome_callableloci.bed -itargetsbed=/HG003/IonExome/AmpliseqExome.20141120_effective_regions.bed -ichrom=Y --destination=/HG003/IonExome/Integration_prepare_ion_output/














