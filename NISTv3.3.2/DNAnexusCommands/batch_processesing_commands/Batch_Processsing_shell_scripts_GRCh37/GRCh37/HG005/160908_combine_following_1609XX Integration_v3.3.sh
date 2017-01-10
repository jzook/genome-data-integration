#Combine individual chromosome 160908 integration v3.3 output files for HG005 for use in comparison and manual curation
#combined for chromosomes 1-22
#All run 9/8/16 bed combine ran for 1 min and vcf combine ran for 8-11 min (2cores)

#Combine Chromosome 1-22
#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v3.3_callablemultinter_gt0 --destination=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v3.3_filteredsites --destination=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/*_highconf.bed"`; do bed_inputs="-ibeds=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v3.3_highconf --destination=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/*v3.3_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v3.3_all  --destination=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/  --instance-type=mem2_hdd2_x2 

#*_highconf_vcf.gz -- High-confidence integrated variants.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v3.3_highconf --destination=/HG005/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/ --instance-type=mem2_hdd2_x2