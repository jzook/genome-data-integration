#Combine individual chromosome 160831 integration v3.3 output files for HG003 for use in comparison and manual curation
#All run 8/31/16 bed combine ran for 1 min and vcf combine ran for 6-8 min (2cores)

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG003_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v3.3_callablemultinter_gt0 --destination=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG003_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v3.3_filteredsites --destination=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/*_highconf.bed"`; do bed_inputs="-ibeds=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG003_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v3.3_highconf --destination=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/*v3.3_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG003_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v3.3_all  --destination=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/  --instance-type=mem2_hdd2_x2 

#*_highconf_vcf.gz -- High-confidence integrated variants.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG003_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v3.3_highconf --destination=/HG003/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/ --instance-type=mem2_hdd2_x2