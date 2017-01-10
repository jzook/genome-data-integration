#Combine individual chromosome 160628 integration v3.2.2 output files for HG003 for use in comparison and manual curation
#All run 7/8/16 bed combine ran for 1 min and vcf combine ran for 6-8 min (2cores)

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG003_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_callablemultinter_gt0 --destination=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG003_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_filteredsites --destination=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/*_highconf.bed"`; do bed_inputs="-ibeds=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG003_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_highconf --destination=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/*v3.2.2_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG003_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_all  --destination=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/  --instance-type=mem2_hdd2_x2 

#*_highconf_vcf.gz -- High-confidence integrated variants.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG003_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_highconf --destination=/HG003/Integration_v3.2.2_output/160708_Ill150-250-MP_CG_Ion_v3.2.2/ --instance-type=mem2_hdd2_x2