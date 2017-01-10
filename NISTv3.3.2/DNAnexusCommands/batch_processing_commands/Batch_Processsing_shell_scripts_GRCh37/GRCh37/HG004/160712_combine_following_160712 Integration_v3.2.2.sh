#Combine individual chromosome 160712/13 integration v3.2.2 output files for HG004 for use in comparison and manual curation
#Ran combine twice: CHR1-22 and CHR1-X (1-22+X)
#All run 7/12/16 bed combine ran for 1 min and vcf combine ran for 8-11 min (2cores)

#Combine Chromosomes 1-22 + X
#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-X_v3.2.2_callablemultinter_gt0 --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-X_v3.2.2_filteredsites --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*_highconf.bed"`; do bed_inputs="-ibeds=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-X_v3.2.2_highconf --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*v3.2.2_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-X_v3.2.2_all  --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/  --instance-type=mem2_hdd2_x2 

#*_highconf_vcf.gz -- High-confidence integrated variants.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-X_v3.2.2_highconf --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/ --instance-type=mem2_hdd2_x2



#Combine Chromosome 1-22
#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_callablemultinter_gt0 --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_filteredsites --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*_highconf.bed"`; do bed_inputs="-ibeds=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_highconf --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*v3.2.2_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_all  --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/  --instance-type=mem2_hdd2_x2 

#*_highconf_vcf.gz -- High-confidence integrated variants.
#Increased memory by adding an additional core
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG004_GIAB_highconf_IllFB-IllGATKHC-CG-Ion_CHROM1-22_v3.2.2_highconf --destination=/HG004/Integration_v3.2.2_output/160712_Ill150-250-MP_CG_Ion_v3.2.2/ --instance-type=mem2_hdd2_x2