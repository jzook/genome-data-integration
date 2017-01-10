#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_ALLCHROM_v3.2a_callablemultinter_gt0 --destination=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/*_filteredsites.bed"`; do bed_inputs="-ibeds=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_ALLCHROM_v3.2a_filteredsites --destination=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/*_highconf.bed"`; do bed_inputs="-ibeds=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_ALLCHROM_v3.2a_highconf --destination=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
vcf_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/*v3.2a_all.vcf.gz"`; do vcf_inputs="-ivcfs=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_ALLCHROM_v3.2a_all  --destination=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/

#*_highconf_vcf.gz -- High-confidence integrated variants.
vcf_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_ALLCHROM_v3.2a_highconf --destination=/NA12878/Integration_v3.2_output/160510_Ill_CG_Ion_Solid_v3.2a/