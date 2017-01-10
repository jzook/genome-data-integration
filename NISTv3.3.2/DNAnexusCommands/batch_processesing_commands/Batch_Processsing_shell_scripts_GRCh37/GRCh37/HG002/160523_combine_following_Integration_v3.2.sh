#Combined individual chromosome integration output files for HG002 for use in comparison and manual curation

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2_callablemultinter_gt0 --destination=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2_filteredsites --destination=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/*_highconf.bed"`; do bed_inputs="-ibeds=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2_highconf --destination=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
vcf_inputs=""
for l in `dx ls "GIAB:/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/*v3.2_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2_all  --destination=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/

#*_highconf_vcf.gz -- High-confidence integrated variants.
vcf_inputs=""
for l in `dx ls "GIAB:/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2_highconf --destination=/HG002/Integration_v3.2_output/160520_Ill_CG_Ion_Solid_v3.2/