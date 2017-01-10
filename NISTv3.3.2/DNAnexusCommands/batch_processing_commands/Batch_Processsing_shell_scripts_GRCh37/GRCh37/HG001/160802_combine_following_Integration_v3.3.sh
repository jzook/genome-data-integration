#Combine all chromosome integration v.3.3 outputs

#HG001 combined Chr1-X
#Run 8/2/16 ran 1 min for beds and 10 min for vcfs

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_CHROM1-X_v3.3_callablemultinter_gt0 --destination=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/*_filteredsites.bed"`; do bed_inputs="-ibeds=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_CHROM1-X_v3.3_filteredsites --destination=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/*_highconf.bed"`; do bed_inputs="-ibeds=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_CHROM1-X_v3.3_highconf --destination=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#need to increase memory to 2 cores
vcf_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/*v3.3_all.vcf.gz"`; do vcf_inputs="-ivcfs=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_CHROM1-X_v3.3_all  --destination=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/ --instance-type=mem2_hdd2_x2

#*_highconf_vcf.gz -- High-confidence integrated variants.
vcf_inputs=""
for l in `dx ls "GIAB:/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid_CHROM1-X_v3.3_highconf --destination=/NA12878/Integration_v3.3_output/160801_CG-IllFB-IllGATKHC-Ion-Solid-10X_v3.3/