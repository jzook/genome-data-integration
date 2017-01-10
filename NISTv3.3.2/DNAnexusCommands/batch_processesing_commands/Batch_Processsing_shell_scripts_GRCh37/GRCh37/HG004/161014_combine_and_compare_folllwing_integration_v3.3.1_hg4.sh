#Combine all chromosome integration v.3.3.1 outputs

# NEW COMBINE APPLET TO DEAL WITH CHR PREFIX ONLY FOR GRCh37 vcf-combineallchrom'

#HG004 combined Chr1-22 -- female but just doing 1-22 for phasing
#Run 10/14 ran 1 min for beds and 10 min for vcfs

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.1_callablemultinter_gt0 --destination=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.1_filteredsites --destination=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/*_highconf.bed"`; do bed_inputs="-ibeds=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG004_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.1_highconf --destination=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/

#### MUST USE vcf-combineallchrom with GRCh37,  withchr is for GRCh38 only

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#need to increase memory to 2 cores
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/*v3.3_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG004_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.1_all  --destination=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/ --instance-type=mem2_hdd2_x2

#*_highconf_vcf.gz -- High-confidence integrated variants.
#need to increase memory to 2 cores
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG004_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.1_highconf --destination=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/ --instance-type=mem2_hdd2_x2

#no compare to PG, there is only PG for HG001

#compare to v3.3
#Run 10/14/16  ran
dx run -y GIAB:/Workflow/compare-callsets-anyref -ivcfhighconfgz=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/HG004_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.1_highconf.vcf.gz -ivcfhighconftbi=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/HG004_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.1_highconf.vcf.gz.tbi -ibedhighconf=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/HG004_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.1_highconf.bed -ivcftestgz=/HG004/GRCh37/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/HG004_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v3.3_highconf.vcf.gz -ivcftesttbi=/HG004/GRCh37/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/HG004_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v3.3_highconf.vcf.gz.tbi -ibedtest=/HG004/GRCh37/Integration_v3.3_output/160831_CG-IllFB-IllGATKHC-Ion-10X_v3.3/HG004_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v3.3_highconf.bed -isvbed=/HG004/GRCh37/HG002_HG003_HG004_allsvs_merged.bed -iprefix=HG003_highconf_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1_b37_compare_v3.3_b37 -irtgsdf=/assets/rtgsdf37.tar.gz --destination=/HG004/GRCh37/Integration_v.3.3.1_output/161014_CG-IllFB-IllGATKHC-Ion-10X_v.3.3.1/comp3.3.1_b37_to_v3.3_b37/
