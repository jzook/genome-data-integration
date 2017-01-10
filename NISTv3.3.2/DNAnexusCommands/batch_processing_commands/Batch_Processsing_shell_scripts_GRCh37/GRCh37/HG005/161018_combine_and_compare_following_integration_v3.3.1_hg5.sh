#Combine all chromosome integration v.3.3.1 outputs

# NEW COMBINE APPLET TO DEAL WITH CHR PREFIX ONLY FOR GRCh37 vcf-combineallchrom'

#HG005 combined Chr1-22 -- female but just doing 1-22 for phasing
#Run 10/18 ran 1 min for beds and 10 min for vcfs

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.1_callablemultinter_gt0 --destination=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.1_filteredsites --destination=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/*_highconf.bed"`; do bed_inputs="-ibeds=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.1_highconf --destination=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/

#### MUST USE vcf-combineallchrom with GRCh37,  withchr is for GRCh38 only

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#need to increase memory to 2 cores
vcf_inputs=""
for l in `dx ls "GIAB:/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/*v3.3.1_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.1_all  --destination=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/ --instance-type=mem2_hdd2_x2

#*_highconf_vcf.gz -- High-confidence integrated variants.
#need to increase memory to 2 cores
vcf_inputs=""
for l in `dx ls "GIAB:/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3 --destination=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/ --instance-type=mem2_hdd2_x2

#no compare to PG, there is only PG for HG001

#compare to v3.3
#Run 10/18/16  ran 28 min
dx run -y GIAB:/Workflow/compare-callsets-anyref -ivcfhighconfgz=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.1_highconf.vcf.gz -ivcfhighconftbi=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.1_highconf.vcf.gz.tbi -ibedhighconf=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.1_highconf.bed -ivcftestgz=/HG005/GRCh37/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v3.3_highconf.vcf.gz -ivcftesttbi=/HG005/GRCh37/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v3.3_highconf.vcf.gz.tbi -ibedtest=/HG005/GRCh37/Integration_v3.3_output/160908_CG-IllFB-IllGATKHC-Ion_Solid_v3.3/HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v3.3_highconf.bed -isvbed=/HG005/GRCh37/HG005_FB_GATKHC_CG_allsvs_merged.bed -iprefix=HG005_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1_b37_compare_v3.3_b37 -irtgsdf=/assets/rtgsdf37.tar.gz --destination=/HG005/GRCh37/Integration_v.3.3.1_output/161018_CG-IllFB-IllGATKHC-Ion-SOLID_v.3.3.1/comp3.3.1_b37_to_v3.3_b37/

# !!!!!!individual integration files erroneously included 10X in file name, 10X callsets were not used in 10/18 integration.  I tried to 
# !!!!!! naming on DNA nexus but could not get the script below to work.  The combined files do have the correct naming.
for file in `dx ls GIAB:/HG005/GRCh37/Integration_v.3.3.1_output/test/*`; do dx mv -- "$file" "${file/HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID HG005_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-SOLID}"; done
HG002_NA24385_SRR1767407_IonXpress_020_rawlib_24030.1.bam