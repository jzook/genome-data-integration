#Combine all chromosome integration v.3.3.1 outputs

# NEW COMBINE APPLET TO DEAL WITH CHR PREFIX ONLY FOR GRCH38 vcf-combineallchromwithchr'

#HG001 combined Chr1-X
#Run 10/5-10/6 ran 1 min for beds and 10 min for vcfs

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_callablemultinter_gt0 --destination=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_filteredsites --destination=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/*_highconf.bed"`; do bed_inputs="-ibeds=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf --destination=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/


START HERE
#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#need to increase memory to 2 cores
vcf_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/*v3.3_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchromwithchr $vcf_inputs -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_all  --destination=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/ --instance-type=mem2_hdd2_x2

#*_highconf_vcf.gz -- High-confidence integrated variants.
vcf_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchromwithchr $vcf_inputs -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf --destination=/HG001/GRCh38/Integration_v3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/

#compare to PG 
#Run 10/6/16  ran 18 min
dx run -y GIAB:/Workflow/compare-callsets-anyref -ivcfhighconfgz=/HG001/GRCh38/Integration_v.3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.vcf.gz -ivcfhighconftbi=/HG001/GRCh38/Integration_v.3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.vcf.gz.tbi -ibedhighconf=/HG001/GRCh38/Integration_v.3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.bed
-ivcftestgz=/HG001/GRCh38/PlatinumGenomes/PG2016-1.0_NA12878_b38.vcf.gz
-ivcftesttbi=/HG001/GRCh38/PlatinumGenomes/PG2016-1.0_NA12878_b38.vcf.gz.tbi
-ibedtest=/HG001/GRCh38/PlatinumGenomes/PG2016-1.0_NA12878_b38_ConfidentRegions_slop50.bed -isvbed=/HG001/GRCh38/remapped_PacBio_MetaSV_svclassify_mergedSVs.bed -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_1_v.3.3.1_b38_compare_PGslop50 -irtgsdf=/assets/rtgsdf38.tar.gz --destination=/HG001/GRCh38/Integration_v.3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/comp3.3.1_b38_PGslop50/