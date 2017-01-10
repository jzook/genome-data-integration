#Combine all chromosome integration v.3.3.2 outputs

# NEW COMBINE APPLET TO DEAL WITH CHR PREFIX ONLY FOR GRCH38 vcf-combineallchromwithchr'

#HG001 combined Chr1-X
#Run 11/8/16 ran 1 min for beds and 10 min for vcfs

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_callablemultinter_gt0 --destination=/HG001/GRCh38/Integration_v3..3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_filteredsites --destination=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/*_highconf.bed"`; do bed_inputs="-ibeds=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf --destination=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/


START HERE
#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#need to increase memory to 2 cores
vcf_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/*v.3.3.2_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchromwithchr $vcf_inputs -iprefix=HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_all  --destination=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/ --instance-type=mem2_hdd2_x2

#*_highconf_vcf.gz -- High-confidence integrated variants.
vcf_inputs=""
for l in `dx ls "GIAB:/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchromwithchr $vcf_inputs -iprefix=HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf --destination=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/

#*_annotated_vcf.gz
-iprefix=HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_annotated
--instance-type=mem2_hdd2_x2

#compare to PG
#Run 11/14/16  ran 30 min
dx run -y GIAB:/Workflow/compare-callsets-anyref -ivcfhighconfgz=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz -ivcfhighconftbi=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz.tbi -ibedhighconf=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.bed -ivcftestgz=/HG001/GRCh38/PlatinumGenomes/PG2016-1.0_NA12878_b38.vcf.gz -ivcftesttbi=/HG001/GRCh38/PlatinumGenomes/PG2016-1.0_NA12878_b38.vcf.gz.tbi -ibedtest=/HG001/GRCh38/PlatinumGenomes/PG2016-1.0_NA12878_b38_ConfidentRegions_slop50.bed -isvbed=/HG001/GRCh38/remapped_PacBio_MetaSV_svclassify_mergedSVs.bed -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2_b38_compare_PGslop50_b38 -irtgsdf=/assets/rtgsdf38.tar.gz --destination=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/comp3.3.2_b38_to_PG/


#compare to v3.3.1 b38
#Run 11/14/16  ran 30 min
dx run -y GIAB:/Workflow/compare-callsets-anyref -ivcfhighconfgz=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz -ivcfhighconftbi=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz.tbi -ibedhighconf=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.bed -ivcftestgz=/HG001/GRCh38/Integration_v.3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.vcf.gz -ivcftesttbi=/HG001/GRCh38/Integration_v.3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.vcf.gz.tbi -ibedtest=/HG001/GRCh38/Integration_v.3.3.1_output/161004_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.1/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.bed -isvbed=/HG001/GRCh38/remapped_PacBio_MetaSV_svclassify_mergedSVs.bed -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2_b38_compare_v3.3.1_b38 -irtgsdf=/assets/rtgsdf38.tar.gz --destination=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/comp3.3.2_b38_to_3.3.1_b38/


#compare b38 v3.3.2 HC to b37-->b38 lifted over v3.3.2
#Run 11/14/16 
dx run -y GIAB:/Workflow/compare-callsets-anyref -ivcfhighconfgz=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz -ivcfhighconftbi=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz.tbi -ibedhighconf=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.bed -ivcftestgz=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/comp3.3.2_b38_TO_b37tob38_liftover_highconfidenceV3.3.2/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_LIFTOVER_to_h38.vcf.gz -ivcftesttbi=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/comp3.3.2_b38_TO_b37tob38_liftover_highconfidenceV3.3.2/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_LIFTOVER_to_h38.vcf.gz.tbi -ibedtest=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/comp3.3.2_b38_TO_b37tob38_liftover_highconfidenceV3.3.2/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_LIFTOVER_to_h38.bed -isvbed=/HG001/GRCh38/remapped_PacBio_MetaSV_svclassify_mergedSVs.bed -iprefix=HG001_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_1_v.3.3.1_b38_compare_b37tob38_liftover_highconfidence_calls -irtgsdf=/assets/rtgsdf38.tar.gz --destination=/HG001/GRCh38/Integration_v.3.3.2_output/161107_CG-IllFB-IllGATKHC-Ion-10X-SOLID_v.3.3.2/comp3.3.2_b38_TO_b37tob38_liftover_highconfidenceV3.3.2/
