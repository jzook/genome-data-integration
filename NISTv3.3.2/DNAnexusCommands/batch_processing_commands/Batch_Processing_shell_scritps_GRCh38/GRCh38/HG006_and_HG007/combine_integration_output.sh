#!/usr/bin/bash
# Combine CHROM VCF and BED file integration pipeline output

## Input Arguments
PREFIX #integation pipeline output root name
DIR #integration pipeline output directory
CHR #for GRCh37 vs. GRCh38

## BED Files ###################################################################
#### BED File types
#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
#*_highconf.bed  --High-confidence BED file with at least 1 technology callable

for bed in callablemultinter_gt0 filteredsites highconf; do
    bed_inputs=""

    for i in {1:22}; do
        bed_inputs="-ibeds=${PREFIX}_${i}_${bed}.bed ${bed_inputs}"
    done

    dx run -y GIAB:/Workflow/bed-combineallchrom${CHR} \
      $bed_inputs \
      -iprefix=${PREFIX}_CHROM1-22_v.3.3.2_${bed} \
      --destination=${DIR}
done


## VCF Files ###################################################################
#### VCF File Types
#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
#*_highconf_vcf.gz -- High-confidence integrated variants.
#*_annotated.vcf.gz

for vcf in all highconf_ annotated; do
    vcf_inputs=""

    for i in {1:22}; do
        vcf_inputs="-ivcfs=${PREFIX}_${i}_${vcf}.vcf.gz ${vcf_inputs}"
    done

    dx run -y GIAB:/Workflow/vcf-combineallchrom${CHR} \
      $vcf_inputs \
      -iprefix=${PREFIX}_CHROM1-22_v.3.3.2_${vcf} \
      --destination=${DIR} \
      --instance-type=mem2_hdd2_x2
done
