#!/usr/bin/env bash
## Variant Calling and integration prep for Chinese parents 10X data

## Paths to reference genomes
GENOME37=/assets/hg19_10X_210.fa.gz
GENOME37IDX=/assets/hg19_10X_210.fa.fai
REF37=/assets/hg19_10X_2.1.0.fasta-index.tar.gz
GENOME38=/assets/GRCh38_10X_2.1.0.fa.gz
GENOME38IDX=/assets/GRCh38_10X_2.1.0.fa.fai
REF38=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz


## HG006 #######################################################################
# ######### GRCh37
# bash giab_10X_variant_calling_pipeline.sh \
#     --hg HG006 \
#     --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_06202017/GM24694_ChineseTrioFather_hg19/GM24694_ChineseTrioFather_hg19_phased_possorted_bam.bam \
#     --refid GRCh37 \
#     --genome ${GENOME37} \
#     --genomeidx ${GENOME37IDX} \
#     --ref ${REF37} \
#     --rmsg NA24694 \
#     --hasY
#
# ######### GRCh38
# bash giab_10X_variant_calling_pipeline.sh \
#     --hg HG006 \
#     --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_06202017/GM24631_ChineseTrioFather_GRCh38/GM24631_ChineseTrioFather_GRCh38_phased_possorted_bam.bam \
#     --refid GRCh38 \
#     --genome ${GENOME38} \
#     --genomeidx ${GENOME38IDX} \
#     --ref ${REF38} \
#     --rmsg NA24694 \
#     --hasY

## HG007 #######################################################################
######### GRCh37
# bash giab_10X_variant_calling_pipeline.sh \
#     --hg HG007 \
#     --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_06202017/GM24695_ChineseTrioMother_hg19/GM24695_ChineseTrioMother_hg19_phased_possorted_bam.bam \
#     --refid GRCh37 \
#     --genome ${GENOME37} \
#     --genomeidx ${GENOME37IDX} \
#     --ref ${REF37} \
#     --rmsg NA24695

######### GRCh38
bash giab_10X_variant_calling_pipeline.sh \
    --hg HG007 \
    --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_06202017/GM24631_ChineseTrioMother_GRCh38/GM24631_ChineseTrioMother_GRCh38_phased_possorted_bam.bam \
    --refid GRCh38 \
    --genome ${GENOME38} \
    --genomeidx ${GENOME38IDX} \
    --ref ${REF38} \
    --rmsg NA24695
