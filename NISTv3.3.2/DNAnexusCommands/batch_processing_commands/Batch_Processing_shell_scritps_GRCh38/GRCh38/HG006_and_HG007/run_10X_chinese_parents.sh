#!/usr/bin/env bash
## Variant Calling and integration prep for Chinese parents 10X data

## Paths to reference genomes
GENOME37=NEEDTODEFINE
GENOME37IDX=NEEDTODEFINE
REF37=NEEDTODEFINE
GENOME38=/assets/GRCh38_10X_2.1.0.fa.gz
GENOME38IDX=/assets/GRCh38_10X_2.1.0.fa.fai
REF38=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz


## HG006 #######################################################################
######### GRCh37
bash giab_10X_variant_calling_pipeline.sh \
    --hg HG006 \
    --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh37/NA24694_LongRanger_phased_possorted.bam \
    --refid GRCh37 \
    --genome ${GENOME37} \
    --genomeidx ${GENOME37IDX}
    --ref ${REF37}
    --rmsg GM24694

######### GRCh38
bash giab_10X_variant_calling_pipeline.sh \
    --hg HG006 \
    --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh38/NA24694.GRCh38.phased_possorted_bam.bam \
    --refid GRCh38 \
    --genome ${GENOME37} \
    --genomeidx ${GENOME37IDX}
    --ref ${REF37}
    --rmsg GM24694


## HG007 #######################################################################
######### GRCh37
bash giab_10X_variant_calling_pipeline.sh \
    --hg HG007 \
    --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh37/NA24695_LongRanger_phased_possorted.bam \
    --refid GRCh37 \
    --genome ${GENOME37} \
    --genomeidx ${GENOME37IDX}
    --ref ${REF37}
    --rmsg GM24695

######### GRCh38
bash giab_10X_variant_calling_pipeline.sh \
    --hg HG007 \
    --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh38/NA24695.GRCh38.phased_possorted_bam.bam \
    --refid GRCh38 \
    --genome ${GENOME38} \
    --genomeidx ${GENOME38IDX}
    --ref ${REF38}
    --rmsg GM24695
