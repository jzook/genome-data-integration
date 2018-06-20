#!/usr/bin/bash
## Running the Complete genomics pipelines for Chinese Parents

## Common parameters
LIFTOVERCHAIN=../Liftover_scripts/hg19ToHg38.over.chain
QUERYREF=ucsc.hg19.fasta
TARGETREF=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
GENOMEWARP=/scratch/nolson/giab/genomewarp/target

## HG006 Father
bash giab_complete_genomics_pipeline.sh \
    --hg HG006 \
    --vcfurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/CompleteGenomics_HanTrio_ExtractedFromCoriellCells_SmallVariants_CGAtools_08082014/HG006_NA24694-huCA017E_father/ASM/vcfBeta-GS000037476-ASM.vcf.bz2 \
    --cgvcf vcfBeta-GS000037476-ASM \
    --liftover ${LIFTOVERCHAIN} \
    --queryref ${QUERYREF} \
    --tergetref ${TARGETREF} \
    --genomewarp ${GENOMEWAPR}

## HG007 Mother
bash giab_complete_genomics_pipeline.sh \
    --hg HG007 \
    --vcfurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/CompleteGenomics_HanTrio_ExtractedFromCoriellCells_SmallVariants_CGAtools_08082014/HG007_NA24695-hu38168_mother/ASM/vcfBeta-GS000037477-ASM.vcf.bz2 \
    --cgvcf vcfBeta-GS000037477-ASM \
    --liftover ${LIFTOVERCHAIN} \
    --queryref ${QUERYREF} \
    --tergetref ${TARGETREF} \
    --genomewarp ${GENOMEWAPR}
