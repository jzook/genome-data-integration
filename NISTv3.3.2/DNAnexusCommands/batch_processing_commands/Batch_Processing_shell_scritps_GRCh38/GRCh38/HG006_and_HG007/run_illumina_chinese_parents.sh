#!/usr/bin/bash

### Running the giab_illumina_variant_calling_pipeline on the
### chinese parent Illumina datasets

## Father
# 100X HiSeq
bash giab_illumina_variant_calling_pipeline.sh \
  --hg HG006 \
  --dataset Hiseq100X \
  --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam \
  --baiurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam.bai

# MatePair
bash giab_illumina_variant_calling_pipeline.sh \
 --hg HG006 \
 --dataset 6Kb_MatePair \
 --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NIST_Stanford_Illumina_6kb_matepair/bams/GRCh38/HG006.sorted.bam \
 --baiurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NIST_Stanford_Illumina_6kb_matepair/bams/GRCh38/HG006.sorted.bam.bai


## Mother
# 100X HiSeq
bash giab_illumina_variant_calling_pipeline.sh \
  --hg HG007 \
  --dataset Hiseq100X \
  --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NA24695_Mother_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG007.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam \
  --baiurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NA24695_Mother_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG007.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam.bai

# MatePair
bash giab_illumina_variant_calling_pipeline.sh \
  --hg HG007 \
  --dataset 6Kb_MatePair \
  --bamurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NIST_Stanford_Illumina_6kb_matepair/bams/GRCh38/HG007.sorted.bam \
  --baiurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NIST_Stanford_Illumina_6kb_matepair/bams/GRCh38/HG007.sorted.bam.bai
