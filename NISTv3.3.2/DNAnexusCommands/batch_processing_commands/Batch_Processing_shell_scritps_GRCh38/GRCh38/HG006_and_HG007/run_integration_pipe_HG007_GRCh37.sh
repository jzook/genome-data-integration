#!/usr/bin/bash
## Running Integration Pipeline GRCh37 HG007

for i in {1..22}; do

    ## Defining instance size - first 15 chroms need more memory
    INSTANCE=mem3_hdd2_x2
    if [[ ${i} -gt 15 ]]; then
      INSTANCE=mem2_hdd2_x2
    fi

    dx run -y GIAB:/Workflow/nist-integration-v3.3 \
        -ivcfs=GIAB:/HG007/GRCh37/10XGenomics/Integration_prepare_10X_output_v3.3/HG007_${i}_GRCh37_10X_sentieonHCbyhaplo.vcf.gz \
        -ivcfs=GIAB:/HG007/GRCh37/Complete_Genomics/Integration_prepare_cg_output/vcfBeta-GS000037477-ASM_${i}.vcf.gz \
        -ivcfs=GIAB:/HG007/GRCh37/Illumina/Illumina_GRCh37_Hiseq100X/FreeBayes_output/HG007_${i}_GRCh37_novoalign_Hiseq100X_FB.vcf.gz \
        -ivcfs=GIAB:/HG007/GRCh37/Illumina/Illumina_GRCh37_Hiseq100X/Sentieon_output/HG007_${i}_GRCh37_novoalign_Hiseq100X_sentieonHC.vcf.gz \
        -ivcfs=GIAB:/HG007/GRCh37/Illumina/Illumina_GRCh37_6Kb_MatePair/FreeBayes_output/HG007_${i}_GRCh37_novoalign_6Kb_MatePair_FB.vcf.gz \
        -ivcfs=GIAB:/HG007/GRCh37/Illumina/Illumina_GRCh37_6Kb_MatePair/Sentieon_output/HG007_${i}_GRCh37_novoalign_6Kb_MatePair_sentieonHC.vcf.gz \
        -ibeds=GIAB:/HG007/GRCh37/10XGenomics/Integration_prepare_10X_output_v3.3/HG007_${i}_GRCh37_10X_sentieonHCbyhaplo_callable.bed \
        -ibeds=GIAB:/HG007/GRCh37/Complete_Genomics/Integration_prepare_cg_output/vcfBeta-GS000037477-ASM_callable_${i}.bed \
        -ibeds=GIAB:/HG007/GRCh37/Illumina/Illumina_GRCh37_Hiseq100X/Integration_prepare_sentieon_v.3.3.2/HG007_${i}_GRCh37_novoalign_Hiseq100X_sentieonHC_gvcf_callable.bed \
        -ibeds=GIAB:/HG007/GRCh37/Illumina/Illumina_GRCh37_Hiseq100X/CallableLoci_output/HG007_${i}_GRCh37_novoalign_Hiseq100X_callableloci.bed \
        -ibeds=GIAB:/HG007/GRCh37/Illumina/Illumina_GRCh37_6Kb_MatePair/Integration_prepare_sentieon_v.3.3.2/HG007_${i}_GRCh37_novoalign_6Kb_MatePair_sentieonHC_gvcf_callable.bed \
        -ibeds=GIAB:/HG007/GRCh37/Illumina/Illumina_GRCh37_6Kb_MatePair/CallableLoci_output/HG007_${i}_GRCh37_novoalign_6Kb_MatePair_callableloci.bed \
        -iannotations=/Annotation_files/GATK_Annotations_160509.txt \
        -iannotations=/Annotation_files/Freebayes_Annotations_160408.txt \
        -iannotations=/Annotation_files/CG_Annotations_160408.txt \
        -icallsettable=callset_tables/180730_HG007_callset_tables/HG007_Datasets_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_GRCh37_${i}.txt \
        -ifiltbeds=/HG005/GRCh37/HG005_HG006_HG007_FB_GATKHC_CG_MetaSV_allsvs_merged.bed \
        -ifiltbeds=/filtbeds/GRCh37/AllRepeats_lt51bp_gt95identity_merged_slop5.bed.gz \
        -ifiltbeds=/filtbeds/GRCh37/AllRepeats_51to200bp_gt95identity_merged_slop5.bed.gz \
        -ifiltbeds=/filtbeds/GRCh37/AllRepeats_gt200bp_gt95identity_merged_sort.bed \
        -ifiltbeds=/filtbeds/GRCh37/oldfiles/hg19_self_chain_split_withalts_gt10k.bed.gz \
        -ifiltbeds=/filtbeds/GRCh37/SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz \
        -ifiltbeds=/filtbeds/GRCh37/superdupsmerged_all_sort.bed \
        -ifiltbeds=/filtbeds/GRCh37/mm-2-merged.bed \
        -ichrom=${i} \
        -iprefix=HG007_GIAB_GRCh37_highconf_CG-IllFB-IllSNT-10X_${i}_v3.3 \
        --destination=GIAB:/HG007/GRCh37/Integration_v3.3_output/180723_CG-IllFB-IllSNT-10X_v3.3/ \
        --instance-type=${INSTANCE}
done
