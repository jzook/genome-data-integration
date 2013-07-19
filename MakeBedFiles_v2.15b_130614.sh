#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=4:00:00
#$ -l h_vmem=4G
#$ -N MakeBedFiles_v2.15b
#$ -pe thread 5
#$ -t 1



BASE=/home/justin.zook
GATK=$BASE/GATK/GenomeAnalysisTK-2.5-2-gf57256b

TEMP_LIST=$BASE/scripts/NA12878_datasets_abbr.txt



NUM_TASKS=`cat $TEMP_LIST | wc -l`
if [ $SGE_TASK_ID -gt $NUM_TASKS ]; then
    echo "Error: No array job for task id $SGE_TASK_ID"
    exit 1
fi


# Debug
echo -e "NUM_TASKS: $NUM_TASKS\n"


LINE=`head -n $SGE_TASK_ID $TEMP_LIST | tail -n 1`


 
echo "Running job $JOB_NAME, $JOB_ID on $HOSTNAME"

#    perl ~/scripts/CombineChromVcf_FDA.pl /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_VQSRv2.15 & 
#    perl ~/scripts/CombineChromVcf_FDA.pl /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15 & 
#    perl ~/scripts/CombineChromVcf_FDA.pl /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.15 &
#    perl ~/scripts/CombineChromVcf_FDA.pl /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_Het_VQSRv2.15 
#    perl ~/scripts/CombineChromVcf_FDA.pl /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomVar_VQSRv2.15 &
#    perl ~/scripts/CombineChromVcf_FDA.pl /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomRef_VQSRv2.15 
    
#    java -jar ~/GATK/bcbio.variation-0.0.9-SNAPSHOT-standalone.jar variant-prep --keep-ref /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.15_all.vcf /scratch/justin.zook/references/human_g1k_v37.fasta
#    java -jar ~/GATK/bcbio.variation-0.0.9-SNAPSHOT-standalone.jar variant-prep --keep-ref /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf /scratch/justin.zook/references/human_g1k_v37.fasta
#    java -jar ~/GATK/bcbio.variation-0.0.9-SNAPSHOT-standalone.jar variant-prep --keep-ref /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_Het_VQSRv2.15_all.vcf /scratch/justin.zook/references/human_g1k_v37.fasta
#    java -jar ~/GATK/bcbio.variation-0.0.9-SNAPSHOT-standalone.jar variant-prep --keep-ref /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomVar_VQSRv2.15_all.vcf /scratch/justin.zook/references/human_g1k_v37.fasta
#    java -jar ~/GATK/bcbio.variation-0.0.9-SNAPSHOT-standalone.jar variant-prep --keep-ref /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomRef_VQSRv2.15_all.vcf /scratch/justin.zook/references/human_g1k_v37.fasta
#    java -jar ~/GATK/bcbio.variation-0.0.9-SNAPSHOT-standalone.jar variant-prep --keep-ref /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_VQSRv2.15_all.vcf /scratch/justin.zook/references/human_g1k_v37.fasta
    

#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomRef_VQSRv2.15_all.vcf \
-U LENIENT_VCF_PROCESSING -selectType SNP \
 -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomRef_VQSRv2.15_all_snp.vcf

#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.15_all.vcf \
-U LENIENT_VCF_PROCESSING -selectType SNP \
 -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.15_all_snp.vcf


#grep '^#\|PASS' /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomRef_VQSRv2.15_all_snp.vcf > /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomRef_VQSRv2.15_all_snp_PASS.vcf 

#add new dataset's callableloci if necessary
#perl $BASE/scripts/CombineChromVcf_nohead_FDA.pl /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign .frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly
#cat /scratch/justin.zook/NA12878/bed/union11callableonlymerged_all_sort.bed /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_all.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonly.bed
#sort -k 1,1 -k2,2 -k3,3 -n /scratch/justin.zook/NA12878/bed/union12callableMQonly.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonly_sort.bed
# $BASE/bedtools-2.17.0/bin/mergeBed -i /scratch/justin.zook/NA12878/bed/union12callableMQonly_sort.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged.bed
#grep -nm 1 '^1' /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged.bed
#6040:1	1	10000
#15981 /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged.bed
#1246396:1	1	10000
#wc -l /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged.bed
# 1265898 bed/union12callableMQonlymerged.bed
# exit
#head -6039 /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_XYMT.bed
#tail -9942 /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_1-22.bed
#sort -k 1,1 -k2,2n -k3,3n /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_XYMT.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonly_XYMT_sort.bed
# $BASE/bedtools-2.17.0/bin/mergeBed -i /scratch/justin.zook/NA12878/bed/union12callableMQonly_XYMT_sort.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_XYMT_sort.bed
#cat /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_1-22.bed /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_XYMT_sort.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_all_sort.bed


java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf -U LENIENT_VCF_PROCESSING -L /scratch/justin.zook/NA12878/bed/unionall12.cov3.bed  | grep -v '^#' | wc -l > /scratch/justin.zook/NA12878/bed/unionall12.cov3_v2.15b_countvar.txt & 

#add in sites called by consensus method
 $BASE/bedtools-2.17.0/bin/complementBed -i /scratch/justin.zook/NA12878/bed/unionall12.cov3.bed -g $BASE/bedtools-2.17.0/genomes/human.b37.genome > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_comp_v2.15b.bed
$BASE/bedtools-2.17.0/bin/subtractBed -a /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_comp_v2.15b.bed -b /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HomRef_VQSRv2.15_all_snp_PASS.vcf | $BASE/bedtools-2.17.0/bin/subtractBed -a stdin -b /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.15_all_snp.vcf > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_comp_addcert_v2.15b.bed
$BASE/bedtools-2.17.0/bin/complementBed -i /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_comp_addcert_v2.15b.bed -g $BASE/bedtools-2.17.0/genomes/human.b37.genome | grep -v '^Y' > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_v2.15b.bed

~/bedtools-2.17.0/bin/flankBed -i /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all_filtered.vcf -b 10 -g ~/bedtools-2.17.0/genomes/human.b37.genome | less

java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf \
-U LENIENT_VCF_PROCESSING -L /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_v2.15b.bed \
 | grep -v '^#' | wc -l > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_v2.15b_countvar.txt &

#grep -v 'PASS' /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf > /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all_filtered.vcf 

$BASE/bedtools-2.17.0/bin/subtractBed -a /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_v2.15b.bed -b /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all_filtered.vcf > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_v2.15b.bed

java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf \
-U LENIENT_VCF_PROCESSING -L /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_v2.15b.bed \
 | grep -v '^#' | wc -l > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_comp_addcert_nouncert_v2.15b_countvar.txt &

#subtract simple repeat regions that aren't covered completely by >4 reads in at least one dataset
$BASE/bedtools-2.17.0/bin/subtractBed -a /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_v2.15b.bed -b /scratch/justin.zook/NA12878/bed/All.simplerepeatsnocov.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_v2.15b.bed

java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf \
-U LENIENT_VCF_PROCESSING -L /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_v2.15b.bed \
 | grep -v '^#' | wc -l > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_v2.15b_countvar.txt &

#subtract known seg dups
$BASE/bedtools-2.17.0/bin/subtractBed -a /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_v2.15b.bed -b /scratch/justin.zook/NA12878/bed/superdupsmerged_all_sort.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_v2.15b.bed

java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf \
-U LENIENT_VCF_PROCESSING -L /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_v2.15b.bed \
 | grep -v '^#' | wc -l > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_v2.15b_countvar.txt &

#subtract regions manually determined to be in decoy
$BASE/bedtools-2.17.0/bin/subtractBed -a /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_v2.15b.bed -b /scratch/justin.zook/NA12878/bed/decoy.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_v2.15b.bed

java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf \
-U LENIENT_VCF_PROCESSING -L /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_v2.15b.bed \
 | grep -v '^#' | wc -l > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_v2.15b_countvar.txt &

#subtract regions in CNVs
$BASE/bedtools-2.17.0/bin/subtractBed -a /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_v2.15b.bed -b /scratch/justin.zook/NA12878/bed/dbvar_NA12878_CNVs_merged_all_sort.bed > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_noCNVs_v2.15b.bed

java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_allcall_UGHapMerge_HetHomVarAll_VQSRv2.15_all.vcf \
-U LENIENT_VCF_PROCESSING -L /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_noCNVs_v2.15b.bed \
 | grep -v '^#' | wc -l > /scratch/justin.zook/NA12878/bed/union12callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_noCNVs_v2.15b_countvar.txt &




