#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=15:00:00
#$ -l h_vmem=4G
#$ -N VcfClassify_2.15
#$ -pe thread 1
#$ -t 2-21

BASE=/home/justin.zook
GATK=$BASE/GATK/GenomeAnalysisTK-2.5-2-gf57256b

TEMP_LIST=$BASE/GATK/runs/ChromSizes.txt
TEMP_LIST1=$BASE/GATK/runs/ChromNames.txt
TEMP_DSET=$BASE/scripts/NA12878_datasets_abbr.txt



NUM_TASKS=`cat $TEMP_LIST | wc -l`
if [ $SGE_TASK_ID -gt $NUM_TASKS ]; then
    echo "Error: No array job for task id $SGE_TASK_ID"
    exit 1
fi


NUM_TASKS1=`cat $TEMP_LIST1 | wc -l`
if [ $SGE_TASK_ID -gt $NUM_TASKS1 ]; then
    echo "Error: No array job for task id $SGE_TASK_ID"
    exit 1
fi


# Debug
echo -e "NUM_TASKS: $NUM_TASKS\n"
echo -e "NUM_TASKS: $NUM_TASKS1\n"

# Chromosome range to call in
LINE=`head -n $SGE_TASK_ID $TEMP_LIST | tail -n 1` 
# Chromosome names
LINE1=`head -n $SGE_TASK_ID $TEMP_LIST1 | tail -n 1`


# export PATH=/home/mikem/GATK/GenomeAnalysisTK-2.1-11-g13c0244/resources:$PATH
 
echo "Running job $JOB_NAME, $JOB_ID on $HOSTNAME"

NUM_TASKS2=`cat $TEMP_DSET | wc -l`
for i in {1..100}
do
#if [ $i -gt $NUM_TASKS2 ]; then
if [ $i -gt 0 ]; then
    break
fi
DSET=`head -n $i $TEMP_DSET | tail -n 1`
#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -U LENIENT_VCF_PROCESSING -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Het_map_recal90_all.vcf $LINE -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Het_map_recal90_"$LINE1".vcf &
#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -U LENIENT_VCF_PROCESSING -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Het_SSE_recal90_all.vcf $LINE -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Het_SSE_recal90_"$LINE1".vcf &
#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -U LENIENT_VCF_PROCESSING -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Het_Align_recal90_all.vcf $LINE -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Het_Align_recal90_"$LINE1".vcf &
#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -U LENIENT_VCF_PROCESSING -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Het_AB_recal90_all.vcf $LINE -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Het_AB_recal90_"$LINE1".vcf
#sleep 30
#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -U LENIENT_VCF_PROCESSING -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Hom_map_recal90_all.vcf $LINE -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Hom_map_recal90_"$LINE1".vcf &
#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -U LENIENT_VCF_PROCESSING -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Hom_SSE_recal90_all.vcf $LINE -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Hom_SSE_recal90_"$LINE1".vcf &
#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -U LENIENT_VCF_PROCESSING -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Hom_Align_recal90_all.vcf $LINE -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Hom_Align_recal90_"$LINE1".vcf &
#java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -U LENIENT_VCF_PROCESSING -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Hom_AB_recal90_all.vcf $LINE -o /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_"$DSET"call_UGHapMerge_Hom_AB_recal90_"$LINE1".vcf
#sleep 30
ll $i rrr
done

#sleep 2m

perl ~/scripts/VcfClassifyUsingVQSR_snpindel_v2.15_FDA.pl /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_ $LINE1 HSWG ill ILLWG ill XIll ill 454WG 454 ILLCLIA ill IllPCRFree ill XPSolWGLS sol IonEx ion HSWEx ill ILLWEx ill CG cg PlatGen ill

