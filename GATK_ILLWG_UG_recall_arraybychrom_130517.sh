#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=20:00:00
#$ -l h_vmem=5G
#$ -N GATK_ILLWG_UGrecall
#$ -pe thread 1
#$ -t 15

BASE=/home/justin.zook

TEMP_LIST=$BASE/GATK/runs/ChromSizes.txt
TEMP_LIST1=$BASE/GATK/runs/ChromNames.txt



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

java -jar -Xmx4g  ~/bcbio.variation-0.0.5-standalone.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T UnifiedGenotyper -glm BOTH -I /scratch/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU.sorted.bam \
	$LINE -alleles:consensus /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_"$LINE1".vcf -U LENIENT_VCF_PROCESSING -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -mbq 10 \
	-A AlleleBalance -A BaseQualityRankSumTest -A DepthOfCoverage -A FisherStrand -A HaplotypeScore -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZeroFraction -A MappingQualityZero -A QualByDepth -A ReadPosRankSumTest -A RMSMappingQuality -A ClippingRankSumTest -A ReadPosEndDist -A MeanNeighboringBaseQuality \
	-o /scratch/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_UGrecall130517_"$LINE1".vcf -nt 1


