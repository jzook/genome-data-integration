#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=40:00:00
#$ -l h_vmem=4G
#$ -N GATK_PlatGen_UG
#$ -pe thread 1
#$ -t 2-21

BASE=/home/justin.zook
GATK=$BASE/GATK/GenomeAnalysisTK-2.5-2-gf57256b

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

java -jar -Xmx3g  $GATK/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T UnifiedGenotyper -glm BOTH \
	-I /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".bam \
	-L $LINE -U LENIENT_VCF_PROCESSING  -stand_call_conf 2 -stand_emit_conf 2 -mbq 10 -o /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_UG2.5_"$LINE1".vcf
