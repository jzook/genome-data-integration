#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=60:00:00
#$ -l h_vmem=6G
#$ -N GATK_NCIIonWG_UG
#$ -pe thread 1
#$ -t 2-21

BASE=/home/justin.zook
GATK=$BASE/GATK/GenomeAnalysisTK-2.7-2-g6bda569
export JAVA_HOME=/home/justin.zook/jre1.7.0_25
export PATH=$JAVA_HOME/bin:$PATH

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

# $BASE/samtools-0.1.18/samtools index /scratch/justin.zook/NA12878/MiSeq250/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130520.bam
java -jar -Xmx3g  $GATK/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T UnifiedGenotyper -glm BOTH \
	-I /scratch/justin.zook/NA12878/NCIIonWG/NCIIonWG_"$LINE1".bam \
	-L $LINE -U LENIENT_VCF_PROCESSING  -stand_call_conf 2 -stand_emit_conf 2 -mbq 10 -o /scratch/justin.zook/NA12878NCIIonWG/NCIIonWG_UG2.7_"$LINE1".vcf
