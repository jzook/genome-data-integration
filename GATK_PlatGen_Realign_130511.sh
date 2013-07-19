#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=250:00:00
#$ -l h_vmem=180G
#$ -N Realign_PlatGen
#$ -pe thread 1
#$ -t 12-13

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


LINE=`head -n $SGE_TASK_ID $TEMP_LIST | tail -n 1`
LINE1=`head -n $SGE_TASK_ID $TEMP_LIST1 | tail -n 1`


# export PATH=/home/mikem/GATK/GenomeAnalysisTK-2.1-11-g13c0244/resources:$PATH
 
echo "Running job $JOB_NAME, $JOB_ID on $HOSTNAME"


#java -jar -Xmx21g $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $BASE/references/hs37d5.fa \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174324_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174325_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174326_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174327_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174328_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174329_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174330_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174331_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174332_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174333_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174334_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174335_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174336_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174337_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174338_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174339_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174340_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174341_bwa_b37_RG.sorted.bam \
	-L $LINE -o /scratch/justin.zook/NA12878/PlatGen/PlatGen_"$LINE1".intervals

java -jar -Xmx21g $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $BASE/references/hs37d5.fa \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174324_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174325_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174326_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174327_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174328_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174329_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174330_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174331_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174332_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174333_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174334_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174335_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174336_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174337_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174338_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174339_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174340_bwa_b37_RG.sorted.bam \
	-I /scratch/justin.zook/NA12878/PlatGen/ERR174341_bwa_b37_RG.sorted.bam \
	-L $LINE -maxInMemory 500000 -maxReads 80000 -targetIntervals /scratch/justin.zook/NA12878/PlatGen/PlatGen_"$LINE1".intervals -o /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".bam



