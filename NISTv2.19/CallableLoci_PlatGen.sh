#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=5:00:00
#$ -l h_vmem=5G
#$ -N CallableLoci_PlatGen
#$ -pe thread 1
#$ -t 1-21



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

 
echo "Running job $JOB_NAME, $JOB_ID on $HOSTNAME"

#java -jar -Xmx3g  $GATK/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T CallableLoci \
#	-I  /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".bam \
#	 -summary /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".callable.summary \
#	 -o  /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed \
#	 -L $LINE -U LENIENT_VCF_PROCESSING -frlmq 0.1 -mlmq 10 -mbq 10 -minDepth 5 --minDepthForLowMAPQ 5 -mmq 10
        
 grep CALLABLE /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | cut -f 1-3 > /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed 
 #  sed 's/ /'$'\t''/g' /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed > /scratch/justin.zook/NA12878/PlatGen/PlatGen_realign_"$LINE1".frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable_tab.bed
   
    
    



