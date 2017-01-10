#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=24:00:00
#$ -l h_vmem=11G
#$ -N GATK_VQSR_Het
#$ -pe thread 1
#$ -t 2-14


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

#    export PERL5LIB=/usr/bin/perl
    perl ~/scripts/runVariantRecalHetRef_UGHapMerge_FDA_130523.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103 $LINE all $SGE_TASK_ID 
    
    



