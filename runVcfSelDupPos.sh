#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=1:00:00
#$ -l h_vmem=1G
#$ -N runVcfSelDup
#$ -pe thread 1
#$ -t 17

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

 
echo "Running job $JOB_NAME, $JOB_ID on $HOSTNAME"
perl VcfSelDupPos.pl /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_dups_"$LINE1"_bad.vcf

#remove variants that have an unassembled event with another allele since GATK doesn't like these
grep -v ',<U' /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_dups_"$LINE1"_bad.vcf > /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_dups_"$LINE1".vcf

rm /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_dups_"$LINE1"_bad.vcf
