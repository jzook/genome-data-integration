#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=4:00:00
#$ -l h_vmem=5G
#$ -N GATK_CombVarChromUG
#$ -pe thread 2
#$ -t 4

BASE=/home/justin.zook
GATK=$BASE/GATK/GenomeAnalysisTK-2.6-4-g3e5ff60

export  java_HOME=/home/justin.zook/jre1.7.0_25
export PATH=$java_HOME/bin:$PATH

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


#  java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:AllUG130719 /scratch/justin.zook/NA12878/UGCombined_snpindel_130517_"$LINE1".vcf \
 -V:MS250UG /scratch/justin.zook/NA12878/MiSeq250/MiSeq250_UG2.5_"$LINE1".vcf \
 -o /scratch/justin.zook/NA12878/UGCombined_snpindel_130719_"$LINE1".vcf -U LENIENT_VCF_PROCESSING
 
  
#  java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:ILLCLIA /scratch/justin.zook/NA12878/ILLCLIA/CEPH_NA12878_haplo2.6_"$LINE1".vcf \
 -V:HSWEx /scratch/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_haplo2.6_"$LINE1".vcf \
 -V:HSWG /scratch/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_haplo2.6_"$LINE1".vcf \
 -V:IllPCRFree /scratch/justin.zook/NA12878/IllPCRFree/ERR09157x_dedup_realign_haplo2.6_"$LINE1".vcf \
 -V:IonEx /scratch/justin.zook/NA12878/IonEx/IonExome_haplo2.6_"$LINE1".vcf \
 -V:XIllWG /scratch/justin.zook/NA12878/XIllWG/LP6005036-DNA_A01_haplo2.6_"$LINE1".vcf \
 -V:ILLWG /scratch/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_haplo2.6_"$LINE1".vcf \
 -V:ILLWEx /scratch/justin.zook/NA12878/ILLWEx/NA12878.all.WExS.ILLUMINA.bwa.CEU_haplo2.6_"$LINE1".vcf \
 -V:XSol4WG /scratch/justin.zook/NA12878/XSol4WG/NIST_XPrizeSOLiD_Lifescope_haplo2.6_"$LINE1".vcf \
 -V:PlatGen /scratch/justin.zook/NA12878/PlatGen/PlatGen_haplo2.6_"$LINE1".vcf \
 -V:MS250 /scratch/justin.zook/NA12878/MiSeq250/MiSeq250_haplo2.6_"$LINE1".vcf \
 -priority ILLCLIA,HSWG,HSWEx,IllPCRFree,IonEx,XIllWG,ILLWG,ILLWEx,XSol4WG,PlatGen,MS250 \
 -o /scratch/justin.zook/NA12878/HaploCombined_snpindel_130719_"$LINE1".vcf -U LENIENT_VCF_PROCESSING
 
#  java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:UG /scratch/justin.zook/NA12878/UGCombined_snpindel_130719_"$LINE1".vcf \
 -V:Haplo /scratch/justin.zook/NA12878/HaploCombined_snpindel_130719_"$LINE1".vcf \
 -o /scratch/justin.zook/NA12878/UGHaploCombined_snpindel_130719_"$LINE1".vcf -U LENIENT_VCF_PROCESSING
 
#  java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/UGHaploCombined_snpindel_130719_"$LINE1".vcf \
-selectType INDEL -selectType MIXED -selectType MNP -selectType SYMBOLIC \
 -o /scratch/justin.zook/NA12878/UGHaploCombined_indel_130719_"$LINE1".vcf  -U LENIENT_VCF_PROCESSING

#only for X,Y,MT
#cp /scratch/justin.zook/NA12878/HaploCombined_snpindel_130719_"$LINE1".vcf \
  /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130719_"$LINE1".vcf 
 

#For others
# java -jar -Xmx4g $BASE/GATK/bcbio.variation-0.0.9-standalone.jar variant-prep /scratch/justin.zook/NA12878/cortex/vcfs/IllPCRFree_cortex_k55.vcf_union_BC_calls_k55.decomp.noMISMAPPEDorMAPQ_"$LINE1".vcf /scratch/justin.zook/references/human_g1k_v37.fasta 
  java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:UG /scratch/justin.zook/NA12878/UGCombined_snpindel_130719_"$LINE1".vcf \
 -V:Haplo /scratch/justin.zook/NA12878/HaploCombined_snpindel_130719_"$LINE1".vcf \
 -V:Cortexk55IllPCRFree /scratch/justin.zook/NA12878/cortex/vcfs/IllPCRFree_cortex_k55.vcf_union_BC_calls_k55.decomp.noMISMAPPEDorMAPQ_"$LINE1"-fullprep.vcf \
 -o /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130719_"$LINE1".vcf -U LENIENT_VCF_PROCESSING
 
 # java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:UG /scratch/justin.zook/NA12878/UGCombined_snpindel_130719_"$LINE1".vcf \
 -V:Haplo /scratch/justin.zook/NA12878/HaploCombined_snpindel_130719_"$LINE1".vcf \
 -V:Cortexk55IllPCRFree /scratch/justin.zook/NA12878/cortex/vcfs/IllPCRFree_cortex_k55.vcf_union_BC_calls_k55.decomp.noMISMAPPEDorMAPQorbadref_"$LINE1".vcf \
 -o /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130719_"$LINE1".vcf -U LENIENT_VCF_PROCESSING

 java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130719_"$LINE1".vcf \
-selectType INDEL -selectType MIXED -selectType MNP -selectType SYMBOLIC \
 -o /scratch/justin.zook/NA12878/UGHaploCortexCombined_indel_130719_"$LINE1".vcf  -U LENIENT_VCF_PROCESSING

