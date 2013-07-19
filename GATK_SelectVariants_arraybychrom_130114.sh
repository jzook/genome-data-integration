#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=1:00:00
#$ -l h_vmem=11G
#$ -N GATK_SelVarChromUG
#$ -pe thread 1
#$ -t 2-20 

BASE=/home/justin.zook
GATK=$BASE/GATK/GenomeAnalysisTK-2.3-4-g57ea19f

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


java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/CEPH_NA12878.gatk.conf2.snps.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/CEPH_NA12878.gatk.conf2.snps.raw_"$LINE1".vcf

##java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/CEPH_NA12878.realigned-reorder.gatk.conf2.indel.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/CEPH_NA12878.realigned-reorder.gatk.conf2.indel.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/CEPH_NA12878_b37_RG.sorted.gatk.conf2.indel.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/CEPH_NA12878_b37_RG.sorted.gatk.conf2.indel.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf2.mbq10.snpindel.LAlign.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf2.mbq10.snpindel.LAlign.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf2.mbq10.snpindel.LAlign.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf2.mbq10.snpindel.LAlign.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/CG_NA12878_all_nochr_sort.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/CG_NA12878_all_nochr_sort_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/ERR09157x_dedup_realign.gatk.conf2.snps.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/ERR09157x_dedup_realign.gatk.conf2.snps.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/ERR09157x_dedup_realign.sorted.gatk.conf2.indel.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/ERR09157x_dedup_realign.sorted.gatk.conf2.indel.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/IonExome.gatk.conf2.indel.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/IonExome.gatk.conf2.indel.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/IonExome.gatk.conf2.snps.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/IonExome.gatk.conf2.snps.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/LP6005036-DNA_A01.gatk.conf2.snps.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/LP6005036-DNA_A01.gatk.conf2.snps.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/LP6005036-DNA_A01_b37_RG.sorted.gatk.conf2.indel.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/LP6005036-DNA_A01_b37_RG.sorted.gatk.conf2.indel.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/NA12878-7k-illumina-fosmid.gatk.conf2.mbq10.snpindel.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/NA12878-7k-illumina-fosmid.gatk.conf2.mbq10.snpindel.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/NA12878-7k-solid-rgfix.gatk.conf2.mbq10.snpindel.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/NA12878-7k-solid-rgfix.gatk.conf2.mbq10.snpindel.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/NA12878.all.ILLUMINA.bwa.CEU.sorted.gatk.conf2.mbq10.snpindel.LAlign.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/NA12878.all.ILLUMINA.bwa.CEU.sorted.gatk.conf2.mbq10.snpindel.LAlign.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/NA12878.all.LS454.ssaha2.CEU.sorted.gatk.conf2.mbq10.snpindel.LAlign.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/NA12878.all.LS454.ssaha2.CEU.sorted.gatk.conf2.mbq10.snpindel.LAlign.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/XPSolWGLS.gatk.conf2.mbq10.snpindel.raw.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/XPSolWGLS.gatk.conf2.mbq10.snpindel.raw_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/UGvcfs/Omni25_genotypes_NA12878.b37_SNPs.vcf $LINE -o /scratch/justin.zook/NA12878/UGvcfs/Omni25_genotypes_NA12878.b37_SNPs_"$LINE1".vcf

java -jar -Xmx10g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta -V /scratch/justin.zook/NA12878/cortex/vcfs/IllPCRFree_cortex_k55.vcf_union_BC_calls_k55.decomp.noMISMAPPEDorMAPQ.sorted.vcf $LINE -o /scratch/justin.zook/NA12878/cortex/vcfs/IllPCRFree_cortex_k55.vcf_union_BC_calls_k55.decomp.noMISMAPPEDorMAPQ_"$LINE1".vcf -U LENIENT_VCF_PROCESSING


