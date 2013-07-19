#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=3:00:00
#$ -l h_vmem=5G
#$ -N GATK_CombVarChromUG
#$ -pe thread 1
#$ -t 15

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


#java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:ILLCLIAsnp /scratch/justin.zook/NA12878/UGvcfs/CEPH_NA12878.gatk.conf2.snps.raw_"$LINE1".vcf \
 -V:ILLCLIAindel /scratch/justin.zook/NA12878/UGvcfs/CEPH_NA12878_b37_RG.sorted.gatk.conf2.indel.raw_"$LINE1".vcf \
 -V:HSWEx /scratch/justin.zook/NA12878/UGvcfs/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf2.mbq10.snpindel.LAlign.raw_"$LINE1".vcf \
 -V:HSWG /scratch/justin.zook/NA12878/UGvcfs/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.gatk.conf2.mbq10.snpindel.LAlign.raw_"$LINE1".vcf \
 -V:CG /scratch/justin.zook/NA12878/UGvcfs/CG_NA12878_all_nochr_sort_"$LINE1".vcf \
 -V:IllPCRFreesnp /scratch/justin.zook/NA12878/UGvcfs/ERR09157x_dedup_realign.gatk.conf2.snps.raw_"$LINE1".vcf \
 -V:IllPCRFreeindel /scratch/justin.zook/NA12878/UGvcfs/ERR09157x_dedup_realign.sorted.gatk.conf2.indel.raw_"$LINE1".vcf \
 -V:IonExindel /scratch/justin.zook/NA12878/UGvcfs/IonExome.gatk.conf2.indel.raw_"$LINE1".vcf \
 -V:IonExsnp /scratch/justin.zook/NA12878/UGvcfs/IonExome.gatk.conf2.snps.raw_"$LINE1".vcf \
 -V:XIllWGsnp /scratch/justin.zook/NA12878/UGvcfs/LP6005036-DNA_A01.gatk.conf2.snps.raw_"$LINE1".vcf \
 -V:XIllWGindel /scratch/justin.zook/NA12878/UGvcfs/LP6005036-DNA_A01_b37_RG.sorted.gatk.conf2.indel.raw_"$LINE1".vcf \
 -V:XFosIll /scratch/justin.zook/NA12878/UGvcfs/NA12878-7k-illumina-fosmid.gatk.conf2.mbq10.snpindel.raw_"$LINE1".vcf \
 -V:XFosSol /scratch/justin.zook/NA12878/UGvcfs/NA12878-7k-solid-rgfix.gatk.conf2.mbq10.snpindel.raw_"$LINE1".vcf \
 -V:ILLall /scratch/justin.zook/NA12878/UGvcfs/NA12878.all.ILLUMINA.bwa.CEU.sorted.gatk.conf2.mbq10.snpindel.LAlign.raw_"$LINE1".vcf \
 -V:454all /scratch/justin.zook/NA12878/UGvcfs/NA12878.all.LS454.ssaha2.CEU.sorted.gatk.conf2.mbq10.snpindel.LAlign.raw_"$LINE1".vcf \
 -V:XSol4WG /scratch/justin.zook/NA12878/UGvcfs/XPSolWGLS.gatk.conf2.mbq10.snpindel.raw_"$LINE1".vcf \
 -V:OMNI /scratch/justin.zook/NA12878/UGvcfs/Omni25_genotypes_NA12878.b37_SNPs_"$LINE1".vcf \
 -V:PlatGen /scratch/justin.zook/NA12878/PlatGen/PlatGen_UG2.5_"$LINE1".vcf \
 -priority ILLCLIAsnp,ILLCLIAindel,HSWG,HSWEx,CG,IllPCRFreesnp,IllPCRFreeindel,IonExindel,IonExsnp,XIllWGsnp,XIllWGindel,XFosIll,XFosSol,ILLall,454all,XSol4WG,OMNI,PlatGen \
 -o /scratch/justin.zook/NA12878/UGCombined_snpindel_130517_"$LINE1".vcf
 
#java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:ILLCLIA /scratch/justin.zook/NA12878/ILLCLIA/CEPH_NA12878_haplo2.5_"$LINE1".vcf \
 -V:HSWEx /scratch/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_haplo2.5_"$LINE1".vcf \
 -V:HSWG /scratch/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_haplo2.5_"$LINE1".vcf \
 -V:IllPCRFree /scratch/justin.zook/NA12878/IllPCRFree/ERR09157x_dedup_realign_haplo2.5_"$LINE1".vcf \
 -V:IonEx /scratch/justin.zook/NA12878/IonEx/IonExome_haplo2.5_"$LINE1".vcf \
 -V:XIllWG /scratch/justin.zook/NA12878/XIllWG/LP6005036-DNA_A01_haplo2.5_"$LINE1".vcf \
 -V:ILLWG /scratch/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_haplo2.5_"$LINE1".vcf \
 -V:ILLWEx /scratch/justin.zook/NA12878/ILLWEx/NA12878.all.WExS.ILLUMINA.bwa.CEU_haplo2.5_"$LINE1".vcf \
 -V:XSol4WG /scratch/justin.zook/NA12878/XSol4WG/NIST_XPrizeSOLiD_Lifescope_haplo2.5_"$LINE1".vcf \
 -V:PlatGen /scratch/justin.zook/NA12878/PlatGen/PlatGen_haplo2.5_"$LINE1".vcf \
 -priority ILLCLIA,HSWG,HSWEx,IllPCRFree,IonEx,XIllWG,ILLWG,ILLWEx,XSol4WG,PlatGen \
 -o /scratch/justin.zook/NA12878/HaploCombined_snpindel_130517_"$LINE1".vcf
# -V:XFosIll /scratch/justin.zook/NA12878/XFos/NA12878-7k-illumina-fosmid_haplo2.5_"$LINE1".vcf \
# -V:XFosSol /scratch/justin.zook/NA12878/XFos/NA12878-7k-solid-rgfix_haplo2.5_"$LINE1".vcf \
 
#java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:UG /scratch/justin.zook/NA12878/UGCombined_snpindel_130517_"$LINE1".vcf \
 -V:Haplo /scratch/justin.zook/NA12878/HaploCombined_snpindel_130517_"$LINE1".vcf \
 -o /scratch/justin.zook/NA12878/UGHaploCombined_snpindel_130517_"$LINE1".vcf
 
#java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/UGHaploCombined_snpindel_130517_"$LINE1".vcf \
-selectType INDEL -selectType MIXED -selectType MNP -selectType SYMBOLIC \
 -o /scratch/justin.zook/NA12878/UGHaploCombined_indel_130517_"$LINE1".vcf

#only for X,Y,MT
#cp /scratch/justin.zook/NA12878/UGCombined_snpindel_130517_"$LINE1".vcf \
  /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_"$LINE1".vcf -U LENIENT_VCF_PROCESSING
 
#For others
java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:UG /scratch/justin.zook/NA12878/UGCombined_snpindel_130517_"$LINE1".vcf \
 -V:Haplo /scratch/justin.zook/NA12878/HaploCombined_snpindel_130517_"$LINE1".vcf \
 -V:Cortexk55IllPCRFree /scratch/justin.zook/NA12878/cortex/vcfs/IllPCRFree_cortex_k55.vcf_union_BC_calls_k55.decomp.noMISMAPPEDorMAPQ_"$LINE1".vcf \
 -o /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_"$LINE1".vcf -U LENIENT_VCF_PROCESSING
 
java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T CombineVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
 -V:UG /scratch/justin.zook/NA12878/UGCombined_snpindel_130517_"$LINE1".vcf \
 -V:Haplo /scratch/justin.zook/NA12878/HaploCombined_snpindel_130517_"$LINE1".vcf \
 -V:Cortexk55IllPCRFree /scratch/justin.zook/NA12878/cortex/vcfs/IllPCRFree_cortex_k55.vcf_union_BC_calls_k55.decomp.noMISMAPPEDorMAPQorbadref_"$LINE1".vcf \
 -o /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_"$LINE1".vcf -U LENIENT_VCF_PROCESSING



 
java -jar -Xmx4g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/justin.zook/references/human_g1k_v37.fasta \
-V /scratch/justin.zook/NA12878/UGHaploCortexCombined_snpindel_130517_"$LINE1".vcf \
-selectType INDEL -selectType MIXED -selectType MNP -selectType SYMBOLIC \
 -o /scratch/justin.zook/NA12878/UGHaploCortexCombined_indel_130517_"$LINE1".vcf -U LENIENT_VCF_PROCESSING



