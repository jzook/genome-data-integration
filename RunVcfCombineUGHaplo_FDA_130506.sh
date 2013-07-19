#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=2:00:00
#$ -l h_vmem=10G
#$ -N RunVcfCombineUGHaplo
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


LINE=`head -n $SGE_TASK_ID $TEMP_LIST | tail -n 1`
LINE1=`head -n $SGE_TASK_ID $TEMP_LIST1 | tail -n 1`


# export PATH=/home/mikem/GATK/GenomeAnalysisTK-2.1-11-g13c0244/resources:$PATH
 
echo "Running job $JOB_NAME, $JOB_ID on $HOSTNAME"

    gunzip /scratch/justin.zook/NA12878/ILLCLIA/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/IonEx/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/XIllWG/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/ILLWEx/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/ILLWG/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/HSWEx/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/HSWG/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/IllPCRFree/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/454WG/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/XSol4WG/*_"$LINE1"*.gz
    gunzip /scratch/justin.zook/NA12878/CG/*_"$LINE1"*.gz

    ln /scratch/justin.zook/NA12878/ILLCLIA/CEPH_NA12878_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_ILLCLIAcall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/IonEx/IonExome_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_IonExcall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/XIllWG/LP6005036-DNA_A01_b37_RG_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_XIllcall_UGdups_"$LINE1".vcf
  
    ln /scratch/justin.zook/NA12878/ILLWEx/NA12878.all.WExS.ILLUMINA.bwa.CEU_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_ILLWExcall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_HSWExcall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_HSWGcall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/IllPCRFree/ERR09157x_dedup_realign_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_IllPCRFreecall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/454WG/NA12878.all.WGS.LS454.ssaha2.CEU_UGrecall130402_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_454WGcall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/XSol4WG/NIST_XPrizeSOLiD_Lifescope_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_XPSolWGLScall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_ILLWGcall_UGdups_"$LINE1".vcf

    ln /scratch/justin.zook/NA12878/PlatGen/PlatGen_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_PlatGencall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/CG/GS00362-DNA_C04_all_b37_UGrecall130517_dups_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_CGcall_UGdups_"$LINE1".vcf
    
    ln /scratch/justin.zook/NA12878/CG/GS00362-DNA_C04_all_b37_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_CGcall_UGHapMerge_"$LINE1".vcf
    
mv /scratch/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_haplotrigger2.4_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_haplotrigger2.5_130517_"$LINE1".vcf
 
     
	 #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/ILLCLIA/CEPH_NA12878_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/ILLCLIA/CEPH_NA12878_b37_RG_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_ILLCLIAcall_UGHapMerge_"$LINE1".vcf
	  
	 
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/IonEx/IonExome_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/IonEx/IonExome_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_IonExcall_UGHapMerge_"$LINE1".vcf
    
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/XIllWG/LP6005036-DNA_A01_b37_RG_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/XIllWG/LP6005036-DNA_A01_b37_RG_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_XIllcall_UGHapMerge_"$LINE1".vcf
    
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/ILLWEx/NA12878.all.WExS.ILLUMINA.bwa.CEU_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/ILLWEx/NA12878.all.WExS.ILLUMINA.bwa.CEU_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_ILLWExcall_UGHapMerge_"$LINE1".vcf
    
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_ILLWGcall_UGHapMerge_"$LINE1".vcf
    
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_HSWExcall_UGHapMerge_"$LINE1".vcf
    
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_HSWGcall_UGHapMerge_"$LINE1".vcf
    
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/IllPCRFree/ERR09157x_dedup_realign_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/IllPCRFree/ERR09157x_dedup_realign_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_IllPCRFreecall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/454WG/NA12878.all.WGS.LS454.ssaha2.CEU_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/454WG/NA12878.all.WGS.LS454.ssaha2.CEU_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_454WGcall_UGHapMerge_"$LINE1".vcf
    
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/XSol4WG/NIST_XPrizeSOLiD_Lifescope_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/XSol4WG/NIST_XPrizeSOLiD_Lifescope_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_XPSolWGLScall_UGHapMerge_"$LINE1".vcf
    
    #perl ~/scripts/VcfCombineUGHaplo_v0.3.pl /scratch/justin.zook/NA12878/PlatGen/PlatGen_UGrecall130517_"$LINE1".vcf /scratch/justin.zook/NA12878/PlatGen/PlatGen_haplotrigger2.5_130517_"$LINE1".vcf /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_PlatGencall_UGHapMerge_"$LINE1".vcf
    
   perl ~/scripts/VcfHighConfUGHaploMulti_HomJoint_1.2_FDA.pl /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_ "$LINE1" HSWG ILLWG XIll 454WG ILLCLIA IllPCRFree XPSolWGLS IonEx HSWEx ILLWEx CG PlatGen
    
    gzip -f /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_*Hom*_"$LINE1".vcf
    
    gzip -f /scratch/justin.zook/NA12878/Integration130517/AllFDAdatasets_130517_*Het*_"$LINE1".vcf
    
    



