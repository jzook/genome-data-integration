#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=4:00:00
#$ -l h_vmem=10G
#$ -N RunVcfCombineUGHaplo
#$ -pe thread 1
#$ -t 1-21


BASE=/home/justin.zook
GATK=$BASE/GATK/GenomeAnalysisTK-2.6-2-gf57256b

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

    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/ILLCLIA/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/IonEx/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/XIllWG/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/ILLWEx/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/ILLWG/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/HSWEx/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/HSWG/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/IllPCRFree/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/454WG/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/XSol4WG/*_"$LINE1"*.gz
    #gunzip /projects/scratch-data-backup/justin.zook/NA12878/CG/*_"$LINE1"*.gz

if [ 0 -gt 1 ]; then

    ln /projects/scratch-data-backup/justin.zook/NA12878/ILLCLIA/CEPH_NA12878_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ILLCLIAcall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/IonEx/IonExome_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_IonExcall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/XIllWG/LP6005036-DNA_A01_b37_RG_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_XIllcall_UGdups_"$LINE1".vcf
  
    ln /projects/scratch-data-backup/justin.zook/NA12878/ILLWEx/NA12878.all.WExS.ILLUMINA.bwa.CEU_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ILLWExcall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_HSWExcall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_HSWGcall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/IllPCRFree/ERR09157x_dedup_realign_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_IllPCRFreecall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/454WG/NA12878.all.WGS.LS454.ssaha2.CEU_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_454WGcall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/XSol4WG/NIST_XPrizeSOLiD_Lifescope_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_XPSolWGLScall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ILLWGcall_UGdups_"$LINE1".vcf

    ln /projects/scratch-data-backup/justin.zook/NA12878/PlatGen/PlatGen_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_PlatGencall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/MiSeq250/MiSeq250_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ILL250call_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/NCIIonWG/NCIIonWG_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_NCIIonWGcall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/CG/GS00362-DNA_C04_all_b37_UGrecall130719_dups_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_CGcall_UGdups_"$LINE1".vcf
    
    ln /projects/scratch-data-backup/justin.zook/NA12878/CG/GS00362-DNA_C04_all_b37_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_CGcall_UGHapMerge_"$LINE1".vcf
     
     
	 perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/ILLCLIA/CEPH_NA12878_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/ILLCLIA/CEPH_NA12878_b37_RG_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ILLCLIAcall_UGHapMerge_"$LINE1".vcf
	  
	 
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/IonEx/IonExome_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/IonEx/IonExome_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_IonExcall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/XIllWG/LP6005036-DNA_A01_b37_RG_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/XIllWG/LP6005036-DNA_A01_b37_RG_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_XIllcall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/ILLWEx/NA12878.all.WExS.ILLUMINA.bwa.CEU_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/ILLWEx/NA12878.all.WExS.ILLUMINA.bwa.CEU_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ILLWExcall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/ILLWG/NA12878.all.WGS.ILLUMINA.bwa.CEU_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ILLWGcall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_HSWExcall_UGHapMerge_"$LINE1".vcf
#    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130719_1a.vcf /projects/scratch-data-backup/justin.zook/NA12878/HSWEx/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117_haplotrigger2.6_130719_1a.vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_HSWExcall_UGHapMerge_1a.vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/HSWG/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_HSWGcall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/IllPCRFree/ERR09157x_dedup_realign_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/IllPCRFree/ERR09157x_dedup_realign_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_IllPCRFreecall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/454WG/NA12878.all.WGS.LS454.ssaha2.CEU_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/454WG/NA12878.all.WGS.LS454.ssaha2.CEU_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_454WGcall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/XSol4WG/NIST_XPrizeSOLiD_Lifescope_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/XSol4WG/NIST_XPrizeSOLiD_Lifescope_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_XPSolWGLScall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/PlatGen/PlatGen_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/PlatGen/PlatGen_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_PlatGencall_UGHapMerge_"$LINE1".vcf
    
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/MiSeq250/MiSeq250_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/MiSeq250/MiSeq250_haplotrigger2.6_130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ILL250call_UGHapMerge_"$LINE1".vcf
  
    perl ~/scripts/VcfCombineUGHaplo_v0.4_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/NCIIonWG/NCIIonWG_UGrecall130719_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/NCIIonWG/NCIIonWG_haplotrigger2.7_130909_"$LINE1".vcf /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_NCIIonWGcall_UGHapMerge_"$LINE1".vcf
    
fi  
   perl ~/scripts/VcfHighConfUGHaploMulti_HomJoint_1.2_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_ "$LINE1" HSWG ILLWG XIll 454WG ILLCLIA IllPCRFree XPSolWGLS IonEx HSWEx ILLWEx CG PlatGen ILL250 NCIIonWG
#exit    
    gzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_*Hom*_"$LINE1".vcf
    
    gzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_*Het*_"$LINE1".vcf
    
#run this after all have finished to combine csvs    
#cat /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_1a.csv /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_1b.csv /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_2a.csv /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_2b.csv | gzip -c > /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_1_2.csv.gz     


