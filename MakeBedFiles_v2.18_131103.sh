#$ -cwd
#$ -S /bin/sh
#$ -j y
#$ -l h_rt=7:00:00
#$ -l h_vmem=20G
#$ -N MakeBedFiles_v2.18
#$ -pe thread 1
#$ -t 1



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


#sleep 1h

 
echo "Running job $JOB_NAME, $JOB_ID on $HOSTNAME"

if [ 0 -gt 1 ]; then
    perl ~/scripts/CombineChromVcf_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_VQSRv2.18_2mindatasets_5minYesNoRatio 
    perl ~/scripts/CombineChromVcf_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio  
    perl ~/scripts/CombineChromVcf_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio 
    perl ~/scripts/CombineChromVcf_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_Het_VQSRv2.18_2mindatasets_5minYesNoRatio 
    perl ~/scripts/CombineChromVcf_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomVar_VQSRv2.18_2mindatasets_5minYesNoRatio 
    perl ~/scripts/CombineChromVcf_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio 
   
 ~/freebayes/freebayes_0.9.9.2/vcflib/vcfallelicprimitives /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf  | perl ~/scripts/vcfsorter.pl /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.dict - > /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives.vcf 
 ~/freebayes/freebayes_0.9.9.2/vcflib/vcfallelicprimitives /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf  | perl ~/scripts/vcfsorter.pl /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.dict - > /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives.vcf 
 perl ~/scripts/VcfRemoveDup.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives 
 perl ~/scripts/VcfRemoveDup.pl /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives 
#exit
#fi
    # java -jar ~/GATK/bcbio.variation-0.1.0-SNAPSHOT-standalone.jar variant-prep --keep-ref /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta
    # java -jar ~/GATK/bcbio.variation-0.1.0-SNAPSHOT-standalone.jar variant-prep --keep-ref /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta
    # java -jar ~/GATK/bcbio.variation-0.1.0-SNAPSHOT-standalone.jar variant-prep --keep-ref /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_Het_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta
    # java -jar ~/GATK/bcbio.variation-0.1.0-SNAPSHOT-standalone.jar variant-prep --keep-ref /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomVar_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta
    # java -jar ~/GATK/bcbio.variation-0.1.0-SNAPSHOT-standalone.jar variant-prep --keep-ref /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta
    # java -jar ~/GATK/bcbio.variation-0.1.0-SNAPSHOT-standalone.jar variant-prep --keep-ref /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta


 java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf \
-U LENIENT_VCF_PROCESSING -selectType SNP \
 -o /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all_snp.vcf

 java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf \
-U LENIENT_VCF_PROCESSING -selectType SNP \
 -o /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all_snp.vcf


grep '^#\|PASS' /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all_snp.vcf > /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all_snp_PASS.vcf 

#add new dataset's callableloci if necessary
#grep CALLABLE /Volumes/SSD960/bed/GS00362-DNA_C04_all_b37.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/GS00362-DNA_C04_all_b37.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/LP6005036-DNA_A01_b37_RG.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/LP6005036-DNA_A01_b37_RG.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/XPSolWG_combined.bam.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/XPSolWG_combined.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/Ionexome.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/Ionexome.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/ERR09157x_dedup_realign.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/ERR09157x_dedup_realign.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/CEPH_NA12878_b37_RG.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/CEPH_NA12878_b37_RG.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/NA12878.all.WExS.ILLUMINA.bwa.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/NA12878.all.WExS.ILLUMINA.bwa.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/NA12878.all.WGS.ILLUMINA.bwa.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/NA12878.all.WGS.ILLUMINA.bwa.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/NA12878.all.WGS.LS454.ssaha2.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/NA12878.all.WGS.LS454.ssaha2.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &
#grep CALLABLE /Volumes/SSD960/bed/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callable.bed | grep -v GL | sed 's/ /'$'\t''/g' | cut -f 1-3 > /Volumes/SSD960/bed/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed &

#combine bed files
#cat /Volumes/SSD960/bed/GS00362-DNA_C04_all_b37.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/LP6005036-DNA_A01_b37_RG.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/XPSolWG_combined.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/Ionexome.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/ERR09157x_dedup_realign.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/CEPH_NA12878_b37_RG.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/NA12878.all.WExS.ILLUMINA.bwa.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/NA12878.all.WGS.ILLUMINA.bwa.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/NA12878.all.WGS.LS454.ssaha2.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /Volumes/SSD960/bed/PlatGen_realign_all.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed > /Volumes/SSD960/bed/union13callableonly.bed

#sort -k 1,1 -k2,2 -k3,3 -n /Volumes/SSD960/bed/union13callableonly.bed > /Volumes/SSD960/bed/union13callableonly_sort.bed
#mergeBed -i /Volumes/SSD960/bed/union13callableonly_sort.bed > /Volumes/SSD960/bed/union13callableonlymerged.bed

#grep -nm 1 '^1' /Volumes/SSD960/bed/union13callableonlymerged.bed
#594604:1	10001	10053
#wc -l /Volumes/SSD960/bed/union13callableonlymerged.bed
#1467073 /Volumes/SSD960/bed/union13callableonlymerged.bed

#head -594603 /Volumes/SSD960/bed/union13callableonlymerged.bed > /Volumes/SSD960/bed/union13callableonly_XYMT.bed
#tail -872470 /Volumes/SSD960/bed/union13callableonlymerged.bed > /Volumes/SSD960/bed/union13callableonlymerged_1-22.bed
#sort -k 1,1 -k2,2n -k3,3n /Volumes/SSD960/bed/union13callableonly_XYMT.bed > /Volumes/SSD960/bed/union13callableonly_XYMT_sort.bed
#mergeBed -i /Volumes/SSD960/bed/union13callableonly_XYMT_sort.bed > /Volumes/SSD960/bed/union13callableonlymerged_XYMT_sort.bed
#cat /Volumes/SSD960/bed/union13callableonlymerged_1-22.bed /Volumes/SSD960/bed/union13callableonlymerged_XYMT_sort.bed > /Volumes/SSD960/bed/union13callableonlymerged_all_sort.bed

#make new union of all sites covered by at least 3 datasets
#perl $BASE/scripts/CombineChromVcf_nohead_FDA.pl /projects/scratch-data-backup/justin.zook/NA12878/MiSeq250/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130520n .frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly
#cat /projects/scratch-data-backup/justin.zook/NA12878/bed/CEPH_NA12878_b37_RG.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.20120117.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/ERR09157x_dedup_realign.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/Ionexome.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/LP6005036-DNA_A01_b37_RG.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/NA12878.all.WExS.ILLUMINA.bwa.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/NA12878.all.WGS.ILLUMINA.bwa.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/NA12878.all.WGS.LS454.ssaha2.CEU.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/GS00362-DNA_C04_all_b37.sorted.bam.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/bed/PlatGen_realign_all.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed /projects/scratch-data-backup/justin.zook/NA12878/MiSeq250/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130520n_all.frlmq01.mlmq10.mbq10.mdflmq5.dp5.mmq10.callableonly.bed | cut -f 1-3 | grep -v '^Y' | sort -k 1,1 -k2,2 -k3,3 -n  > /projects/scratch-data-backup/justin.zook/NA12878/bed/unionall13.bed
# $BASE/bedtools-2.17.0/bin/genomeCoverageBed -i /projects/scratch-data-backup/justin.zook/NA12878/bed/unionall13.bed -bg -g $BASE/bedtools-2.17.0/genomes/human.b37.genome > /projects/scratch-data-backup/justin.zook/NA12878/bed/unionall13.cov.bed
#awk '$4 > 2' /projects/scratch-data-backup/justin.zook/NA12878/bed/unionall13.cov.bed | cut -f 1-3 | $BASE/bedtools-2.17.0/bin/mergeBed -i stdin > /projects/scratch-data-backup/justin.zook/NA12878/bed/unionall13.cov3.bed



# java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta -V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf -U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/unionall13.cov3.bed  | grep -v '^#' | wc -l > /projects/scratch-data-backup/justin.zook/NA12878/bed/unionall13.cov3_v2.18_2mindatasets_5minYesNoRatio_countvar.txt & 

#add in sites called by consensus method
 $BASE/bedtools-2.17.0/bin/complementBed -i /projects/scratch-data-backup/justin.zook/NA12878/bed/unionall13.cov3.bed -g $BASE/bedtools-2.17.0/genomes/human.b37.genome > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_comp_v2.18_2mindatasets_5minYesNoRatio.bed
$BASE/bedtools-2.17.0/bin/subtractBed -a /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_comp_v2.18_2mindatasets_5minYesNoRatio.bed -b /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all_snp_PASS.vcf | $BASE/bedtools-2.17.0/bin/subtractBed -a stdin -b /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all_snp.vcf > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_comp_addcert_v2.18_2mindatasets_5minYesNoRatio.bed
$BASE/bedtools-2.17.0/bin/complementBed -i /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_comp_addcert_v2.18_2mindatasets_5minYesNoRatio.bed -g $BASE/bedtools-2.17.0/genomes/human.b37.genome | grep -v '^Y' > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_v2.18_2mindatasets_5minYesNoRatio.bed



 java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf \
-U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_v2.18_2mindatasets_5minYesNoRatio.bed \
 | grep -v '^#' | wc -l > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_v2.18_2mindatasets_5minYesNoRatio_countvar.txt &

grep -v 'PASS' /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf > /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered.vcf 

#add 50bp flanking regions to Het and HomVar uncertain calls
grep -v '0/0' /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered.vcf > /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered_noHomRef.vcf 

#fi
$BASE/bedtools-2.17.0/bin/subtractBed -a ~/bedtools-2.17.0/genomes/human.b37.genome.bed -b /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered_noHomRef.vcf | ~/bedtools-2.17.0/bin/complementBed -i stdin -g ~/bedtools-2.17.0/genomes/human.b37.genome | ~/bedtools-2.17.0/bin/flankBed -i stdin -b 50 -g ~/bedtools-2.17.0/genomes/human.b37.genome > /projects/scratch-data-backup/justin.zook/NA12878/bed/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered_noHomRef_flank50.bed

$BASE/bedtools-2.17.0/bin/subtractBed -a /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_v2.18_2mindatasets_5minYesNoRatio.bed -b /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered.vcf | $BASE/bedtools-2.17.0/bin/subtractBed -a stdin -b /projects/scratch-data-backup/justin.zook/NA12878/bed/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered_noHomRef_flank50.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_v2.18_2mindatasets_5minYesNoRatio.bed



# java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf \
-U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_v2.18_2mindatasets_5minYesNoRatio.bed \
 | grep -v '^#' | wc -l > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_comp_addcert_nouncert_v2.18_2mindatasets_5minYesNoRatio_countvar.txt &

#subtract simple repeat regions that aren't covered completely by >4 reads in at least one dataset
$BASE/bedtools-2.17.0/bin/subtractBed -a /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_v2.18_2mindatasets_5minYesNoRatio.bed -b /projects/scratch-data-backup/justin.zook/NA12878/bed/All.simplerepeatsnocov.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_v2.18_2mindatasets_5minYesNoRatio.bed

# java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf \
-U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_v2.18_2mindatasets_5minYesNoRatio.bed \
 | grep -v '^#' | wc -l > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_v2.18_2mindatasets_5minYesNoRatio_countvar.txt &

#subtract known seg dups
$BASE/bedtools-2.17.0/bin/subtractBed -a /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_v2.18_2mindatasets_5minYesNoRatio.bed -b /projects/scratch-data-backup/justin.zook/NA12878/bed/superdupsmerged_all_sort.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_v2.18_2mindatasets_5minYesNoRatio.bed

# java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf \
-U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_v2.18_2mindatasets_5minYesNoRatio.bed \
 | grep -v '^#' | wc -l > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_v2.18_2mindatasets_5minYesNoRatio_countvar.txt &

#subtract regions manually determined to be in decoy
$BASE/bedtools-2.17.0/bin/subtractBed -a /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_v2.18_2mindatasets_5minYesNoRatio.bed -b /projects/scratch-data-backup/justin.zook/NA12878/bed/decoy.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_v2.18_2mindatasets_5minYesNoRatio.bed

# java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf \
-U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_v2.18_2mindatasets_5minYesNoRatio.bed \
 | grep -v '^#' | wc -l > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_v2.18_2mindatasets_5minYesNoRatio_countvar.txt &


#subtract regions within 10bp of RepeatSeq STRs
$BASE/bedtools-2.17.0/bin/subtractBed -a /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_v2.18_2mindatasets_5minYesNoRatio.bed -b /projects/scratch-data-backup/justin.zook/NA12878/bed/hg19.max5.regions.pad10.corr.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_v2.18_2mindatasets_5minYesNoRatio.bed

 java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives_nodup.vcf \
-U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_v2.18_2mindatasets_5minYesNoRatio.bed \
 -o /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs.vcf
 ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs.vcf

#subtract regions in CNVs
$BASE/bedtools-2.17.0/bin/subtractBed -a /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_v2.18_2mindatasets_5minYesNoRatio.bed -b /projects/scratch-data-backup/justin.zook/NA12878/bed/dbvar_NA12878_CNVs_merged_all_sort.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed

 java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives_nodup.vcf \
-U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed \
-o /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf
  ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf

grep -v 'bias=none' /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf > /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all_bias.vcf

 java -jar -Xmx3g $GATK/GenomeAnalysisTK.jar -T SelectVariants -R /projects/scratch-data-backup/justin.zook/references/human_g1k_v37.fasta \
-V /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all_bias.vcf \
-U LENIENT_VCF_PROCESSING -L /projects/scratch-data-backup/justin.zook/NA12878/bed/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed \
-o /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all_bias_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf
  ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all_bias_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf

#make vcf containing all het and hom sites supported by 2+ platforms, without considering "confidence" otherwise
~/freebayes/freebayes_0.9.9.2/vcflib/vcffilter -f "platforms > 1" /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives_nodup.vcf | grep -v '0/0\|0|0' | sed 's/Uncertain/PASS/' > /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_2platforms_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives.vcf 

  ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_2platforms_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives.vcf
  ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives_nodup.vcf
  ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomRef_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf
  ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_Het_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf
  ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HomVar_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf
  ~/tabix-0.2.6/bgzip -f /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_VQSRv2.18_2mindatasets_5minYesNoRatio_all.vcf

fi

#annotate all variants with bed files
~/tabix-0.2.6/tabix /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_2platforms_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives.vcf.gz

sort -k 1,1 -k2,2n -k3,3n /projects/scratch-data-backup/justin.zook/NA12878/bed/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered_noHomRef_flank50.bed | ~/bedtools-2.17.0/bin/mergeBed -i stdin | awk '{$4="Complex"; print $0}' OFS=$'\t' > /projects/scratch-data-backup/justin.zook/NA12878/bed/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered_noHomRef_flank50.annotate.bed
awk '{$4="SimpleRepeat"; print $0}' OFS=$'\t' /projects/scratch-data-backup/justin.zook/NA12878/bed/All.simplerepeatsnocov.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/All.simplerepeatsnocov.annotate.bed
awk '{$4="SegDups"; print $0}' OFS=$'\t' /projects/scratch-data-backup/justin.zook/NA12878/bed/superdupsmerged_all_sort.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/superdupsmerged_all_sort.annotate.bed
awk '{$4="Decoy"; print $0}' OFS=$'\t' /projects/scratch-data-backup/justin.zook/NA12878/bed/decoy.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/decoy.annotate.bed
awk '{$4="RepeatSeq"; print $0}' OFS=$'\t' /projects/scratch-data-backup/justin.zook/NA12878/bed/hg19.max5.regions.pad10.corr.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/hg19.max5.regions.pad10.corr.annotate.bed
awk '{$4="dbVarSV"; print $0}' OFS=$'\t' /projects/scratch-data-backup/justin.zook/NA12878/bed/dbvar_NA12878_CNVs_merged_all_sort.bed > /projects/scratch-data-backup/justin.zook/NA12878/bed/dbvar_NA12878_CNVs_merged_all_sort.annotate.bed 

~/freebayes/freebayes_0.9.9.2/vcflib/vcfannotate -k filter -b /projects/scratch-data-backup/justin.zook/NA12878/bed/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_VQSRv2.18_2mindatasets_5minYesNoRatio_all_filtered_noHomRef_flank50.annotate.bed /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_2platforms_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives.vcf.gz |  ~/freebayes/freebayes_0.9.9.2/vcflib/vcfannotate -k filter -b /projects/scratch-data-backup/justin.zook/NA12878/bed/All.simplerepeatsnocov.annotate.bed | ~/freebayes/freebayes_0.9.9.2/vcflib/vcfannotate -k filter -b /projects/scratch-data-backup/justin.zook/NA12878/bed/superdupsmerged_all_sort.annotate.bed | ~/freebayes/freebayes_0.9.9.2/vcflib/vcfannotate -k filter -b /projects/scratch-data-backup/justin.zook/NA12878/bed/decoy.annotate.bed | ~/freebayes/freebayes_0.9.9.2/vcflib/vcfannotate -k filter -b /projects/scratch-data-backup/justin.zook/NA12878/bed/hg19.max5.regions.pad10.corr.annotate.bed | ~/freebayes/freebayes_0.9.9.2/vcflib/vcfannotate -k filter -b /projects/scratch-data-backup/justin.zook/NA12878/bed/dbvar_NA12878_CNVs_merged_all_sort.annotate.bed | ~/tabix-0.2.6/bgzip -fc >  /projects/scratch-data-backup/justin.zook/NA12878/Integration131103/AllFDAdatasets_131103_allcall_UGHapMerge_HetHomVarAll_2platforms_VQSRv2.18_2mindatasets_5minYesNoRatio_all.primitives.filterannotate.vcf.gz

