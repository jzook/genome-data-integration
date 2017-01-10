#!/bin/sh
# *******************************************
# Script to perform joint calling on 3 samples
# with fastq files named sample<i>_1.fastq.gz
# and sample<i>_2.fastq.gz
# *******************************************

# Update with the fullpath location and suffix of the fastq files
fastq_folder="/home/pipeline/samples"
fastq_1_suffix="1.fastq.gz"
fastq_2_suffix="2.fastq.gz" #If using Illumina paired data
platform="ILLUMINA" #platform

# Update with the location of the reference data files
regions="/home/pipeline/ref_hg19/TruSeq_exome_targeted_regions.hg19.bed.chr"
fasta="/home/pipeline/ref_hg19/ucsc.hg19.fasta"
dbsnp="/home/pipeline/ref_hg19/dbsnp_138.hg19.vcf.gz"
known_Mills_indels="/home/pipeline/ref_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
known_1000G_indels="/home/pipeline/ref_hg19/1000G_phase1.indels.hg19.sites.vcf.gz"

# Update with the location of the Sentieon software package and license file
release_dir=/home/release/sentieon-genomics-201611
export SENTIEON_LICENSE=/home/Licenses/regression.lic

# Other settings
nt=16 #number of threads to use in computation
workdir="$PWD/test_joint" #Determine where the output files will be stored


# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1

# ******************************************
# 0. Process all samples independently
# ******************************************
#update with the prefix of the fastq files
for sample in sample1 sample2 sample3; do
  group=$sample
  mkdir $workdir/$sample
  cd $workdir/$sample

  # ******************************************
  # 1. Mapping reads with BWA-MEM, sorting
  # ******************************************
  $release_dir/bin/bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/${sample}_$fastq_1_suffix $fastq_folder/${sample}_$fastq_2_suffix | $release_dir/bin/sentieon util sort -r $fasta -o sorted.bam -t $nt --sam2bam -i -
  
  # ******************************************
  # 2. Metrics
  # ******************************************
  $release_dir/bin/sentieon driver -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
  $release_dir/bin/sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt
  
  # ******************************************
  # 3. Remove Duplicate Reads
  # ******************************************
  $release_dir/bin/sentieon driver  -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
  $release_dir/bin/sentieon driver  -t $nt -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam 
  
  # ******************************************
  # 4. Indel realigner
  # ******************************************
  $release_dir/bin/sentieon driver -r $fasta  -t $nt -i deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels realigned.bam
  
  # ******************************************
  # 5. Base recalibration
  # ******************************************
  $release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table
  $release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
  $release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
  $release_dir/bin/sentieon plot bqsr -o recal_plots.pdf recal.csv
  
  # ******************************************
  # 6b. HC Variant caller for GVCF
  # ******************************************
  $release_dir/bin/sentieon driver -r $fasta --interval $regions -t $nt -i realigned.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=30 --emit_mode gvcf output-hc.g.vcf.gz
  gvcf_inputs="$gvcf_inputs -v $workdir/$sample/output-hc.g.vcf.gz"
done

# ******************************************
# Perform the joint calling 
# ******************************************
$release_dir/bin/sentieon driver -r $fasta --algo GVCFtyper $gvcf_inputs -d $dbsnp $workdir/output-jc.vcf.gz
