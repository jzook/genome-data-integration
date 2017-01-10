#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using an exome sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_folder="/home/pipeline/samples"
fastq_1="1.fastq.gz"
fastq_2="2.fastq.gz" #If using Illumina paired data
sample="sample_name"
group="read_group_name"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/pipeline/ref_hg19/ucsc.hg19.fasta"
dbsnp="/home/pipeline/ref_hg19/dbsnp_138.hg19.vcf.gz"
known_Mills_indels="/home/pipeline/ref_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
known_1000G_indels="/home/pipeline/ref_hg19/1000G_phase1.indels.hg19.sites.vcf.gz"
interval_file="/home/pipeline/ref_hg19/SeqCap_EZ_Exome_v2.hg19.interval"

# Update with the location of the Sentieon software package and license file
release_dir=/home/release/sentieon-genomics-201611
export SENTIEON_LICENSE=/home/Licenses/regression.lic

# Other settings
nt=16 #number of threads to use in computation
workdir="$PWD/test_interval" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir


# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
$release_dir/bin/bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$fastq_1 $fastq_folder/$fastq_2 | $release_dir/bin/sentieon util sort -r $fasta -o sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2. Metrics
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt --interval $interval_file -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$release_dir/bin/sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads
# ******************************************
$release_dir/bin/sentieon driver  -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
$release_dir/bin/sentieon driver  -t $nt --interval $interval_file -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam 

# ******************************************
# 4. Indel realigner
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels --interval_list $interval_file realigned.bam

# ******************************************
# 5. Base recalibration
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt --interval $interval_file -i realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt --interval $interval_file -i realigned.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
$release_dir/bin/sentieon plot bqsr -o recal_plots.pdf recal.csv

# ******************************************
# 6a. UG Variant caller
# ******************************************
$release_dir/bin/sentieon driver -r $fasta --interval $interval_file -t $nt -i realigned.bam -q recal_data.table --algo Genotyper -d $dbsnp --var_type BOTH --emit_conf=10 --call_conf=30 output-ug.vcf.gz

# ******************************************
# 6b. HC Variant caller
# ******************************************
$release_dir/bin/sentieon driver -r $fasta --interval $interval_file -t $nt -i realigned.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=30 output-hc.vcf.gz

# ******************************************
# 5b. ReadWriter to output recalibrated bam
# This stage is optional as variant callers
# can perform the recalibration on the fly
# using the before recalibration bam plus
# the recalibration table
# This stage should not include interval
# option to prevent read loss
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo ReadWriter recaled.bam

