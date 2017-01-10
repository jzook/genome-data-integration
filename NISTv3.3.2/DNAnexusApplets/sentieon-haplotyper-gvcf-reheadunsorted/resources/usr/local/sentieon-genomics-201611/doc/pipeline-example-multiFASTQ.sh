#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with more than one
# set of input fastq files (4 in this example)
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_folder="/home/pipeline/samples"
num_sets="4"
# fastq should be named $fastq_prefix$set_number$fastq_suffix_1
# for this case set1_1.fastq.gz, set1_2.fastq.gz, set2_1.fastq.gz, set2_2.fastq.gz, set3_1.fastq.gz, set3_2.fastq.gz
fastq_prefix="set"
fastq_suffix_1="_1.fastq.gz"
fastq_suffix_2="_2.fastq.gz"
group_prefix="read_group_name"
sample="sample_name"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/pipeline/ref_hg19/ucsc.hg19.fasta"
dbsnp="/home/pipeline/ref_hg19/dbsnp_138.hg19.vcf.gz"
known_Mills_indels="/home/pipeline/ref_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
known_1000G_indels="/home/pipeline/ref_hg19/1000G_phase1.indels.hg19.sites.vcf.gz"

# Update with the location of the Sentieon software package and license file
release_dir=/home/release/sentieon-genomics-201611
export SENTIEON_LICENSE=/home/Licenses/regression.lic

# Other settings
nt=16 #number of threads to use in computation
workdir="$PWD/test_miltiFASTQ" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir


# ******************************************
# 1. Mapping each set of input fastq with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
for i in $(seq 1 $num_sets); do
 bam_input="$bam_input -i sorted_set$i.bam"
 $release_dir/bin/bwa mem -M -R "@RG\tID:${group_prefix}_$i\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$fastq_prefix$i$fastq_suffix_1 $fastq_folder/$fastq_prefix$i$fastq_suffix_2 | $release_dir/bin/sentieon util sort -r $fasta -o sorted_set$i.bam -t $nt --sam2bam -i -
done

# ******************************************
# 2. Metrics on the multiple sorted BAM files
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt $bam_input --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$release_dir/bin/sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads on the multiple sorted BAM files
# ******************************************
$release_dir/bin/sentieon driver  -t $nt $bam_input --algo LocusCollector --fun score_info score.txt
$release_dir/bin/sentieon driver  -t $nt $bam_input --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam 

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
# 6a. UG Variant caller
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Genotyper -d $dbsnp --var_type BOTH --emit_conf=10 --call_conf=30 output-ug.vcf.gz

# ******************************************
# 6b. HC Variant caller
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=30 output-hc.vcf.gz

