#!/bin/sh
# *******************************************
# Script to perform TNscope variant calling
# using a matched paired Tumor+normal sample with fastq
# files named normal_1.fastq.gz, normal_2.fastq.gz
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_folder="/home/pipeline/samples"
tumor_fastq_1="tumor_1.fastq.gz"
tumor_fastq_2="tumor_2.fastq.gz" #If using Illumina paired data
tumor_sample="tumor_sample_name"
tumor_group="tumor_read_group_name"
normal_fastq_1="normal_1.fastq.gz"
normal_fastq_2="normal_2.fastq.gz" #If using Illumina paired data
normal_sample="normal_sample_name"
normal_group="normal_read_group_name"
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
workdir="$PWD/test_TNscope" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1a. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
$release_dir/bin/bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$tumor_fastq_1 $fastq_folder/$tumor_fastq_2 | $release_dir/bin/sentieon util sort -o tumor_sorted.bam -t $nt --sam2bam -i -
# ******************************************
# 1b. Mapping reads with BWA-MEM, sorting for normal sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
$release_dir/bin/bwa mem -M -R "@RG\tID:$normal_group\tSM:$normal_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$normal_fastq_1 $fastq_folder/$normal_fastq_2 | $release_dir/bin/sentieon util sort -o normal_sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_sorted.bam --algo MeanQualityByCycle tumor_mq_metrics.txt --algo QualDistribution tumor_qd_metrics.txt --algo GCBias --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt
$release_dir/bin/sentieon plot metrics -o tumor_metrics-report.pdf gc=tumor_gc_metrics.txt qd=tumor_qd_metrics.txt mq=tumor_mq_metrics.txt isize=tumor_is_metrics.txt
# ******************************************
# 2b. Metrics for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_sorted.bam --algo MeanQualityByCycle normal_mq_metrics.txt --algo QualDistribution normal_qd_metrics.txt --algo GCBias --summary normal_gc_summary.txt normal_gc_metrics.txt --algo AlignmentStat --adapter_seq '' normal_aln_metrics.txt --algo InsertSizeMetricAlgo normal_is_metrics.txt
$release_dir/bin/sentieon plot metrics -o normal_metrics-report.pdf gc=normal_gc_metrics.txt qd=normal_qd_metrics.txt mq=normal_mq_metrics.txt isize=normal_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor sample
# ******************************************
$release_dir/bin/sentieon driver  -t $nt -i tumor_sorted.bam --algo LocusCollector --fun score_info tumor_score.txt
$release_dir/bin/sentieon driver  -t $nt -i tumor_sorted.bam --algo Dedup --rmdup --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam 
# ******************************************
# 3b. Remove Duplicate Reads for normal sample
# ******************************************
$release_dir/bin/sentieon driver  -t $nt -i normal_sorted.bam --algo LocusCollector --fun score_info normal_score.txt
$release_dir/bin/sentieon driver  -t $nt -i normal_sorted.bam --algo Dedup --rmdup --score_info normal_score.txt --metrics normal_dedup_metrics.txt normal_deduped.bam 

# ******************************************
# 4a. Indel realigner for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tumor_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tumor_realigned.bam
# ******************************************
# 4b. Indel realigner for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i normal_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels normal_realigned.bam

# ******************************************
# 5a. Base recalibration for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
$release_dir/bin/sentieon plot bqsr -o tumor_recal_plots.pdf tumor_recal.csv
# ******************************************
# 5b. Base recalibration for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels normal_recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_realigned.bam -q normal_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels normal_recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before normal_recal_data.table --after normal_recal_data.table.post normal_recal.csv
$release_dir/bin/sentieon plot bqsr -o normal_recal_plots.pdf normal_recal.csv

# ******************************************
# 6. Corealignment of tumor and normal
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tumor_realigned.bam -i normal_realigned.bam -q tumor_recal_data.table -q normal_recal_data.table --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tn_corealigned.bam

# ******************************************
# 7. Somatic and Structural variant calling
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tn_corealigned.bam --algo TNscope --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp output_tnscope.vcf.gz
