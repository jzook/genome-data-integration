#!/bin/sh
# *******************************************
# Script to perform TN seq variant calling
# using a Tumor sample with fastq files
# named tumor_1.fastq.gz, tumor_2.fastq.gz
# and data from a Panel of Normals and Cosmic DB
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_folder="/home/pipeline/samples"
tumor_fastq_1="tumor_1.fastq.gz"
tumor_fastq_2="tumor_2.fastq.gz" #If using Illumina paired data
tumor_sample="tumor_sample_name"
tumor_group="tumor_read_group_name"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/pipeline/ref_b37/Homo_sapiens_assembly19.fasta"
dbsnp="/home/pipeline/ref_b37/dbsnp_138.b37.vcf.gz"
known_Mills_indels="/home/pipeline/ref_b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
known_1000G_indels="/home/pipeline/ref_b37/1000G_phase1.indels.b37.vcf.gz"

# Update with the location of the panel of normal and CosmicDB vcf files
panel_of_normal_TNsnv="/home/pipeline/ref_b37/b37_panel_of_normal.vcf"
# We recommend that you create the panel of normal file with the corresponding algorithm that you plan to use for the somatic mutation calling. 
panel_of_normal_TNhaplotyper="/home/pipeline/ref_b37/b37_panel_of_normal.vcf"
cosmic_db="/home/pipeline/ref_b37/b37_cosmic_v54_120711.vcf.gz" 

# Update with the location of the Sentieon software package and license file
release_dir=/home/release/sentieon-genomics-201611
export SENTIEON_LICENSE=/home/Licenses/regression.lic

# Other settings
nt=16 #number of threads to use in computation
workdir="$PWD/test_TNseq_tumoronly" #Determine where the output files will be stored

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
$release_dir/bin/bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000  $fasta $fastq_folder/$tumor_fastq_1 $fastq_folder/$tumor_fastq_2 | $release_dir/bin/sentieon util sort -o tumor_sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_sorted.bam --algo MeanQualityByCycle tumor_mq_metrics.txt --algo QualDistribution tumor_qd_metrics.txt --algo GCBias --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt
$release_dir/bin/sentieon plot metrics -o tumor_metrics-report.pdf gc=tumor_gc_metrics.txt qd=tumor_qd_metrics.txt mq=tumor_mq_metrics.txt isize=tumor_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor sample
# ******************************************
$release_dir/bin/sentieon driver  -t $nt -i tumor_sorted.bam --algo LocusCollector --fun score_info tumor_score.txt
$release_dir/bin/sentieon driver  -t $nt -i tumor_sorted.bam --algo Dedup --rmdup --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam 

# ******************************************
# 4a. Indel realigner for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tumor_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tumor_realigned.bam

# ******************************************
# 5a. Base recalibration for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
$release_dir/bin/sentieon plot bqsr -o tumor_recal_plots.pdf tumor_recal.csv

# ******************************************
# 7. Somatic Variant Calling
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo TNsnv --tumor_sample $tumor_sample --pon $panel_of_normal_TNsnv --cosmic $cosmic_db --dbsnp $dbsnp --call_stats_out output-call.stats output-tnsnv.vcf.gz
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo TNhaplotyper --tumor_sample $tumor_sample --pon $panel_of_normal_TNhaplotyper --cosmic $cosmic_db --dbsnp $dbsnp output-tnhaplotyper.vcf.gz
