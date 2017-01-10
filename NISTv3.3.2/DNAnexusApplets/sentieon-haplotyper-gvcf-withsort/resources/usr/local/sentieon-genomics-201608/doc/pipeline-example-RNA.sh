#!/bin/sh
# *******************************************
# Script to perform RNA variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_folder="/home/pipeline/samples"
fastq_1="rna1.fastq.gz"
fastq_2="rna2.fastq.gz"
sample="sample_name"
group="read_group_name"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/pipeline/ref_hg19/ucsc.hg19.fasta"
dbsnp="/home/pipeline/ref_hg19/dbsnp_138.hg19.vcf.gz"
known_Mills_indels="/home/pipeline/ref_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
known_1000G_indels="/home/pipeline/ref_hg19/1000G_phase1.indels.hg19.sites.vcf.gz"
#uncomment if START should generate a new genomeDir
star_fasta="/home/pipeline/ref_hg19/ucsc.hg19.fasta.genomeDir"

# Update with the location of the Sentieon software package and license file
release_dir=/home/release/sentieon-genomics-201606
export SENTIEON_LICENSE=/home/Licenses/regression.lic
star_binary="/home/release/other_tools/STAR/bin/Linux_x86_64_static/STAR"

# Other settings
nt=16 #number of threads to use in computation
workdir="$HOME/testcaseRNA" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1. Mapping reads with STAR
# ******************************************
if [ -z "$star_fasta" ]; then
  star_fasta="genomeDir"
  # The genomeDir generation could be reused
  mkdir $star_fasta
  $star_binary --runMode genomeGenerate --genomeDir $star_fasta --genomeFastaFiles $fasta --runThreadN $nt
fi
#perform the actual alignment and sorting
$star_binary --twopassMode Basic --genomeDir $star_fasta --runThreadN $nt --outSAMtype BAM SortedByCoordinate --twopass1readsN -1 --sjdbOverhang 75 --readFilesIn $fastq_folder/$fastq_1 $fastq_folder/$fastq_2 --readFilesCommand zcat --outSAMattrRGline ID:$group SM:$sample PL:$platform
mv Aligned.sortedByCoord.out.bam sorted.bam
$release_dir/bin/sentieon util index sorted.bam

# ******************************************
# 2. Metrics
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$release_dir/bin/sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads
# ******************************************
$release_dir/bin/sentieon driver  -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
$release_dir/bin/sentieon driver  -t $nt -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam 

# ******************************************
# 4. Split reads at Junction
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo RNASplitReadsAtJunction --reassign_mapq 255:60 split.bam

# ******************************************
# 5. Indel realigner
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i split.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels realigned.bam

# ******************************************
# 6. Base recalibration
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
$release_dir/bin/sentieon plot bqsr -o recal_plots.pdf recal.csv

# ******************************************
# 7. HC Variant caller for RNA
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Haplotyper -d $dbsnp --trim_soft_clip --emit_conf=20 --call_conf=20 output-hc-rna.vcf.gz

