#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs 
#
dx download "$gvcf1" -o input1.vcf.gz
dx download "$gvcftbi1" -o input1.vcf.gz.tbi
dx download "$bed1" -o input1.bed
dx download "$gvcf2" -o input2.vcf.gz
dx download "$gvcftbi2" -o input2.vcf.gz.tbi
dx download "$bed2" -o input2.bed
dx download "$ref" -o ref.tar.gz
#ls -l
#pwd

tar zxvf *.tar.gz #unzip ref files
grep '^@SQ' genome.dict  | awk 'NR < 26' | sed 's/.*SN:\(.*\)\tLN:\(.*\)\tM5:.*/\1\t\2/' > human.genome
grep '^@SQ' genome.dict  | awk 'NR < 26' | sed 's/.*SN:\(.*\)\tLN:\(.*\)\tM5:.*/\1\t1\t\2/' > human.genome.bed

mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`
java="java -Xmx${mem_in_mb}m"

mv /usr/bin/GenomeAnalysisTK.jar .
mv /usr/bin/process_10X_gvcf.pl .

#
# Run GATK GenotypeGVCFs after renaming sample name in header
#
gunzip -c input1.vcf.gz | awk 'BEGIN {FS = OFS = "\t"} {if($1 ~ /^#CHROM/) $10="HP1"; print }' > input1_rename.vcf
gunzip -c input2.vcf.gz | awk 'BEGIN {FS = OFS = "\t"} {if($1 ~ /^#CHROM/) $10="HP2"; print }' > input2_rename.vcf

$java -jar GenomeAnalysisTK.jar -nt `nproc` -T GenotypeGVCFs -R genome.fa -o output.vcf -V input1_rename.vcf -V input2_rename.vcf -L $chrom -stand_emit_conf 2 -stand_call_conf 2

#Find diploid genotypes at sites where both haplotypes have homozygous genotypes
grep -v '0/1\|1/2\|,\*\|\./\.' output.vcf | sed 's/\(.*\tGT:\).*\(.\/\).*\/\(.\).*/\1GQ\t\2\3:99/' | cut -f1-10 > output_diploid.vcf

bgzip output_diploid.vcf
tabix output_diploid.vcf.gz

#Find regions from callableloci that are callable in both haplotypes
grep CALLABLE input1.bed > input1_CALLABLE.bed
grep CALLABLE input2.bed > input2_CALLABLE.bed

multiIntersectBed -i input1_CALLABLE.bed input2_CALLABLE.bed > intersect_CALLABLE.bed

awk '{if($4>1) sum+=$3-$2} END {print sum}' intersect_CALLABLE.bed
awk '{if ($4>1) print}' intersect_CALLABLE.bed > intersect_CALLABLEinboth.bed

#Find callable regions that are <1kb even after merging regions within 50bp, since these are questionable
mergeBed -i intersect_CALLABLEinboth.bed -d 50 | awk '{if ($3-$2<2000) print}' > intersect_CALLABLEinboth_merged50_lt1000.bed


#Find regions that have low quality or not homozygous from gvcf in both haplotypes
perl process_10X_gvcf.pl input1_rename.vcf $maxcov
perl process_10X_gvcf.pl input2_rename.vcf $maxcov

 awk '{sum+=$3;sum-=$2} END {print sum}' input1_rename_notcallable.bed 
 awk '{sum+=$3;sum-=$2} END {print sum}' input2_rename_notcallable.bed 

multiIntersectBed -i input1_rename_notcallable.bed input2_rename_notcallable.bed | mergeBed -i stdin > input_rename_notcallable_eitherhaplo.bed
 awk '{sum+=$3;sum-=$2} END {print sum}' input_rename_notcallable_eitherhaplo.bed

#Find final callable regions by subtracting gvcf low quality regions, heterozygous sites in either haplotype, and small callable regions from callableloci callable regions
 subtractBed -a intersect_CALLABLEinboth.bed  -b input_rename_notcallable_eitherhaplo.bed | subtractBed -a stdin  -b intersect_CALLABLEinboth_merged50_lt1000.bed | cut -f1-3 >  intersect_CALLABLEinboth_notingvcf10Xbed.bed
 awk '{sum+=$3;sum-=$2} END {print sum}' intersect_CALLABLEinboth_notingvcf10Xbed.bed
 

# Take complement to find not callable regions
awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || $1 ~ /^#/) print $1,0,$2 }' human.genome > human.chrom.genome.bed
wc -l human.chrom.genome.bed
subtractBed -a human.chrom.genome.bed -b intersect_CALLABLEinboth_notingvcf10Xbed.bed > intersect_CALLABLEinboth_notingvcf10Xbed_complement.bed

awk '{ sum+=$3; sum-=$2 } END { print sum }' intersect_CALLABLEinboth_notingvcf10Xbed_complement.bed


#
# Upload results
#
file_idvcf=`dx upload output_diploid.vcf.gz -o "$prefix".vcf.gz --brief`
dx-jobutil-add-output "vcfgz" "$file_idvcf"
file_idvcftbi=`dx upload output_diploid.vcf.gz.tbi -o "$prefix".vcf.gz.tbi --brief`
dx-jobutil-add-output "tbi" "$file_idvcftbi"
file_idnotcallablebed=`dx upload intersect_CALLABLEinboth_notingvcf10Xbed_complement.bed -o "$prefix"_notcallable.bed --brief`
dx-jobutil-add-output "outnotcallablebed" "$file_idnotcallablebed"
file_idcallablebed=`dx upload intersect_CALLABLEinboth_notingvcf10Xbed.bed -o "$prefix"_callable.bed --brief`
dx-jobutil-add-output "outcallablebed" "$file_idcallablebed"

