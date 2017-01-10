#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs 
#
dx download "$vcf_in" -o input.vcf
dx download "$callablelocibed" -o callableloci.bed
dx download "$targetsbed" -o targets.bed
#ls -l
#pwd

#
# Set up options
#

echo "$vcf_in"
echo "$chrom"

#
# Take bzipped vcfbeta from Complete Genomics and create filtered vcf and callable bed for integration
#
# To unzip vcf file
wc -l input.vcf

sed 's/chr//g' targets.bed | awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || $1 ~ /^#/) print }'  > targets_chrom.bed

sed 's/chr//g' input.vcf | awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || $1 ~ /^#/) print }' > input_chrom.vcf
wc -l input_chrom.vcf

awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || $1 ~ /^#/) print }' /human.b37.genome > /human.b37.chrom.genome

# To filter out variants with multi-allelic genotypes and create a bed file, since these are less accurate
grep '^#\|0/1\|1/1' input_chrom.vcf > output_Filter_chrom.vcf
grep -v '^#\|0/1\|1/1' input_chrom.vcf  | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' | slopBed -i stdin -b 50 -g /human.b37.chrom.genome | sort -k1,1 -k2,2n | mergeBed -i stdin > output_chrom_Filteredmultiallelic_slop50.bed


# To find overlap between callable loci and effective regions files
grep "CALLABLE" callableloci.bed | awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || $1 ~ /^#/) print }' > callableloci_callable.bed
intersectBed -a callableloci_callable.bed -b targets_chrom.bed | subtractBed -a stdin -b output_chrom_Filteredmultiallelic_slop50.bed | sort -k1,1 -k2,2n > output_chrom_callable_intersect.bed
complementBed -i output_chrom_callable_intersect.bed -g /human.b37.chrom.genome > output_chrom_callable_intersect_notcallable.bed


bgzip output_Filter_chrom.vcf
tabix -p vcf output_Filter_chrom.vcf.gz

#
# Upload results
#
name=`dx describe "$vcf_in" --name`
name="${name%.vcf}"
file_idvcf=`dx upload output_Filter_chrom.vcf.gz -o "$name"_"$chrom".vcf.gz --brief`
dx-jobutil-add-output "outvcfgz" "$file_idvcf"
file_idvcftbi=`dx upload output_Filter_chrom.vcf.gz.tbi -o "$name"_"$chrom".vcf.gz.tbi --brief`
dx-jobutil-add-output "outvcftbi" "$file_idvcftbi"
file_idnotcallablebed=`dx upload output_chrom_callable_intersect_notcallable.bed -o "$name"_notcallable_"$chrom".bed --brief`
dx-jobutil-add-output "outnotcallablebed" "$file_idnotcallablebed"
file_idcallablebed=`dx upload output_chrom_callable_intersect.bed -o "$name"_callable_"$chrom".bed --brief`
dx-jobutil-add-output "outcallablebed" "$file_idcallablebed"

