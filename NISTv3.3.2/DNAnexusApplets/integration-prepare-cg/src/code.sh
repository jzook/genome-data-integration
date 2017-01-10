#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs 
#
dx download "$vcf_in" -o input.vcf.bz2
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
bzip2 -d input.vcf.bz2
wc -l input.vcf

python /vcfBeta_to_VCF_simple.py -i input.vcf -o input_simplify.vcf
# To remove NO CALL lines and lines where genotype (GT) is '.' or '1' ('1' lines are SVs in breakend format) 
perl /cg_vcf_filter_v2.pl -vcf input_simplify.vcf -o input_Filter.vcf
wc -l input_Filter.vcf

# To get variants from chromosome specified
awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || $1 ~ /^#/) print }' input_Filter.vcf > output_Filter_chrom.vcf
wc -l output_Filter_chrom.vcf


#find not callable regions using CG's script (exclude partial calls, no-calls, and no-call regions from callable regions)
awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || $1 ~ /^#/) print }' /human.b37.genome > /human.b37.chrom.genome
python /vcf2bed.py input.vcf /example_of_no_ref_regions_input_file.bed output_notcallable.bed -print_no-calls
wc -l output_notcallable.bed
awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom) print $1 "\t" $2 "\t" $3}' output_notcallable.bed | slopBed -i stdin -b 50 -g /human.b37.genome | sort -k1,1 -k2,2n | mergeBed -i stdin > output_notcallable_slop50merged_chrom.bed
complementBed -i output_notcallable_slop50merged_chrom.bed -g /human.b37.chrom.genome > output_callable_slop50merged_chrom.bed


bgzip output_Filter_chrom.vcf
tabix -p vcf output_Filter_chrom.vcf.gz

#
# Upload results
#
name=`dx describe "$vcf_in" --name`
name="${name%.vcf.bz2}"
file_idvcf=`dx upload output_Filter_chrom.vcf.gz -o "$name"_"$chrom".vcf.gz --brief`
dx-jobutil-add-output "outvcfgz" "$file_idvcf"
file_idvcftbi=`dx upload output_Filter_chrom.vcf.gz.tbi -o "$name"_"$chrom".vcf.gz.tbi --brief`
dx-jobutil-add-output "outvcftbi" "$file_idvcftbi"
file_idnotcallablebed=`dx upload output_notcallable_slop50merged_chrom.bed -o "$name"_notcallable_"$chrom".bed --brief`
dx-jobutil-add-output "outnotcallablebed" "$file_idnotcallablebed"
file_idcallablebed=`dx upload output_callable_slop50merged_chrom.bed -o "$name"_callable_"$chrom".bed --brief`
dx-jobutil-add-output "outcallablebed" "$file_idcallablebed"

