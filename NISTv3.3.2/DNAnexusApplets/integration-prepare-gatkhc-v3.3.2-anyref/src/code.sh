#!/bin/bash

#v3.3.1 - don't include gvcf lines covering single base and only <NON_REF> in the ALT column since these lines often have low GQ after an indel with high GQ that looks unambiguous

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs 
#
dx download "$gvcf" -o input.vcf.gz
dx download "$gvcftbi" -o input.vcf.gz.tbi
dx download "$ref" -o ref.tar.gz
#ls -l
#pwd
mv /usr/bin/*.pl .

tar zxvf *.tar.gz #unzip ref files
grep '^@SQ' genome.dict  | awk 'NR < 26' | sed 's/.*SN:\(.*\)\tLN:\(.*\)\tM5:.*/\1\t\2/' > human.genome
grep '^@SQ' genome.dict  | awk 'NR < 26' | sed 's/.*SN:\(.*\)\tLN:\(.*\)\tM5:.*/\1\t1\t\2/' > human.genome.bed
#ls -l
#pwd


#

gunzip input.vcf.gz

perl process_Illumina_gvcf.pl input.vcf $maxcov

 awk '{sum+=$3;sum-=$2} END {print sum}' input_notcallable.bed 


awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || $1 ~ /^#/) print $1,0,$2 }' human.genome > human.chrom.genome.bed

wc -l human.chrom.genome.bed
# Take complement to find callable regions
subtractBed -a human.chrom.genome.bed -b input_notcallable.bed > input_callable.bed

awk '{ sum+=$3; sum-=$2 } END { print sum }' input_callable.bed


#
# Upload results
#
name=`dx describe "$gvcf" --name`
name="${name%.vcf.gz}"
file_idnotcallablebed=`dx upload input_notcallable.bed -o "$name"_notcallable.bed --brief`
dx-jobutil-add-output "outnotcallablebed" "$file_idnotcallablebed"
file_idcallablebed=`dx upload input_callable.bed -o "$name"_callable.bed --brief`
dx-jobutil-add-output "outcallablebed" "$file_idcallablebed"

