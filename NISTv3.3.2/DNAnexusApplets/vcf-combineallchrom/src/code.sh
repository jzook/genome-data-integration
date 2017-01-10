#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs (~/in/vcfs/0/*, ~/in/vcfs/1/*, etc...)
#
dx-download-all-inputs --parallel


##Move all files from file arrays to the same data folder
mkdir data
mv ~/in/vcfs/*/* data

#gunzip data/*.vcf.gz

##Find names of processed vcf files and make sorted combined vcf using vcflib vcfoverlay
for v in data/*.vcf.gz; do
  zgrep -v '0/0' "$v" > "$v".vcf 
done
vars=""
for v in data/*.vcf; do
  vars="$vars $v"
done
vcfoverlay $vars | sed 's/^X/23/;s/^Y/24/;s/^MT/25/' | sort -k1,1n -k2,2n | sed 's/^23/X/;s/^24/Y/;s/^25/MT/' | grep -v ^# > outa.vcf
grep ^# `ls data/*.vcf | head -1 | cut -f1` > head.vcf
cat head.vcf outa.vcf | bgzip -c > out.vcf.gz
tabix -p vcf out.vcf.gz



#
# Upload results
#
file_idvcf=`dx upload out.vcf.gz -o "$prefix".vcf.gz --brief`
dx-jobutil-add-output "vcfout" "$file_idvcf"
file_idvcftbi=`dx upload out.vcf.gz.tbi -o "$prefix".vcf.gz.tbi --brief`
dx-jobutil-add-output "vcfouttbi" "$file_idvcftbi"
