#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs (~/in/beds/0/*, ~/in/beds/1/*, etc...)
#
dx-download-all-inputs --parallel


##Move all files from file arrays to the same data folder
mkdir data
mv ~/in/beds/*/* data

##Find names of processed vcf files and make sorted combined vcf using vcflib vcfoverlay
vars=""
for v in data/*.bed; do
  vars="$vars $v"
done
cat $vars | sed 's/^X/23/;s/^Y/24/;s/^MT/25/' | sort -k1,1n -k2,2n | sed 's/^23/X/;s/^24/Y/;s/^25/MT/' | grep -v ^# > out.bed



#
# Upload results
#
file_idbed=`dx upload out.bed -o "$prefix".bed --brief`
dx-jobutil-add-output "bedout" "$file_idbed"
