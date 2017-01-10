#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs
#
# ~/in/sorted_bam/*
# ~/in/sorted_bai/*
#
dx-download-all-inputs --parallel

mv ~/in/ref/* .
tar zxvf *.tar.gz #unzip ref files

#
# Set up options
#

# Calculate 80% of memory size, for java
mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`
java="java -Xmx${mem_in_mb}m"

# Move BAI to same folder as BAM
mv ~/in/sorted_bai/* ~/in/sorted_bam/

if [[ "$targets" != "" ]]; then
  extra_options="$extra_options -L ./in/targets/*"
fi

#
# Run GATK
#
$java -jar /GenomeAnalysisTK.jar -nct `nproc` -T HaplotypeCaller -R genome.fa -o output.vcf -I "$sorted_bam_path" -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -pairHMM VECTOR_LOGLESS_CACHING $extra_options

bgzip output.vcf
tabix -p vcf output.vcf.gz

#
# Upload results
#
if [[ "$output_prefix" == "" ]]; then
  output_prefix="$sorted_bam_prefix"
fi

mkdir -p ~/out/vcfgz/ ~/out/tbi/
mv output.vcf.gz ~/out/vcfgz/"$output_prefix".vcf.gz
mv output.vcf.gz.tbi ~/out/tbi/"$output_prefix".vcf.gz.tbi
dx-upload-all-outputs
