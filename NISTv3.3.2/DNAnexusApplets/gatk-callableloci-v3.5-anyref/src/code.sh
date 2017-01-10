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

#dx cat "$DX_PROJECT_CONTEXT_ID:/assets/hs37d5.fasta-index.tar.gz" | tar zxvf - # => genome.fa, genome.fa.fai, genome.dict
mv ~/in/ref/* .
tar zxvf *.tar.gz #unzip ref files

# Move BAI to same folder as BAM
mv ~/in/input_bai/* ~/in/input_bam/

#
# Set up options
#

# Calculate 80% of memory size, for java
mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`
java="java -Xmx${mem_in_mb}m"

# Move BAI to same folder as BAM
#mv ~/in/sorted_bai/* ~/in/sorted_bam/

#
# Run GATK
#
if [[ "$output_prefix" == "" ]]; then
  output_prefix="$sorted_bam_prefix"
fi
summary_path="$output_prefix".summary
bed_path="$output_prefix".bed
$java -jar /GenomeAnalysisTK.jar -T CallableLoci -R genome.fa -I "$input_bam_path" --summary "$summary_path" -o "$bed_path" $extra_options

#
# Upload results
#
mkdir -p ~/out/summary/ ~/out/bed_file/
mv "$summary_path" ~/out/summary/
mv "$bed_path" ~/out/bed_file/
dx-upload-all-outputs
