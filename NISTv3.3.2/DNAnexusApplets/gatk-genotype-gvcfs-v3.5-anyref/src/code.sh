#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs (~/in/vcfs/0/*, ~/in/vcfs/1/*, etc...)
#
dx-download-all-inputs --parallel

#
# Stream and unpack genome
#
#dx cat "$DX_PROJECT_CONTEXT_ID:/assets/hs37d5.fasta-index.tar.gz" | tar zxvf - # => genome.fa, genome.fa.fai, genome.dict
#dx cat "$DX_PROJECT_CONTEXT_ID:/assets/gatk_resources-2.8.b37.tar.gz" | tar zxvf - # => dbsnp_138.b37.vcf.gz, dbsnp_138.b37.vcf.gz.tbi, ...
mv ~/in/ref/* .
tar zxvf *.tar.gz #unzip ref files

#
# Set up options
#

# Calculate 80% of memory size, for java
mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`
java="java -Xmx${mem_in_mb}m"

# Move all files to the same folder
mkdir data
mv ~/in/vcfs/*/* data

# Compile command line
vars=""
for v in data/*vcf.gz; do
  vars="$vars -V $v"
done

#
# Run GATK
#
$java -jar /GenomeAnalysisTK.jar -nt `nproc`  -T GenotypeGVCFs -R genome.fa -o output.vcf $vars $extra_options

bgzip output.vcf
tabix -p vcf output.vcf.gz

#
# Upload results
#
mkdir -p ~/out/vcfgz/ ~/out/tbi/
mv output.vcf.gz ~/out/vcfgz/"$prefix".vcf.gz
mv output.vcf.gz.tbi ~/out/tbi/"$prefix".vcf.gz.tbi
dx-upload-all-outputs
