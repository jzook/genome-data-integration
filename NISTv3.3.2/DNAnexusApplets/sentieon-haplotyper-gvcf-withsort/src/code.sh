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

#
# Stream and unpack genome and other GATK files
#
mkdir genome
gunzip -c $genome_fasta_path > genome/$(basename $genome_fasta_path .gz)  # => genome/<ref>
mv $genome_fastaindex_path genome/  # => genome/<ref>.fai
genome_file=`ls genome/*.fai`     # Locate a file called <ref>.fai
genome_file="${genome_file%.fai}" # Remove the fai suffix to keep the <ref>

mkdir resources
if [ ! "$gatk_resources_path" = "" ]; then
 tar -zxvf $gatk_resources_path -C resources --strip 1
fi

#
# Set up options
#
release_dir=/usr/local/sentieon-genomics-201606.rc1 #TODO: make version configurable
release_dir=$(echo /usr/local/sentieon-genomics-*|awk '{print $1}')
sentieon_procs=$(nproc)
 #set the DNAnexus license env
 dx download "$license_server_file" -o license_server.info
 license_server=$(grep license_server_location license_server.info|awk -F"license_server_location=" '{print $2}')
 auth_token=$(grep auth_token license_server.info|awk -F"auth_token=" '{print $2}')
 export SENTIEON_AUTH_MECH=DX
 export SENTIEON_AUTH_DATA="{\"auth_token_type\": \"Bearer\", \"auth_token\": \"$auth_token\"}"
 export SENTIEON_LICENSE="$license_server"

# Move BAI to same folder as BAM
mv ~/in/sorted_bai/* ~/in/sorted_bam/

#additional inputs
dbsnp=""
if [ ! "$gatk_resources_path" = "" ]; then
 for vcfgz in resources/dbsnp*.vcf.gz; do
  dbsnp="-d $vcfgz" #Sentieon only accepts 1 dbsnp, used to be #dbsnp="$dbsnp -d $vcfgz"
 done
fi

if [ ! "$targets" = "" ]; then
  for bedfile in ./in/targets/*; do
    extra_driver_options="$extra_driver_options --interval $bedfile"
  done
fi

#
# Run Sentieon Haplotyper
#
$release_dir/bin/sentieon util sort -o sorted.bam -t $sentieon_procs -i "$sorted_bam_path"
$release_dir/bin/sentieon driver $extra_driver_options -r $genome_file  -t $sentieon_procs -i sorted.bam --algo Haplotyper --emit_mode gvcf $extra_algo_options $dbsnp output.gvcf.gz

#
# Upload results
#
mkdir -p ~/out/gvcfgz/ ~/out/tbi/
mv output.gvcf.gz ~/out/gvcfgz/"$sorted_bam_prefix".gvcf.gz
mv output.gvcf.gz.tbi ~/out/tbi/"$sorted_bam_prefix".gvcf.gz.tbi

dx-upload-all-outputs # uploads outputs with the pattern ~/out/[OUTPUT_NAME]/[FILE_NAME]
