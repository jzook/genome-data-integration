#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs 
#
#dx-download-all-inputs --parallel
dx download "$vcfhighconfgz" -o vcfhighconfgz.vcf.gz
dx download "$vcfphasedgz" -o vcfphasedgz.vcf.gz
dx download "$vcfhighconftbi" -o vcfhighconfgz.vcf.gz.tbi
dx download "$vcfphasedtbi" -o vcfphasedgz.vcf.gz.tbi
dx download "$rtgsdf" -o rtgsdf.tar.gz
#
# Stream and unpack genome
#
# => RTG reference sdf files
  tar zxvf rtgsdf.tar.gz





#
# Transfer phasing from RTG Segregation Phasing and Illumina Platinum Genomes
# to another call set, preserving the original calls.
#

RTG=rtg

# Path to formatted reference genome (use rtg format to make this)
REF=rtgsdf

# Path to calls to transfer phasing onto
CALLS=vcfhighconfgz.vcf.gz

# Paths to source phasing set
SP=vcfphasedgz.vcf.gz

# Extra vcfeval options needed for phase transfer (for rtg 3.7).
# Note these --X options may change in future releases.
PHASE_TRANSFER_OPTS="--ref-overlap --Xobey-phase=true,true --XXcom.rtg.vcf.eval.custom-path-processor=phase-transfer"

# Header field declaration fors attributes used to preserve original phasing
IGT_HEADER="##FORMAT=<ID=IGT,Number=1,Type=String,Description=\\\"Original input genotype\\\">"
IPS_HEADER="##FORMAT=<ID=IPS,Number=1,Type=String,Description=\\\"Phase set for IGT\\\">"

# Report the amount of phasing in a VCF
function report_phasing() {
    local message=$1
    local vcf=$2
    local phased=$(${RTG} vcfstats "${vcf}" | grep '^Phased' | awk '{print $(NF-1),$NF}')
    echo "${message}: ${phased:-none}"
}

# Remove any unphased heterozygous calls for a particular sample from the
# input.  This is just a safety net to make sure unphased calls are not
# used for transfer.  It drops 11 calls for the SP and none for PG.
# Haploid calls are also omitted since they are irrelevant for phasing.
function remove_unphased_hets() {
    local dirty=$1
    local clean=$2
    local sample=$3
    echo "Filtering out unphased hets from ${dirty}"
    ${RTG} vcffilter -i ${dirty} -o ${clean} --keep-expr \
           "split=${sample}.GT.split(\"/\"); \
            split.length<2 || split[0].equals(split[1]);" \
        || exit 1
}

# Remove any existing phasing for the specified sample, copying the
# existing phasing and phase set information into new field IGT and IPS.
function remove_phasing() {
    local calls=$1
    local outfile=$2
    local sample=$3
    echo "Archiving existing phasing from ${calls}"
    ${RTG} vcffilter -i ${calls} -o - --javascript \
           "ensureFormatHeader(\"${IGT_HEADER}\"); \
            ensureFormatHeader(\"${IPS_HEADER}\"); \
            function record() { \
              ${sample}.IGT=${sample}.GT; \
              ${sample}.IPS=${sample}.PS; \
              ${sample}.GT=${sample}.GT.replace(\"|\",\"/\"); \
            }" \
        | ${RTG} vcfsubset --remove-format PS -i - -o ${outfile} || exit 1
}

# Remove any existing local phasing for the specified sample, copying the
# existing phasing and phase set information into new field IGT and IPS.
function remove_local_phasing() {
    local calls=$1
    local outfile=$2
    local sample=$3
    echo "Archiving existing phasing from ${calls}"
    ${RTG} vcffilter -i ${calls} -o - --javascript \
           "ensureFormatHeader(\"${IGT_HEADER}\"); \
            ensureFormatHeader(\"${IPS_HEADER}\"); \
            function record() { \
              if (${sample}.PS != \"PATMAT\") { \
                ${sample}.IGT=${sample}.GT; \
                ${sample}.IPS=${sample}.PS; \
                ${sample}.GT=${sample}.GT.replace(\"|\",\"/\"); \
              } \
            }" \
        | ${RTG} vcfsubset --remove-format PS -i - -o ${outfile} || exit 1
}


# Perform the actual phase transfer. It expects to be given the source
# of phasing, the call set to be phase, where to write the output VCF,
# and the sample information (which can be a single sample or in the
# form phase-src-sample,phase-onto-sample).
function phase_transfer() {
    local phase_src=$1
    local phase_onto=$2
    local outfile=$3
    local sample_spec=$4
    echo "Using vcfeval to phase transfer from ${phase_src} onto ${phase_onto}"
    phasesp=phase-transfer-sp
    ${RTG} vcfeval ${PHASE_TRANSFER_OPTS} -t ${REF} -b ${phase_src} -c ${phase_onto} -o ${outfile}  --sample ${sample_spec} || exit 1
}




# Remove any unphased heterozygous calls to ensure full phasing of sources
remove_unphased_hets ${SP} clean-sp.vcf.gz ${PhasedSample}

# Prepare call set by copying any existing phasing into IGT, IPS; correct problem with missing IPS field using sed
report_phasing "Original calls" ${CALLS}
if [ "$removephasing" = "yes" ]; then
	remove_phasing ${CALLS} unphased-calls-IPSerr.vcf.gz ${HighConfSample}
#	gunzip -c unphased-calls-IPSerr.vcf.gz | sed 's/\(.\/.$\)/\1:\./' | ${RTG} bgzip - > unphased-calls.vcf.gz && ${RTG} index -f vcf unphased-calls.vcf.gz
	gunzip -c unphased-calls-IPSerr.vcf.gz | awk 'BEGIN {FS = OFS = "\t"} {if(!($10 ~ /:.*:.*:.*:.*:.*:/) && !($1 ~ /^#/)) $10=$10":.";  print }' | ${RTG} bgzip - > unphased-calls.vcf.gz && ${RTG} index -f vcf unphased-calls.vcf.gz
	zgrep -nm 4 IPS  unphased-calls.vcf.gz 
fi
if [ ! "$removephasing" = "yes" ]; then
	remove_local_phasing ${CALLS} unphased-calls.vcf.gz ${HighConfSample}
fi

# Transfer phasing from Segregation Phasing
phase_transfer clean-sp.vcf.gz unphased-calls.vcf.gz phase-transfer-sp ${PhasedSample},${HighConfSample}
report_phasing "Phase Transfered" phase-transfer-sp/calls.vcf.gz


# Strip out unwanted attributes produced during phase transfer
${RTG} vcfsubset -i phase-transfer-sp/calls.vcf.gz -o - --remove-format OGT --remove-info SYNC,CALL,CALL_WEIGHT,PHASE | grep -v '^##CL\|^##RUN-ID' | ${RTG} bgzip - >filtered-phase-transfer.vcf.gz && ${RTG} index -f vcf filtered-phase-transfer.vcf.gz

echo "Final pedigree transfer result in filtered-phase-transfer.vcf.gz"

${RTG} vcffilter -i filtered-phase-transfer.vcf.gz -o - --javascript "ensureFormatHeader(\"##FORMAT=<ID=PS,Number=1,Type=String,Description=\\\"Phase set for GT\\\">\"); function record() {if(INTEGRATION.GT==\"1/1\") { INTEGRATION.IPS=\".\"; INTEGRATION.PS=\"HOMVAR\"; INTEGRATION.GT=\"1|1\";} else {if((INTEGRATION.GT==\"0/1\" || INTEGRATION.GT==\"1/2\" || INTEGRATION.GT==\"2/1\" || INTEGRATION.GT==\"1/0\") ) {if(INTEGRATION.IPS.length>1) {INTEGRATION.PS=INTEGRATION.IPS; INTEGRATION.GT=INTEGRATION.IGT;} else {INTEGRATION.PS=\".\";};} else { if((INTEGRATION.IPS.length<2)) { INTEGRATION.IPS=\".\";} INTEGRATION.PS=\"PATMAT\";};};}" | awk 'BEGIN {FS = OFS = "\t"} {if(!($10 ~ /:.*:.*:.*:.*:.*:/) && !($1 ~ /^#/)) $10=$10":.:.";  print }' | ${RTG} bgzip - > filtered-phase-transfer-allphasing.vcf.gz && ${RTG} index -f vcf filtered-phase-transfer-allphasing.vcf.gz
report_phasing "Phase Transferred plus original phasing" filtered-phase-transfer-allphasing.vcf.gz

${RTG} vcffilter -i filtered-phase-transfer-allphasing.vcf.gz -o - --keep-expr 'INTEGRATION.GT.indexOf("|")>0&&!has(INTEGRATION.PS)' --no-header | wc -l || echo "ignore stderr here"
${RTG} vcffilter -i filtered-phase-transfer-allphasing.vcf.gz -o - --keep-expr 'INTEGRATION.GT.indexOf("|")==-1&&has(INTEGRATION.PS)' --no-header | wc -l || echo "ignore stderr here"
${RTG} vcffilter -i filtered-phase-transfer-allphasing.vcf.gz -o - --keep-expr 'has(INTEGRATION.PS)&&INTEGRATION.PS!="PATMAT"&&INTEGRATION.PS!="HOMVAR"' | ${RTG} vcffilter -i - --javascript 'function record() {print(CHROM+":"+INTEGRATION.PS)}' | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c || echo "ignore stderr here"


#
# Upload results
#

file_idoutputorigannotvcf=`dx upload filtered-phase-transfer-allphasing.vcf.gz -o "$prefix".vcf.gz --brief`
dx-jobutil-add-output "outputorigannotvcf" "$file_idoutputorigannotvcf"

file_idoutputorigannotvcftbi=`dx upload filtered-phase-transfer-allphasing.vcf.gz.tbi -o "$prefix".vcf.gz.tbi --brief`
dx-jobutil-add-output "outputorigannotvcftbi" "$file_idoutputorigannotvcftbi"

