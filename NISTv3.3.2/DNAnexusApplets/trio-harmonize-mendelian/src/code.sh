#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs 
#
#dx-download-all-inputs --parallel
dx download "$vcfchildgz" -o vcfchildgz.vcf.gz
dx download "$vcfdadgz" -o vcfdadgz.vcf.gz
dx download "$vcfmomgz" -o vcfmomgz.vcf.gz
dx download "$vcfchildtbi" -o vcfchildgz.vcf.gz.tbi
dx download "$vcfdadtbi" -o vcfdadgz.vcf.gz.tbi
dx download "$vcfmomtbi" -o vcfmomgz.vcf.gz.tbi
dx download "$childbed" -o child.bed
dx download "$dadbed" -o dad.bed
dx download "$mombed" -o mom.bed
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

#
# Variant harmonization.
#


MERGED=merged.vcf.gz

PRE_A=vcfchildgz.vcf.gz
PRE_B=vcfdadgz.vcf.gz
PRE_C=vcfmomgz.vcf.gz

A=child.vcf.gz
B=dad.vcf.gz
C=mom.vcf.gz

# Base common command line options to vcfeval for accumulation and recoding.
COMMON_OPTS="-t ${REF} --XXcom.rtg.vcf.eval.maximize=calls-min-base --XXcom.rtg.vcf.eval.max-paths 500000 --XXcom.rtg.vcf.eval.custom-variant-factory=all,sample --all-records --ref-overlap"

# Options to enable haploid allele accumulation mode
ACCUM_OPTS="${COMMON_OPTS} --XXcom.rtg.vcf.eval.custom-path-processor=alleles --squash-ploidy"

# Options for the recoding step
RECODE_OPTS="${COMMON_OPTS} --XXcom.rtg.vcf.eval.custom-path-processor=recode"



function relabel {
    gunzip -c "$1" | sed "/^#CHROM/s/${3}/${4}/"| ${RTG} bgzip - >"$2" && ${RTG} index -f vcf "$2" || exit 1
}

# Run accumulation for all SAMPLES specified.  This uses haploid
# accumulation mode, which can capture at most one ALT per site. Hence
# we do two vceval runs per sample in order to handle any cases where
# the sample had a GT such as 1/2.
function accumulatehap {
    local sample=$1
    local vcf=${sample}.vcf.gz
    shift
    echo "Accumulating alleles from ${sample}"
    # First sample outside the loop as we're starting with empty baseline
    [ -r ac-${sample}-1 ] || ${RTG} vcfeval ${ACCUM_OPTS} -o ac-${sample}-1 -b ${EMPTY_VCF} -c ${vcf} --sample ${sample} || echo "error here can be ignored"
    [ -r ac-${sample}-2 ] || ${RTG} vcfeval ${ACCUM_OPTS} -o ac-${sample}-2 -b ac-${sample}-1/alleles.vcf.gz -c ac-${sample}-1/alternate.vcf.gz --sample ${sample} || exit 1
    [ -r ac-${sample}-2/alt-stats ] || gunzip -c ac-${sample}-2/alleles.vcf.gz | grep -v '^#' | cut -d$'\t' -f 5 | tr , ' ' | wc -lw >ac-${sample}-2/alt-stats
    while [ "$1" ]; do
        local prev=${sample}
        sample=$1
        vcf=${sample}.vcf.gz
        echo "Accumulating alleles from ${sample}"
        shift
        if [ ! -r ac-${sample}-2 ]; then
            [ -r ac-${sample}-1 ] || ${RTG} vcfeval ${ACCUM_OPTS} -o ac-${sample}-1 -b ac-${prev}-2/alleles.vcf.gz -c ${vcf} --sample ${sample} || exit 1
            [ -r ac-${sample}-2 ] || ${RTG} vcfeval ${ACCUM_OPTS} -o ac-${sample}-2 -b ac-${sample}-1/alleles.vcf.gz -c ac-${sample}-1/alternate.vcf.gz --sample ${sample} || exit 1
            [ -r ac-${sample}-2/alt-stats ] || gunzip -c ac-${sample}-2/alleles.vcf.gz | grep -v '^#' | cut -d$'\t' -f 5 | tr , ' ' | wc -lw >ac-${sample}-2/alt-stats
        fi
    done

    # Remove worst culprits for complex regions from the population set.
    gunzip -c ${sample}.vcf.gz | sed '/^[^#]/q' | ${RTG} bgzip - >simple-sample.vcf.gz || echo "ignore error here since sed terminates process"
    ${RTG} index -f vcf simple-sample.vcf.gz
    ${RTG} vcfeval ${RECODE_OPTS} -o hard-regions -b ac-${sample}-2/alleles.vcf.gz -c simple-sample.vcf.gz
    grep 'too complex' <hard-regions/vcfeval.log | sed "s/.* \([^:]*\):/\1:/;s/\..*//;y/:-/  /" | sort -k 2,2n | awk 'BEGIN {OFS = "\t"} {print $1,$2-1,$3-1}' >hard-regions.bed || echo "no difficult regions"

    # Write final accumulated population allele VCF.
    ${RTG} vcffilter --exclude-bed hard-regions.bed -i ac-${sample}-2/alleles.vcf.gz -o population-alleles.vcf.gz || cp ac-${sample}-2/alleles.vcf.gz population-alleles.vcf.gz
	${RTG} index -f vcf population-alleles.vcf.gz
}



# Relabel samples names so that different inputs have different sample names
relabel ${PRE_A} ${A} INTEGRATION child
relabel ${PRE_B} ${B} INTEGRATION dad
relabel ${PRE_C} ${C} INTEGRATION mom

# Create a simple file listing the sample names
SAMPLES=samples.txt
cat >${SAMPLES} <<EOF
child
dad
mom
EOF

# Create an empty VCF with appropriate header as starting point for accumulation
EMPTY_VCF=empty.vcf.gz
#gunzip -c ${A} | sed 's/GS.*/dummy/;/^#CHROM/q' | ${RTG} bgzip - >empty.vcf.gz && ${RTG} index -f vcf ${EMPTY_VCF} || exit 1
zgrep '^#' ${A} | ${RTG} bgzip - > empty.vcf.gz && ${RTG} index -f vcf ${EMPTY_VCF} || exit 1

accumulatehap $(cat ${SAMPLES})

POPULATION_VCF=population-alleles.vcf.gz  # Created by accumulate.sh
# Use the same population allele baseline for all samples
RECODE_OPTS="${RECODE_OPTS} -b ${POPULATION_VCF}"

if [ ! -f ${POPULATION_VCF} ]; then
    echo "Population allele VCF does not exist: ${POPULATION_VCF}"
    exit 1
fi
if [ ! -f ${SAMPLES} ]; then
    echo "Sample list file does not exist: ${SAMPLES}"
    exit 1
fi
set -- $(cat ${SAMPLES})
echo "Sample list file contains $# entries"


# Recode all the samples against accumulated baseline
function recode {
    while [ "$1" ]; do
        local sample=$1
        local vcf=${sample}.vcf.gz
        echo "Recoding sample ${sample}"
        [ -r recode-${sample} ] || ${RTG} vcfeval ${RECODE_OPTS} -o recode-${sample} -c ${vcf} --sample ${sample} || exit 1
        shift
    done
}

recode $(cat ${SAMPLES})


#
# Measure Mendelianness in high-confidence regions.

BEDTOOLS=bedtools

HC=child.bed
HCM=dad.bed
HCF=mom.bed

PED=pedigree.ped
cat >${PED} <<EOF
#fam-id ind-id  pat-id  mat-id  sex     status
1	dad	0	0	1	0
1	mom	0	0	2	0
1	child	dad	mom	1	0
EOF

IHC=intersect-hc.bed
intersectBed -a ${HC} -b ${HCM} | intersectBed -a - -b ${HCF} > ${IHC}

# Piece of javascript that will count records without needing to write them
JS="total=0; function record(){total++} function end(){print(total)}"

# Compute and produce a summary line on Mendelian violations for a supplied
# set of bed regions.
function report {
    local message=$1
    local violations_vcf=$2
    local global_vcf=$3
    local regions=$4
    local i=$(${RTG} vcffilter -i ${violations_vcf} --javascript "${JS}" --bed-regions ${regions})
    local t=$(${RTG} vcffilter -i ${global_vcf} --javascript "${JS}" --bed-regions ${regions})
    local p=$(echo "scale=2;100.0*$i/$t" | bc | sed 's/^\./0./')
    echo "${message}: ${i}/${t} (${p}%)"
}

# Run the actual Mendelianness checking
function mendelian {
    local message=$1
    local src_vcf=$2
    local violations_out_vcf=$3
    local phased_out_vcf=$4
    echo "${message}: "$(${RTG} mendelian --pedigree ${PED} -t ${REF} --all-records --lenient --output-inconsistent ${violations_out_vcf} --output-consistent ${phased_out_vcf} --Xphase -i ${src_vcf} )
}

# Measure Mendelianness before harmonization
[[ -r merged.vcf.gz ]] || ${RTG} vcfmerge -o merged.vcf.gz child.vcf.gz dad.vcf.gz mom.vcf.gz
gunzip -c merged.vcf.gz | awk 'BEGIN {FS = OFS = "\t"} {if($10==".") $10="0/0"; print }'  | awk 'BEGIN {FS = OFS = "\t"} {if($11==".") $11="0/0"; print }' | awk 'BEGIN {FS = OFS = "\t"} {if($12==".") $12="0/0"; print }' | ${RTG} bgzip - > merged_dotashomref.vcf.gz && ${RTG} index -f vcf merged_dotashomref.vcf.gz
mendelian "pre-harmonization, genome-wide" merged_dotashomref.vcf.gz inconsistent.vcf.gz childphased.vcf.gz
report "pre-harmonization, high-confidence child" inconsistent.vcf.gz merged_dotashomref.vcf.gz ${HC}
report "pre-harmonization, high-confidence intersection" inconsistent.vcf.gz merged_dotashomref.vcf.gz ${IHC}

# Measure Mendelianness after harmonization
[[ -r recode-merged.vcf.gz ]] || ${RTG} vcfmerge -o - recode-*/sample.vcf.gz | ${RTG} vcfsubset --remove-infos -i - -o recode-merged.vcf.gz
gunzip -c recode-merged.vcf.gz | awk 'BEGIN {FS = OFS = "\t"} {if($10==".") $10="0/0"; print }'  | awk 'BEGIN {FS = OFS = "\t"} {if($11==".") $11="0/0"; print }' | awk 'BEGIN {FS = OFS = "\t"} {if($12==".") $12="0/0"; print }' | ${RTG} bgzip - > recode-merged_dotashomref.vcf.gz && ${RTG} index -f vcf recode-merged_dotashomref.vcf.gz
mendelian "post-harmonization, genome-wide" recode-merged_dotashomref.vcf.gz recode-inconsistent.vcf.gz recode-childphased.vcf.gz
report "post-harmonization, high-confidence child" recode-inconsistent.vcf.gz recode-merged_dotashomref.vcf.gz ${HC}
report "post-harmonization, high-confidence intersection" recode-inconsistent.vcf.gz recode-merged_dotashomref.vcf.gz ${IHC}

# Construct explicit VCF for most stringent inconsistent set
[[ -r recode-inconsistent-hc-all.vcf.gz ]] || ${RTG} vcffilter -i recode-inconsistent.vcf.gz --bed-regions ${IHC} -o recode-inconsistent-hc-all.vcf.gz

${RTG} vcfstats recode-inconsistent-hc-all.vcf.gz

zgrep '0/0+0/0->0/1' recode-inconsistent-hc-all.vcf.gz | wc -l
zgrep -v '^#\|0/0+0/0->0/1' recode-inconsistent-hc-all.vcf.gz | awk 'BEGIN {FS = OFS = "\t"} { print $1, $2-50, $2+length($4)+50} ' | sed 's/^chr//;s/^X/23/;s/^Y/24/;s/^M/25/' | sort -k1,1n -k2,2n | sed 's/^23/X/;s/^24/Y/;s/^25/M/' > recode-inconsistent-hc-all-nodenovo-slop50-a.bed
chrpref=`zgrep ^chr recode-inconsistent-hc-all.vcf.gz | wc -l` || echo "ignore 0 output as error"
if [ $chrpref -gt 0 ]; then
    sed 's/^/chr/' recode-inconsistent-hc-all-nodenovo-slop50-a.bed > recode-inconsistent-hc-all-nodenovo-slop50.bed
fi
if [ $chrpref -lt 1 ]; then
    mv recode-inconsistent-hc-all-nodenovo-slop50-a.bed recode-inconsistent-hc-all-nodenovo-slop50.bed
fi

subtractBed -a child.bed -b recode-inconsistent-hc-all-nodenovo-slop50.bed > child_noinconsistent.bed
subtractBed -a mom.bed -b recode-inconsistent-hc-all-nodenovo-slop50.bed > mom_noinconsistent.bed
subtractBed -a dad.bed -b recode-inconsistent-hc-all-nodenovo-slop50.bed > dad_noinconsistent.bed

#
#phase transfer for child
#

# Path to calls to transfer phasing onto
CALLS=vcfchildgz.vcf.gz

# Paths to source phasing set
SP=recode-childphased.vcf.gz

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
# existing phasing and phase set information into new field IFT and IPS.
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
${RTG} vcfstats ${SP} 
remove_unphased_hets ${SP} clean-sp.vcf.gz child
${RTG} vcfstats clean-sp.vcf.gz

# Prepare call set by copying any existing phasing into IGT, IPS; correct problem with missing IPS field using sed
report_phasing "Original calls" ${CALLS}
${RTG} vcfstats ${CALLS}

	remove_phasing ${CALLS} unphased-calls-IPSerr.vcf.gz INTEGRATION
#	gunzip -c unphased-calls-IPSerr.vcf.gz | sed 's/\(.\/.$\)/\1:\./' | ${RTG} bgzip - > unphased-calls.vcf.gz && ${RTG} index -f vcf unphased-calls.vcf.gz
	gunzip -c unphased-calls-IPSerr.vcf.gz | awk 'BEGIN {FS = OFS = "\t"} {if(!($10 ~ /:.*:.*:.*:.*:.*:/) && !($1 ~ /^#/)) $10=$10":.";  print }' | ${RTG} bgzip - > unphased-calls.vcf.gz && ${RTG} index -f vcf unphased-calls.vcf.gz
	zgrep -nm 4 IPS  unphased-calls.vcf.gz 

# Transfer phasing from Segregation Phasing
phase_transfer clean-sp.vcf.gz unphased-calls.vcf.gz phase-transfer-sp child,INTEGRATION
report_phasing "Phase Transfered" phase-transfer-sp/calls.vcf.gz


# Strip out unwanted attributes produced during phase transfer
${RTG} vcfsubset -i phase-transfer-sp/calls.vcf.gz -o - --remove-format OGT --remove-info SYNC,CALL,CALL_WEIGHT,PHASE | grep -v '^##CL\|^##RUN-ID' | ${RTG} bgzip - >filtered-phase-transfer.vcf.gz && ${RTG} index -f vcf filtered-phase-transfer.vcf.gz

echo "Final trio result in filtered-phase-transfer.vcf.gz"

${RTG} vcffilter -i filtered-phase-transfer.vcf.gz -o - --javascript "ensureFormatHeader(\"##FORMAT=<ID=PS,Number=1,Type=String,Description=\\\"Phase set for GT\\\">\"); function record() {if(INTEGRATION.GT==\"1/1\") { INTEGRATION.IPS=\".\"; INTEGRATION.PS=\"HOMVAR\"; INTEGRATION.GT=\"1|1\";} else {if((INTEGRATION.GT==\"0/1\" || INTEGRATION.GT==\"1/2\" || INTEGRATION.GT==\"2/1\" || INTEGRATION.GT==\"1/0\") ) {if(INTEGRATION.IPS.length>1) {INTEGRATION.PS=INTEGRATION.IPS; INTEGRATION.GT=INTEGRATION.IGT;} else {INTEGRATION.PS=\".\";};} else { if((INTEGRATION.IPS.length<2)) { INTEGRATION.IPS=\".\";} INTEGRATION.PS=\"PATMAT\";};};}" | awk 'BEGIN {FS = OFS = "\t"} {if(!($10 ~ /:.*:.*:.*:.*:.*:/) && !($1 ~ /^#/)) $10=$10":.:.";  print }' | ${RTG} bgzip - > filtered-phase-transfer-allphasing.vcf.gz && ${RTG} index -f vcf filtered-phase-transfer-allphasing.vcf.gz
report_phasing "Phase Transferred plus original phasing" filtered-phase-transfer-allphasing.vcf.gz

${RTG} vcffilter -i filtered-phase-transfer-allphasing.vcf.gz -o - --keep-expr 'INTEGRATION.GT.indexOf("|")>0&&!has(INTEGRATION.PS)' --no-header | wc -l || echo "ignore stderr here"
${RTG} vcffilter -i filtered-phase-transfer-allphasing.vcf.gz -o - --keep-expr 'INTEGRATION.GT.indexOf("|")==-1&&has(INTEGRATION.PS)' --no-header | wc -l || echo "ignore stderr here"
${RTG} vcffilter -i filtered-phase-transfer-allphasing.vcf.gz -o - --keep-expr 'has(INTEGRATION.PS)&&INTEGRATION.PS!="PATMAT"&&INTEGRATION.PS!="HOMVAR"' | ${RTG} vcffilter -i - --javascript 'function record() {print(CHROM+":"+INTEGRATION.PS)}' | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c || echo "ignore stderr here"

#
# Upload results
#

file_idoutputphasedvcf=`dx upload filtered-phase-transfer-allphasing.vcf.gz -o "$prefix"_triophased.vcf.gz --brief`
dx-jobutil-add-output "outputphasedvcf" "$file_idoutputphasedvcf"

file_idoutputphasedvcftbi=`dx upload filtered-phase-transfer-allphasing.vcf.gz.tbi -o "$prefix"_triophased.vcf.gz.tbi --brief`
dx-jobutil-add-output "outputphasedvcftbi" "$file_idoutputphasedvcftbi"

file_idoutputinconsistentvcf=`dx upload recode-inconsistent-hc-all.vcf.gz -o "$prefix"_trioinconsistent.vcf.gz --brief`
dx-jobutil-add-output "outputinconsistentvcf" "$file_idoutputinconsistentvcf"

file_idoutputinconsistentvcftbi=`dx upload recode-inconsistent-hc-all.vcf.gz.tbi -o "$prefix"_trioinconsistent.vcf.gz.tbi --brief`
dx-jobutil-add-output "outputinconsistentvcftbi" "$file_idoutputinconsistentvcftbi"

file_idoutputinconsistentbed=`dx upload recode-inconsistent-hc-all-nodenovo-slop50.bed -o "$prefix"_trioinconsistent_nodenovo_slop50.bed --brief`
dx-jobutil-add-output "outputinconsistentbed" "$file_idoutputinconsistentbed"

name=`dx describe "$childbed" --name`
name="${name%.bed}"
file_idoutputnewchildbed=`dx upload child_noinconsistent.bed -o "$name"_noinconsistent.bed --brief`
dx-jobutil-add-output "outputnewchildbed" "$file_idoutputnewchildbed"

name=`dx describe "$dadbed" --name`
name="${name%.bed}"
file_idoutputnewdadbed=`dx upload dad_noinconsistent.bed -o "$name"_noinconsistent.bed --brief`
dx-jobutil-add-output "outputnewdadbed" "$file_idoutputnewdadbed"

name=`dx describe "$mombed" --name`
name="${name%.bed}"
file_idoutputnewmombed=`dx upload mom_noinconsistent.bed -o "$name"_noinconsistent.bed --brief`
dx-jobutil-add-output "outputnewmombed" "$file_idoutputnewmombed"

