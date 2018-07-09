#!/usr/bin/bash
## Script for running GIAB variant calling pipeline on 10X data
set -v

### Setting default parameter values
PLATFORM=10XGenomics

## Command Line Arguments
# usage()
# {
#     echo "usage: giab_complete_genomics_pipeline.sh PARAMETERS"
# }

while [ "$1" != "" ]; do
    case $1 in
        --hg )              shift
                            HG=$1
                            ;;
        --bamurl )          shift
                            BAMURL=$1
                            ;;
        --refid )           shift
                            REFID=$1
                            ;;
        --genome )          shift
                            GENOME=$1
                            ;;
        --genomeidx )       shift
                            GENOMEIDX=$1
                            ;;
        --ref )             shift
                            REF=$1
                            ;;
        --rmsg )            shift
                            RGSM=$1
                            ;;
        --hasY )            HASY=TRUE
                            ;;
        # -h | --help )       usage
        #                     exit
        #                     ;;
        * )                 shift
                            echo $1 #usage
                            exit 1
    esac
    shift
done

## Defining variables
ROOTDIR=${HG}/${REFID}/${PLATFORM}/

if [ ${REFID} = "GRCh37" ]; then
  GENOME=/assets/hg19_10X_2.1.0.fa.gz
  GENOMEIDX=/assets/hg19_10X_2.1.0.fa.fai
  REF=/assets/hg19_10X_2.1.0.fasta-index.tar.gz
elif [ ${REFID} = "GRCh38" ]; then
  GENOME=/assets/GRCh38_10X_2.1.0.fa.gz
  GENOMEIDX=/assets/GRCh38_10X_2.1.0.fa.fai
  REF=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz
else
  echo "--refid must be GRCh37 or GRCh38"
fi



################################################################################
############## Upload data to DNAnexus
dx mkdir -p ${ROOTDIR}
UPLOADJOBID=$(dx run -y --brief --destination ${ROOTDIR} url_fetcher -iurl=${BAMURL})


################################################################################
############## split haplotypes
BAMPREFIX=${HG}_${REFID}_10X
HAPSPLITJOBID=$(dx run -y --brief --depends-on ${UPLOADJOBID} \
  /Workflow/split_bam_haplotype \
  -ibam=${UPLOADJOBID}:file \
  -iprefix=${BAMPREFIX} \
  --destination=${ROOTDIR})


################################################################################
############## split by chromsome

## Use array for saving jobids
declare -a SPLITJOBIDS
for i in 1 2;
  do
    ## Sort and index Haplotype bam
    SORTID=$(dx run -y --brief --depends-on ${HAPSPLITJOBID} \
      samtools_sort \
      -imappings_bam=${HAPSPLITJOBID}:HP${i}\
      --destination=${ROOTDIR})

    IDXID=$(dx run -y --brief --depends-on ${SORTID} \
      samtools_index \
      -isorted_bam=${SORTID}:sorted_bam \
      --destination=${ROOTDIR})

    ## Split haplotype bam by chromosome - not working
    JOBID=$(dx run -y --brief \
      --depends-on ${SORTID} --depends-on ${IDXID}\
      GIAB:/Workflow/samtools_splitchrom_addrg_withchr \
      -isorted_bam=${SORTID}:sorted_bam \
      -iindex_bai=${IDXID}:index_bai \
      -iprefix=${BAMPREFIX}_HP${i}_ \
      -irgid=HP${i} -irglb=10X -irgpl=illumina -irgpu=all -irgsm=${RGSM} \
      --destination=${ROOTDIR} \
      --instance-type=mem2_hdd2_x1 )
    SPLITJOBIDS+=([${i}]=${JOBID})

done

################################################################################
############## Variant calling and integration

# Use X Y chrom flag to prevent Y Chrom issues for mother
CHROMARRAY=( {1..22} X )
if [[ ${HASY} = true ]]; then
    CHROMARRAY+=(Y)
fi


for i in 20; # ${CHROMARRAY[@]};
  do
    declare -a JOBIDSSNT
    declare -a JOBIDSCL
    ## Variant calling for individual haplotypes
    for j in 1 2;
      do
        ## variables
        PREFIX=${HG}_${j}_10X_HP${i}

        JOBID=$(dx run -y --brief --depends-on ${SPLITJOBIDS[${j}]} \
          GIAB:sentieon-haplotyper-gvcf-reheadunsorted \
          -isorted_bam=${SPLITJOBIDS[${j}]}:bam${i} \
          -isorted_bai=${SPLITJOBIDS[${j}]}:bai${i}\
          -ioutput_prefix=${PREFIX}_sentieonHC_gvcf  \
          -igenome_fasta=${GENOME} \
          -igenome_fastaindex=${GENOMEIDX} \
          -iextra_driver_options="--interval chr${i}" \
          -iextra_algo_options="--call_conf 2 --emit_conf 2" \
          -ilicense_server_file=/Workflow/sentieon_license_server.info \
          --destination ${ROOTDIR}/Sentieon_output/)
          JOBIDSSNT+=([${j}]=${JOBID})

          ## GATK callableLoci
          JOBID=$(dx run -y --brief --depends-on ${SPLITJOBIDS[${j}]} \
            GIAB:/Workflow/GATK_V3.5/gatk-callableloci-v3.5-anyref \
              -iinput_bam=${SPLITJOBIDS[${j}]}:bam${i} \
              -iinput_bai=${SPLITJOBIDS[${j}]}:bai${i} \
              -ioutput_prefix=${PREFIX}_callableloci \
              -iref=${REF} \
              -iextra_options="-L chr${i} -minDepth 20 -mmq 20 -maxDepth 566" \
              --destination ${ROOTDIR}/CallableLoci_output/)
          JOBIDSCL+=([${j}]=${JOBID})
      done

      ## Integration Prepare
      dx run -y \
        --depends-on ${JOBIDSSNT[1]} --depends-on ${JOBIDSSNT[2]} \
        --depends-on ${JOBIDSCL[1]}  --depends-on ${JOBIDSCL[2]}  \
        GIAB:/Workflow/integration-prepare-10X-v3.3-anyref \
        -igvcf1=${JOBIDSSNT[1]}:gvcfgz \
        -igvcftbi1=${JOBIDSSNT[1]}:tbi \
        -igvcf2=${JOBIDSSNT[2]}:gvcfgz \
        -igvcftbi2=${JOBIDSSNT[2]}:tbi \
        -ibed1=${JOBIDSCL[1]}:bed_file \
        -ibed2=${JOBIDSCL[2]}:bed_file \
        -iprefix=${HG}_${i}_${REFID}_10X_sentieonHCbyhaplo \
        -iref=${REF} \
        -ichrom=chr${i} -imaxcov=20 \
        --destination ${ROOTDIR}/Integration_prepare_10X_output_v3.3/
    done
