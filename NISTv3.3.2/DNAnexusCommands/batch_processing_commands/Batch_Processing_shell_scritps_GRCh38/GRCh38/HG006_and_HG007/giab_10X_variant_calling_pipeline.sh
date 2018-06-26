#!/usr/bin/bash
## Script for running GIAB variant calling pipeline on 10X data


### Setting default parameter values
PLATFORM=10XGenomics
GENOME=/assets/GRCh38_10X_2.1.0.fa.gz
GENOMEIDX=/assets/GRCh38_10X_2.1.0.fa.fai
REF=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz


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


## Download bam
BAM=${HG}_${REFID}_10X.bam
wget ${BAMURL} -o ${BAM}

## split haplotypes (done locally)
HP1BAM=${HG}_${REFID}_10X_HP1
HP2BAM=${HG}_${REFID}_10X_HP2
HP1BAI=${HG}_${REFID}_10X_HP1.bai
HP2BAI=${HG}_${REFID}_10X_HP2.bai

bamtools filter -in ${BAM} -out ${HP2BAM} -tag HP:1
bamtools filter -in ${BAM} -out ${HP2BAM} -tag HP:2

## Upload HP1 and HP1
UPLOADJOBID=$(dx upload --destination ${ROOTDIR} \
                  ${HP1BAM} ${HP2BAM} ${HP1BAI} ${HP2BAI})

## Use array for saving jobids
declare -a SPLITJOBIDS
for i in 1 2;
  do


    ## Split bam by chromosome
    JOBID=$(dx run -y --brief \
      GIAB:/Workflow/samtools_reheader_splitchrom_addrg_reord \
      -isorted_bam=${ROOTDIR}/${BAMPREFIX}_${i}.bam \
      -iindex_bai=${ROOTDIR}/${BAMPREFIX}_${i}.bam.bai \
      -irgid=HP1 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=${RGSM} \
      --destination=${ROOTDIR} \
      --instance-type=mem2_hdd2_x1 )
      SPLITJOBIDS+=([${i}]=${JOBID})
done


## Convert loop to loop over chroms then haplotypes
## Use X Y chrom flag to prevent Y Chrom issues for mother
## variant calling sentieon on individual chroms
for i in {1..22} X Y;
  do
    declare -a JOBIDSSNT
    declare -a JOBIDSCL
    ## Variant calling for individual haplotypes
    for j in 1 2;
    ## variables
    PREFIX=${HG}_${j}_10X_HP${i}

    JOBID=$(dx run -y --brief --depends-on ${SPLITJOBID[${j}]} \
      GIAB:sentieon-haplotyper-gvcf \
      -isorted_bam=${SPLITJOBID[${j}]}:bam${i} \
      -isorted_bai=${SPLITJOBID[${j}]}:bai${i}\
      -ioutput_prefix=${PREFIX}_sentieonHC_gvcf  \
      -igenome_fasta=${GENOME} \
      -igenome_fastaindex=${GENOMEIDX} \
      -iextra_driver_options="--interval chr${i}" \
      -iextra_algo_options="--call_conf 2 --emit_conf 2" \
      -ilicense_server_file=/Workflow/sentieon_license_server.info \
      --destination=${ROOTDIR}/Sentieon_output/)
      JOBIDSSNT+=([${j}]=${JOBID})

      ## GATK callableLoci
      JOBID=$(dx run -y --brief --depends-on ${SPLITJOBID[${j}]}  --brief \
        GIAB:/Workflow/GATK_V3.5/gatk-callableloci-v3.5-anyref \
          -iinput_bam=${SPLITJOBID[${j}]}:bam${i} \
          -iinput_bai=${SPLITJOBID[${j}]}:bai${i} \
          -ioutput_prefix=${PREFIX}_callableloci \
          -iref=${REF} \
          -iextra_options="-L chr${j} -minDepth 20 -mmq 20 -maxDepth 566" \
          --destination=${ROOTDIR}/CallableLoci_output/)
      JOBIDSCL+=([${j}]=${JOBID})
  done

  ## Integration Prepare
  dx run -y \
  --depends-on ${JOBIDSSNT[1]} --depdens-on ${JOBIDSSNT[2]} \
  --depends-on ${JOBIDSCL[1]} --depdens-on ${JOBIDSCL[2]} \
  GIAB:/Workflow/integration-prepare-10X-v3.3-anyref \
    -igvcf1=${JOBIDSSNT[1]}:gcvf \
    -igvcftbi1=${JOBIDSSNT[1]}:tbi \
    -igvcf2=${JOBIDSSNT[2]}:gcvf \
    -igvcftbi2=${JOBIDSSNT[2]}:tbi \
    -ibed1=${JOBIDSCL[1]}:callablebed \
    -ibed2=${JOBIDSCL[2]}:callablebed\
    -iprefix=${HG}_${i}_${REFID}_10X_sentieonHCbyhaplo \
    -iref=${REF} \
    -ichrom=chr1 -imaxcov=20 \
    --destination=${ROOTDIR}/Integration_prepare_10X_output_v3.3/
done
