#!/usr/bin/bash
## GIAB integration variant calling pipeline for Illumina data

### Setting default parameter values
PLATFORM="Illumina"
MAPPER="novoalign"


## Command Line Arguments
# usage()
# {
#     echo "usage: giab_illumina_variant_calling_pipeline.sh PARAMETERS"
# }

while [ "$1" != "" ]; do
    case $1 in
        --hg )              shift
                            HG=$1
                            ;;
        --dataset )         shift
                            DATASETID=$1
                            ;;
        --bamurl )          shift
                            BAMURL=$1
                            ;;
        --baiurl )          shift
                            BAIURL=$1
                            ;;
        --refid )           shift
                            REFID=$1
                            ;;
        --mapper )          shift
                            MAPPER=$1
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
        --maxdepth )        shift
                            MAXDEPTH=$1
                            ;;
        # -h | --help )       usage
        #                     exit
        #                     ;;
        * )                 echo "Check Parameters: $1" #usage
                            exit 1
    esac
    shift
done

## Defining parameters
ROOTDIR=${HG}/${REFID}/${PLATFORM}/${PLATFORM}_${REFID}_${DATASETID}

if [ ${REFID} = "GRCh37" ]; then
  GENOME=/assets/hs37d5.fa.gz
  GENOMEIDX=/assets/hs37d5.fa.fai
  REF=/assets/hs37d5.fasta-index.tar.gz
elif [ ${REFID} = "GRCh38" ]; then
  GENOME=/assets/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa.gz
  GENOMEIDX=/assets/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa.fai
  REF=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz
else
  echo "--refid must be GRCh37 or GRCh38"
fi


## Output file prefix
PREFIX=${HG}_${REFID}_${DATASETID}_${MAPPER}

if [[ ${REFID} = "GRCh37" && ${DATASETID} = "6Kb_MatePair"  ]]; then

  IMPORTJOBID=$(dx run -y --brief \
    GIAB:/Workflow/samtools_import_reheader_splitchrom_addrg_reord \
    -iurlbam=${BAMURL} \
    -iurlbai=${BAIURL} \
    -iprefix=${PREFIX} \
    -irgid=1 -irglb=all -irgpl=illumina -irgpu=all -irgsm=${HG} \
    --destination=${ROOTDIR})
else
  ## Download bam and split by chromosome
  CHROMSPLITS=([1]="1to5" [2]="6to12" [3]="13toMT")

  ## Use array for saving jobids
  declare -a IMPORTJOBIDS
  for i in 1 2 3;
    do
      CHR="_withchr"
      if [ ${REFID} = "GRCh37" ]; then
        CHR=""
      fi
      echo ${CHR}

      JOBID=$(dx run -y --brief \
        GIAB:/Workflow/samtools_import_splitchrom_addrg_${CHROMSPLITS[${i}]}${CHR} \
        -iurlbam=${BAMURL} -iurlbai=${BAIURL} \
        -iprefix=${PREFIX} \
        -irgid=1 -irglb=all -irgpl=illumina -irgpu=all -irgsm=${HG} \
        --destination=${ROOTDIR} \
        --instance-type=mem2_hdd2_x4)
        IMPORTJOBIDS+=([${i}]=${JOBID})
  done
fi

## Variant calling run after split chrom finishes
for i in {1..22} MT X Y;
  do

    if [ ${REFID} = "GRCh37" ]; then
      CHROM=${i}
    elif [ ${REFID} = "GRCh38" ]; then
      CHROM=chr${i}

      ### For chrom MT
      if [[ ${i} =~ MT ]]; then
        CHROM=chrM
      fi

    else
      echo "--refid must be GRCh37 or GRCh38"
    fi



    PREFIX=${HG}_${i}_${REFID}_${MAPPER}_${DATASETID}

    ## define BAMJOBID by chromosome and dataset
    if [[ ${REFID} = "GRCh37" && ${DATASETID} = "6Kb_MatePair"  ]]; then
      BAMJOBID=${IMPORTJOBID}
    elif [[ ${i} =~ [MTXY] ]]; then
      BAMJOBID=${IMPORTJOBIDS[3]}
    elif [ ${i} -le 5 ]; then
      BAMJOBID=${IMPORTJOBIDS[1]}
    elif [ ${i} -gt 5 ] && [ ${i} -le 12 ]; then
      BAMJOBID=${IMPORTJOBIDS[2]}
    else
      BAMJOBID=${IMPORTJOBIDS[3]}
    fi


    ## Freebayes variant calling
    # Not sure if -itargets_bed is for 37 and 38
    # dx run -y --depends-on ${BAMJOBID} \
    #   freebayes \
    #   -isorted_bams=${BAMJOBID}:bam${i} \
    #   -ioutput_prefix=${PREFIX}_FB \
    #   -igenotype_qualities=TRUE \
    #   -istandard_filters=FALSE \
    #   -iadvanced_options="-F 0.05 -m 0" \
    #   -igenome_fastagz=${GENOME} \
    #   -itargets_bed="GIAB:/Workflow/Chromosome_bed_files/${i}.bed" \
    #   --destination=${ROOTDIR}/FreeBayes_output/


    ## Callable Loci
    dx run -y --depends-on ${BAMJOBID} \
      GIAB:/Workflow/GATK_V3.5/gatk-callableloci-v3.5-anyref \
      -iinput_bam=${BAMJOBID}:bam${i} \
      -iinput_bai=${BAMJOBID}:bai${i} \
      -ioutput_prefix=${PREFIX}_callableloci \
      -iref=${REF} \
      -iextra_options="-L ${CHROM} -minDepth 20 -mmq 20 -maxDepth ${MAXDEPTH}" \
      --destination=${ROOTDIR}/CallableLoci_output/


    # ## Sentieon
    # JOBIDSNT=$(dx run -y --brief --depends-on ${BAMJOBID} \
    #   GIAB:sentieon-haplotyper-gvcf-reheadunsorted \
    #   -isorted_bam=${BAMJOBID}:bam${i} \
    #   -isorted_bai=${BAMJOBID}:bai${i}\
    #   -ioutput_prefix=${PREFIX}_sentieonHC_gvcf  \
    #   -igenome_fasta=${GENOME} \
    #   -igenome_fastaindex=${GENOMEIDX} \
    #   -iextra_driver_options="--interval ${CHROM}" \
    #   -iextra_algo_options="--call_conf 2 --emit_conf 2" \
    #   -ilicense_server_file=/Workflow/sentieon_license_server.info \
    #   --destination=${ROOTDIR}/Sentieon_output/)
    #
    #   ## Post processing Sentieon variant calls
    #   VCF=${JOBIDSNT}:gvcfgz
    #   VCFIDX=${JOBIDSNT}:tbi
    #
    #   ## GATK genotype gcvf
    #   dx run -y --depends-on ${JOBIDSNT} \
    #     GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5-anyref \
    #     -ivcfs=${VCF} \
    #     -ivcfs=${VCFIDX} \
    #     -iprefix=${PREFIX}_sentieonHC \
    #     -iref=${REF} \
    #     --destination=${ROOTDIR}/Sentieon_output/
    #
    #
    #   ## Integration prep
    #   dx run -y --depends-on ${JOBIDSNT} \
    #     GIAB:/Workflow/integration-prepare-gatkhc-v3.3.2-anyref \
    #     -igvcf=${VCF} \
    #     -igvcftbi=${VCFIDX} \
    #     -iref=${REF} \
    #     -ichrom=${CHROM} \
    #     --destination=${ROOTDIR}/Integration_prepare_sentieon_v.3.3.2/

    ## Cleanup - remove chromosome bams after variant calling finishes
    # dx wait ${JOBIDFB} ${JOBINCL} ${JOBIDSNT} && \
    #   dx rm ${BAMJOBID}:bam${i} ${BAMJOBID}:bai${i}
done
