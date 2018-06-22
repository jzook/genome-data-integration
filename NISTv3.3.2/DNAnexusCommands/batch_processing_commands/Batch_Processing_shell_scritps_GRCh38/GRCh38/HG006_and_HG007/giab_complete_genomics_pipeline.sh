#!/usr/bin/bash
## Script for preparing complete genomics data for integration pipeline

### Setting default parameter values
PLATFORM="Complete_Genomics"
REFID="GRCh37"
WORKDIR="."

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
        --vcfurl )          shift
                            VCFURL=$1
                            ;;
        --cgvcf )           shift
                            CGVCF=$1
                            ;;
        --refid )           shift
                            REFID=$1
                            ;;
        --liftover )        shift
                            LIFTOVERCHAIN=$1
                            ;;
        --queryref )        shift
                            QUERYREF=$1
                            ;;
        --targetref )       shift
                            TARGETREF=$1
                            ;;
        --genomewarp )      shift
                            GENOMEWARP=$1
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

## Defining parameters
ROOTDIR=${HG}/${REFID}/${PLATFORM}
INPUTVCF=${ROOTDIR}/${CGVCF}.vcf.bz2

################################################################################
############## Upload data to DNAnexus
dx mkdir -p ${ROOTDIR} 

UPLOADJOBID=$(dx run -y --brief url_fetcher -iurl=${VCFURL})


dx wait ${UPLOADJOBID} && dx mv GIAB:/${CGVCF}.vcf.bz2 GIAB:/${INPUTVCF}

################################################################################
############## run integration prepare CG

## Use array for saving jobids
declare -a PREPJOBIDS
CHROMARRAY=( {1..22} X ) 
if [ ${HASY} = true ]; then
    CHROMARRAY+=(Y)
fi

test2=(${test1[@]} Y)
for i in ${CHROMARRAY[@]};
  do
    JOBID=$(dx run -y --brief \
        GIAB:/Workflow/integration-prepare-cg \
        -ivcf_in=${INPUTVCF} \
        -ichrom=${i} \
        --destination=${ROOTDIR}/Integration_prepare_cg_output/)
    PREPJOBIDS+=(${JOBID})
done


################################################################################
############## combine GCRh37 vcfs
## Combining input values
DEPENDIDS=""
CHROMVCFS=""
CHROMBEDS=""
for i in ${PREPJOBIDS[@]};
  do
    DEPENDIDS="${DEPENDIDS} --depends-on  ${i}"
    CHROMVCFS="${CHROMVCFS} -ivcfs=${i}:outvcfgz"
    CHROMBEDS="${CHROMBEDS} -ibeds=${i}:outcallablebed"
  done

## Output prefix
COMBINEDPREFIX=${ROOTDIR}/Integration_prepare_cg_output/${HG}_GRCh37_CHROM1-Y_${CGVCF}

############## combine GCRh37 vcfs

COMBINEDVCFJOBID=$(dx run -y --brief ${DEPENDIDS} \
  Workflow/vcf-combineallchrom \
  ${CHROMVCFS} \
  -iprefix=${COMBINEDPREFIX})


############## combine GCRh37 beds
COMBINDBEDJOBID=$(dx run -y --brief ${DEPENDIDS} \
  Workflow/bed-combineallchrom \
  ${CHROMBEDS} \
  -iprefix=${COMBINEDPREFIX})

################################################################################
############## Liftover uses verily genomewarp - Offline

############## Liftover variables
GRCh37ROOT=${HG}_GRCh37_${CGVCF}

# VCF
GRCh37VCFGZ=${GRCh37ROOT}.vcf.gz
GRCh37VCF=${GRCh37ROOT}.vcf
GRCh37HEADER=${GRCh37ROOT}_header.vcf
GRCh37NOHEADER=${GRCh37ROOT}_NOheader.vcf
GRCh37NOHEADERCHROMFIX=${GRCh37ROOT}_NOheader_CHROMfixed.vcf
GRCh37VCFCHROMFIX=${GRCh37ROOT}_CHROMfixed.vcf

# BED
GRCh37BED=${GRCh37ROOT}_callable.bed
GRCh37BEDCHROMFIX=${GRCh37ROOT}_callable_CHROMfixed.bed


############## Download VCF and BED
dx wait ${COMBINEDVCFJOBID} && dx download -o ${GRCh37VCF} ${COMBINEDPREFIX}.vcf.gz
# dx download -o ${GRCh37VCFGZ} ${COMBINEDPREFIX}.vcf.gz
dx wait ${COMBINEDBEDJOBID} && dx download -o ${GRCh37BED} ${COMBINEDPREFIX}.bed
# dx download -o ${GRCh37BED} ${COMBINEDPREFIX}.bed


############# convert CHROM number for GRCh37 highconf files

##Parse out vcf header
gunzip ${GRCh37VCFGZ}
grep ^# ${GRCh37VCF} > ${GRCh37HEADER}

##create no header file
sed '/^#/ d' ${GRCh37VCF} > ${GRCh37NOHEADER}

## replace chrom # with chr# at start of every row in no header file
sed 's/^/chr/' ${GRCh37NOHEADER}  > ${GRCh37NOHEADERCHROMFIX}

## replace chrom # with chr# at start of every row in bed
sed 's/^/chr/' ${GRCh37BED} > ${GRCh37BEDCHROMFIX}

## add header back in
cat ${GRCh37HEADER} ${GRCh37NOHEADERCHROMFIX} > ${GRCh37VCFCHROMFIX}

## Cleanup
rm ${GRCh37VCF} ${GRCh37HEADER} ${GRCh37NOHEADER} ${GRCh37NOHEADERCHROMFIX} ${GRCh37BED}

###### Perform liftover (note: fasta's cannot be zipped)

for i in ${CHROMARRAY[@]};
  do

    # Separate out each autosome and the X
    
    egrep "^(#|(chr)?${i}[[:space:]])" ${GRCh37VCFCHROMFIX} > GRCh37_chr${i}.vcf
    egrep "^(chr)?${i}[[:space:]]" ${GRCh37BEDCHROMFIX} > GRCh37_chr${i}.bed

    LIFTOVERVCF=${HG}_CG_chr${i}_LIFTOVER.vcf.gz
    LIFTOVERBED=${HG}_CG_chr${i}_LIFTOVER.bed


    #liftover (note: fasta's cannot be zipped)
    java -jar -Xmx10g ${GENOMEWARP}/verilylifesciences-genomewarp-1.0.0-runnable.jar \
      --lift_over_chain_path ${LIFTOVERCHAIN} \
      --raw_query_vcf GRCh37_chr${i}.vcf \
      --raw_query_bed GRCh37_chr${i}.bed \
      --ref_query_fasta  ${QUERYREF} \
      --ref_target_fasta  ${TARGETREF} \
      --work_dir ${WORKDIR} \
      --output_variants_file ${LIFTOVERVCF} \
      --output_regions_file ${LIFTOVERBED}

    convGRCh38VCF=convGRCh38.vcf
    GRCh38BED=${HG}_${i}_convGRCh38_CG_${CGVCF}.bed

    #remove extra contigs so only left with chr${i}
    egrep "^(#|(chr)?${i}[[:space:]])" ${LIFTOVERVCF} > ${convGRCh38VCF}
    egrep "^(chr)?${i}[[:space:]]" ${LIFTOVERBED} > ${GRCh38BED}

    ## Sort, index, and zip VCF file
    GRCh38VCF=${HG}_${i}_convGRCh38_CG_${CGVCF}_sorted.vcf.gz
    GRCh38VCFIDX=${HG}_${i}_convGRCh38_CG_${CGVCF}_sorted.vcf.gz.tbi

    sed '/^#/ d' ${convGRCh38VCF} | sort -k2,2n | bgzip -c > ${GRCh38VCF}
    tabix -f ${GRCh38VCF}

    ## upload liftover vcfs and bed to DNAnexus
    dx upload --destination ${HG}/GRCh38/${PLATFORM} ${GRCh38VCF} ${GRCh38VCFIDX} ${GRCh38BED}

    # Cleanup
    rm GRCh37_chr${i}.vcf GRCh37_chr${i}.bed ${LIFTOVERVCF} ${LIFEOVERBED} ${convGRCh38VCF}
done
