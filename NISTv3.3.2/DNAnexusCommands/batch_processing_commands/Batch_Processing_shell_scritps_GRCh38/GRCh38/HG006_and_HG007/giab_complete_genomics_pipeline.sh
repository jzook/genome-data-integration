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
                            cgvcf=$1
        --refid )           shift
                            REFID=$1
                            ;;
        --liftover )        shift
                            LIFEOVERCHAIN=$1
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
        # -h | --help )       usage
        #                     exit
        #                     ;;
        * )                 OTHER=$1 #usage
                            exit 1
    esac
    shift
done

## Defining parameters
ROOTDIR=${HG}/${REFID}/${PLATFORM}
INPUTVCF=${ROOTDIR}/${CGVCF}

################################################################################
############## Upload data to DNAnexus

UPLOADJOBID=$(dx run url_fectcher -iurl=${VCFURL} -ioutput_name=${INPUTVCF})

################################################################################
############## run integration prepare CG

## Use array for saving jobids
declare -a PREPJOBIDS
for i in {1 .. 22} X Y;
  do
    JOBID=$(dx wait ${UPLOADJOBID} && \
      dx run -y \
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
CHROMCALLABLEBEDS=""
CHROMNOTCALLABLEBEDS=""
for i in PREPJOBIDS;
  do
    DEPENDIDS="${DEPENDIDS} --depends-on  ${i}"
    CHROMVCFS="${CHROMVCFS} -ivcfs=${i}:outvcfgz"
    CHROMCALLABLEBEDS="${CHROMCALLABLEBEDS} -ibeds==${i}:outcallablebed"
  done

## Output prefix
COMBINEDPREFIX=${HG}_GRCh37_CHROM1-Y_${CGVCF}

############## combine GCRh37 vcfs

COMBINEDVCFJOBID=$(dx run -y ${DEPENDIDS} \
  Workflow/vcf-combineallchrom \
  ${CHROMVCFS} \
  -iprefix=${COMBINEDPREFIX})


############## combine GCRh37 beds
COMBINDBEDJOBID=$(dx run -y ${DEPENDIDS} \
  Workflow/bed-combineallchrom \
  ${CHROMCALLABLEBEDS} \
  -iprefix=${COMBINEDPREFIX})



################################################################################
############## Liftover uses verily genomewarp - Offline

############## Liftover variables
GRCh37ROOT=${HG}_GRCh37_${CGVCF}

# VCF
GRCh37VCF=${GRCh37ROOT}.vcf
GRCh37HEADER=${GRCh37ROOT}_header.vcf
GRCh37NOHEADER=${GRCh37ROOT}_NOheader.vcf
GRCh37NOHEADERCHROMFIX=${GRCh37ROOT}_NOheader_CHROMfixed.vcf
GRCH37VCFCHROMFIX=${GRCh37ROOT}_CHROMfixed.vcf

# BED
GRCh37BED=${GRCh37ROOT}_callable.bed
GRCh37BEDCHROMFIX=${GRCh37ROOT}_callable_CHROMfixed.bed


############## Download VCF
dx wait ${COMBINEDVCFJOBID} && dx download -o ${GRCh37VCF} ${COMBINEDVCF}
dx wait ${COMBINEDBEDJOBID} && dx download -o ${GRCh37BED} ${COMBINEDBED}




############## convert CHROM number for GRCh37 highconf files

# Parse out vcf header
grep ^# ${GRCh37VCF} > ${GRCh37HEADER}

#create no header file
sed '/^#/ d' ${GRCh37VCF} > ${GRCh37NOHEADER}

# replace chrom # with chr# at start of every row in no header file
sed 's/^/chr/' ${GRCh37NOHEADER}  > ${GRC37NOHEADERCHROMFIX}

# replace chrom # with chr# at start of every row in bed
sed 's/^/chr/' ${GRCh37BED} > ${GRCh37BEDCHROMFIX}

# add header back in
cat ${GRCh37HEADER} ${GRCh37NOHEADERCHROMFIX} > ${GRCH37VCFCHROMFIX}

## Cleanup
rm ${GRCh37VCF} ${GRCh37NOHEADER} ${GRC37NOHEADERCHROMFIX} ${GRCh37BED}

###### Perform liftover (note: fasta's cannot be zipped)

for i in {1..22} X Y;
  do

    # Separate out each autosome and the X
    egrep "^(#|(chr)?${i}[[:space:]])" ${GRCh37VCFCHROMFIX} > GRCh37_chr${i}.vcf
    egrep "^(chr)?${i}[[:space:]]" ${GRCh37BEDCHROMFIX} > GRCh37_chr${i}.bed

    LIFTOVERVCF=${HG}_CG_chr${i}_LIFTOVER.vcf
    LIFTOVERBED=${HG}_CG_chr${i}_LIFTOVER.bed


    #liftover (note: fasta's cannot be zipped)
    java -jar -Xmx10g ${GENOMEWARP}/verilylifesciences-genomewarp-1.0.0-runnable.jar \
      --lift_over_chain_path ${LIFTOVERCHAIN} \
      --raw_query_vcf GRCh37_chr${i}.vcf \
      --raw_query_bed GRCh37_chr${i}.bed \
      --ref_query_fasta  ${QUERYREF}\
      --ref_target_fasta  ${TARGETREF}\
      --work_dir ${WORKDIR} \
      --output_variants_file ${LIFTOVERVCF} \
      --output_regions_file ${LIFEOVERBED}

    #remove extra contigs so only left with chr${i}
    egrep "^(#|(chr)?${i}[[:space:]])" ${LIFTOVERVCF} > ${convGRCh38VCF}
    egrep "^(chr)?${i}[[:space:]]" ${LIFTOVERBED} > ${GRCh38BED}

    ## Sort, index, and zip VCF file
    GRCh38VCF=${HG}_${i}_convGRCh38_CG_${CGVCF}_sorted.vcf.gz
    GRCh38VCFIDX=${HG}_${i}_convGRCh38_CG_${CGVCF}.vcf.tbi

    sed '/^#/ d' ${convGRCh38VCF} | sort -k2,2n | bgzip -c > ${GRCh38VCF}
    tabix -f ${GRCh38VCF}

    ## upload liftover vcfs and bed to DNAnexus
    dx upload --destination ${ROOTDIR} ${GRCh38VCF} ${GRCh38VCFIDX} ${GRCh38BED}

    ## Cleanup
    rm ${LIFTOVERVCF} ${LIFEOVERBED} ${convGRCh38VCF}
done
