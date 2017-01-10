

#HG004 Liftover b37 to b38


##### Complete Genomics ######

## 1) need to convert CHROM number for GRCh37 highconf files
# Run 12/01/16
# Parse out vcf header
grep ^# HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM.vcf > HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_header.vcf

#create no header file
sed '/^#/ d' HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM.vcf > HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_NOheader.vcf

# replace chrom # with chr# at start of every row in no header file
sed 's/^/chr/' HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_NOheader.vcf  > HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_NOheader_CHROMfixed.vcf

# replace chrom # with chr# at start of every row in bed
sed 's/^/chr/' HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_callable.bed > HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_callable_CHROMfixed.bed

# add header back in 
cat HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_header.vcf HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_NOheader_CHROMfixed.vcf > HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_CHROMfixed.vcf

## 2) liftover (note: fasta's cannot be zipped)
#Run 12/2/16
INPUT_VCF=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG4_b37_vcfs_to_b38/cg/HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_CHROMfixed.vcf
INPUT_BED=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG4_b37_vcfs_to_b38/cg/HG004_GRCh37_CHROM1-X_vcfBeta-GS000037262-ASM_callable_CHROMfixed.bed
SHARD_DIR=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG4_b37_vcfs_to_b38/cg/
for c in `seq 1 22` X Y
do
  # Separate out each autosome and the X
  egrep "^(#|(chr)?${c}[[:space:]])" $INPUT_VCF > $SHARD_DIR/'GRCh37_chr'${c}.vcf
  egrep "^(chr)?${c}[[:space:]]" $INPUT_BED > $SHARD_DIR/'GRCh37_chr'${c}.bed
  #liftover (note: fasta's cannot be zipped)
  java -jar -Xmx10g /Applications/bfx_tools/genomewarp-master/target/verilylifesciences-genomewarp-1.0.0-runnable.jar --lift_over_chain_path /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/hg19ToHg38.over.chain --raw_query_vcf GRCh37_chr${c}.vcf --raw_query_bed GRCh37_chr${c}.bed --ref_query_fasta /Volumes/Boron/reference\ files/ucsc.hg19.fasta --ref_target_fasta /Volumes/Boron/reference\ files/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna --work_dir /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG4_b37_vcfs_to_b38/cg/genomewarp_output/ --output_variants_file 'HG004_GRCh37_CHROM_'${c}'_vcfBeta-GS000037262-ASM_LIFTOVER_to_GRCh38.vcf' --output_regions_file 'HG004_GRCh37_CHROM_'${c}'_vcfBeta-GS000037262-ASM_LIFTOVER_to_GRCh38.bed'
  #remove extra contigs so only left with chr${c}
  egrep "^(#|(chr)?${c}[[:space:]])" 'HG004_GRCh37_CHROM_'${c}'_vcfBeta-GS000037262-ASM_LIFTOVER_to_GRCh38.vcf' > $SHARD_DIR/'HG004_'${c}'_convGRCh38_CG_vcfBeta-GS000037262-ASM.vcf'
  egrep "^(chr)?${c}[[:space:]]" 'HG004_GRCh37_CHROM_'${c}'_vcfBeta-GS000037262-ASM_LIFTOVER_to_GRCh38.bed' > $SHARD_DIR/'HG004_'${c}'_convGRCh38_CG_vcfBeta-GS000037262-ASM.bed'
  # Index and zip
  /Applications/bfx_tools/tabix-0.2.6/bgzip -f 'HG004_'${c}'_convGRCh38_CG_vcfBeta-GS000037262-ASM.vcf'
  /Applications/bfx_tools/tabix-0.2.6/tabix -f 'HG004_'${c}'_convGRCh38_CG_vcfBeta-GS000037262-ASM.vcf.gz'
done

##Had issue indexing chr6 vcf.gz (error: file out of order) need to sort first then index.  Sorted vcf.gz and tbi uploaded to DNAnexus.
#gunzip -c HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM.vcf.gz | sort -k2,2n | /Applications/bfx_tools/tabix-0.2.6/bgzip -c > HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_sorted.vcf.gz
#/Applications/bfx_tools/tabix-0.2.6/tabix -f HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_sorted.vcf.gz
#This sort did not work, integration failed. Justin said Sort -k2,2n sorts the second column numerically, and the header often only has one column but the order gets changed, need to remove header first, follow sort procedure below:

#Sort and index chr6 ONLY for CG.  Remove header prior to sort then add back in.
grep ^# HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM.vcf > HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_header.vcf

sed '/^#/ d' HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM.vcf > HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_NOheader.vcf

sort -k2,2n HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_NOheader.vcf > HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_NOheader_sorted.vcf

cat HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_header.vcf HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_NOheader_sorted.vcf > HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_sorted.vcf

/Applications/bfx_tools/tabix-0.2.6/bgzip -c HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_sorted.vcf > HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_sorted.vcf.gz

/Applications/bfx_tools/tabix-0.2.6/tabix -f HG004_6_convGRCh38_CG_vcfBeta-GS000037262-ASM_sorted.vcf.gz

#"sorted" removed from file name on DNAnexus


##### Ion ######

## 1) need to convert CHROM number for GRCh37 vcfs
# Run 12/01/16
# Parse out vcf header
grep ^# HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143.vcf > HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_header.vcf

#create no header file
sed '/^#/ d' HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143.vcf > HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_NOheader.vcf

# replace chrom # with chr# at start of every row in no header file
sed 's/^/chr/' HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_NOheader.vcf  > HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_NOheader_CHROMfixed.vcf

# replace chrom # with chr# at start of every row in bed
sed 's/^/chr/' HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_callable.bed > HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_callable_CHROMfixed.bed

# add header back in 
cat HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_header.vcf HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_NOheader_CHROMfixed.vcf > HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_CHROMfixed.vcf

## 2) liftover (note: fasta's cannot be zipped)
#Run 12/2/16
INPUT_VCF=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG4_b37_vcfs_to_b38/ion/HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_CHROMfixed.vcf
INPUT_BED=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG4_b37_vcfs_to_b38/ion/HG004_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24143_callable_CHROMfixed.bed
SHARD_DIR=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG4_b37_vcfs_to_b38/ion/
for c in `seq 1 22` X Y
do
  # Separate out each autosome and the X
  egrep "^(#|(chr)?${c}[[:space:]])" $INPUT_VCF > $SHARD_DIR/'GRCh37_chr'${c}.vcf
  egrep "^(chr)?${c}[[:space:]]" $INPUT_BED > $SHARD_DIR/'GRCh37_chr'${c}.bed
  #liftover (note: fasta's cannot be zipped)
  java -jar -Xmx10g /Applications/bfx_tools/genomewarp-master/target/verilylifesciences-genomewarp-1.0.0-runnable.jar --lift_over_chain_path /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/hg19ToHg38.over.chain --raw_query_vcf GRCh37_chr${c}.vcf --raw_query_bed GRCh37_chr${c}.bed --ref_query_fasta /Volumes/Boron/reference\ files/ucsc.hg19.fasta --ref_target_fasta /Volumes/Boron/reference\ files/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna --work_dir /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG4_b37_vcfs_to_b38/ion/genomewarp_output/ --output_variants_file 'HG004_GRCh37_CHROM_'${c}'_AmpliseqExome.20141120.NA24143_LIFTOVER_to_GRCh38.vcf' --output_regions_file 'HG004_GRCh37_CHROM_'${c}'_AmpliseqExome.20141120.NA24143_LIFTOVER_to_GRCh38.bed'
  #remove extra contigs so only left with chr${c}
  egrep "^(#|(chr)?${c}[[:space:]])" 'HG004_GRCh37_CHROM_'${c}'_AmpliseqExome.20141120.NA24143_LIFTOVER_to_GRCh38.vcf' > $SHARD_DIR/'HG004_'${c}'_convGRCh38_Ion_AmpliseqExome.20141120.NA24143.vcf'
  egrep "^(chr)?${c}[[:space:]]" 'HG004_GRCh37_CHROM_'${c}'_AmpliseqExome.20141120.NA24143_LIFTOVER_to_GRCh38.bed' > $SHARD_DIR/'HG004_'${c}'_convGRCh38_Ion_AmpliseqExome.20141120.NA24143.bed'
  # Index and zip
  /Applications/bfx_tools/tabix-0.2.6/bgzip -f 'HG004_'${c}'_convGRCh38_Ion_AmpliseqExome.20141120.NA24143.vcf'
  /Applications/bfx_tools/tabix-0.2.6/tabix -f 'HG004_'${c}'_convGRCh38_Ion_AmpliseqExome.20141120.NA24143.vcf.gz'
done


