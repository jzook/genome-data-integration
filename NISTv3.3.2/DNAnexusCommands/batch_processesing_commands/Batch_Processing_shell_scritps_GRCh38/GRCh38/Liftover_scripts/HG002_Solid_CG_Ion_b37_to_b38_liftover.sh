All complete

#HG002 Liftover b37 to b38

##### SOLID ###### 

## 1) need to convert CHROM number for GRCh37 highconf files
#Run 11/28/16
# Parse out vcf header
grep ^# HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC.vcf > HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC_header.vcf

#create no header file
sed '/^#/ d' HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC.vcf > HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC_NOheader.vcf

# replace chrom # with chr# at start of every row in no header file
sed 's/^/chr/' HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC_NOheader.vcf  > HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC_NOheader_CHROMfixed.vcf

#for solid need to pull out "callable" from bed
grep 'CALLABLE' HG002_GRCh37_CHROM1-Y_solid5500_SE75bp_callableloci.bed > HG002_GRCh37_CHROM1-Y_solid5500_SE75bp_callable.bed

# replace chrom # with chr# at start of every row in bed
sed 's/^/chr/' HG002_GRCh37_CHROM1-Y_solid5500_SE75bp_callable.bed > HG002_GRCh37_CHROM1-Y_solid5500_SE75bp_callable_CHROMfixed.bed

# add header back in 
cat HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC_header.vcf HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC_NOheader_CHROMfixed.vcf > HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC_CHROMfixed.vcf

# 2) Process genome warp
#Run 12/01/16
INPUT_VCF=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/solid/HG002_GRCh37_CHROM1-Y_hg19_solid5500_SE75bp_GATKHC_CHROMfixed.vcf
INPUT_BED=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/solid/HG002_GRCh37_CHROM1-Y_solid5500_SE75bp_callable_CHROMfixed.bed
SHARD_DIR=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/solid/
for c in `seq 1 22` X Y
do
	egrep "^(#|(chr)?${c}[[:space:]])" $INPUT_VCF > $SHARD_DIR/'GRCh37_chr'${c}.vcf
	egrep "^(chr)?${c}[[:space:]]" $INPUT_BED > $SHARD_DIR/'GRCh37_chr'${c}.bed
	java -jar -Xmx10g /Applications/bfx_tools/genomewarp-master/target/verilylifesciences-genomewarp-1.0.0-runnable.jar --lift_over_chain_path /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/hg19ToHg38.over.chain --raw_query_vcf GRCh37_chr${c}.vcf --raw_query_bed GRCh37_chr${c}.bed --ref_query_fasta /Volumes/Boron/reference\ files/ucsc.hg19.fasta --ref_target_fasta /Volumes/Boron/reference\ files/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna --work_dir /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/solid/genomewarp_output/ --output_variants_file 'HG002_GRCh37_CHROM_'${c}'_hg19_solid5500_SE75bp_GATKHC_LIFTOVER_to_GRCh38.vcf' --output_regions_file 'HG002_GRCh37_CHROM_'${c}'_hg19_solid5500_SE75bp_GATKHC_LIFTOVER_to_GRCh38.bed'
	egrep "^(#|(chr)?${c}[[:space:]])" 'HG002_GRCh37_CHROM_'${c}'_hg19_solid5500_SE75bp_GATKHC_LIFTOVER_to_GRCh38.vcf' > $SHARD_DIR/'HG002_'${c}'_convGRCh38_solid5500_SE75bp_GATKHC.vcf'
	egrep "^(chr)?${c}[[:space:]]" 'HG002_GRCh37_CHROM_'${c}'_hg19_solid5500_SE75bp_GATKHC_LIFTOVER_to_GRCh38.bed' > $SHARD_DIR/'HG002_'${c}'_convGRCh38_solid5500_SE75bp_GATKHC.bed'
	/Applications/bfx_tools/tabix-0.2.6/bgzip -f 'HG002_'${c}'_convGRCh38_solid5500_SE75bp_GATKHC.vcf'
	/Applications/bfx_tools/tabix-0.2.6/tabix -f 'HG002_'${c}'_convGRCh38_solid5500_SE75bp_GATKHC.vcf.gz'
done


##### Complete Genomics ######

## 1) need to convert CHROM number for GRCh37 highconf files
#Run 11/28/16
# Parse out vcf header
grep ^# HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM.vcf > HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_header.vcf

#create no header file
sed '/^#/ d' HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM.vcf > HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_NOheader.vcf

# replace chrom # with chr# at start of every row in no header file
sed 's/^/chr/' HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_NOheader.vcf  > HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_NOheader_CHROMfixed.vcf

# replace chrom # with chr# at start of every row in bed
sed 's/^/chr/' HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_callable.bed > HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_callable_CHROMfixed.bed

# add header back in 
cat HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_header.vcf HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_NOheader_CHROMfixed.vcf > HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_CHROMfixed.vcf

# 2) Process genome warp
#Run 12/01/16
INPUT_VCF=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/cg/HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_CHROMfixed.vcf
INPUT_BED=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/cg/HG002_GRCh37_CHROM1-Y_vcfBeta-GS000037263-ASM_callable_CHROMfixed.bed
SHARD_DIR=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/cg/
for c in `seq 1 22` X Y
do
  # Separate out each autosome and the X
  egrep "^(#|(chr)?${c}[[:space:]])" $INPUT_VCF > $SHARD_DIR/'GRCh37_chr'${c}.vcf
  egrep "^(chr)?${c}[[:space:]]" $INPUT_BED > $SHARD_DIR/'GRCh37_chr'${c}.bed
  #liftover (note: fasta's cannot be zipped)
  java -jar -Xmx10g /Applications/bfx_tools/genomewarp-master/target/verilylifesciences-genomewarp-1.0.0-runnable.jar --lift_over_chain_path /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/hg19ToHg38.over.chain --raw_query_vcf GRCh37_chr${c}.vcf --raw_query_bed GRCh37_chr${c}.bed --ref_query_fasta /Volumes/Boron/reference\ files/ucsc.hg19.fasta --ref_target_fasta /Volumes/Boron/reference\ files/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna --work_dir /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/cg/genomewarp_output/ --output_variants_file 'HG002_GRCh37_CHROM_'${c}'_vcfBeta-GS000037263-ASM_LIFTOVER_to_GRCh38.vcf' --output_regions_file 'HG002_GRCh37_CHROM_'${c}'_vcfBeta-GS000037263-ASM_LIFTOVER_to_GRCh38.bed'
  #remove extra contigs so only left with chr${c}
  egrep "^(#|(chr)?${c}[[:space:]])" 'HG002_GRCh37_CHROM_'${c}'_vcfBeta-GS000037263-ASM_LIFTOVER_to_GRCh38.vcf' > $SHARD_DIR/'HG002_'${c}'_convGRCh38_CG_vcfBeta-GS000037263-ASM.vcf'
  egrep "^(chr)?${c}[[:space:]]" 'HG002_GRCh37_CHROM_'${c}'_vcfBeta-GS000037263-ASM_LIFTOVER_to_GRCh38.bed' > $SHARD_DIR/'HG002_'${c}'_convGRCh38_CG_vcfBeta-GS000037263-ASM.bed'
  # Index and zip
  /Applications/bfx_tools/tabix-0.2.6/bgzip -f 'HG002_'${c}'_convGRCh38_CG_vcfBeta-GS000037263-ASM.vcf'
  /Applications/bfx_tools/tabix-0.2.6/tabix -f 'HG002_'${c}'_convGRCh38_CG_vcfBeta-GS000037263-ASM.vcf.gz'
done

##Had issue indexing chr6 vcf.gz (error: file out of order) need to sort first then index.  Sorted vcf.gz and tbi uploaded to DNAnexus.
#gunzip -c HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM.vcf.gz | sort -k2,2n | /Applications/bfx_tools/tabix-0.2.6/bgzip -c > HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_sorted.vcf.gz
#/Applications/bfx_tools/tabix-0.2.6/tabix -f HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_sorted.vcf.gz
#This sort did not work, integration failed. Justin said Sort -k2,2n sorts the second column numerically, and the header often only has one column but the order gets changed, need to remove header first, follow sort procedure below:

#Sort and index chr6 ONLY for CG.  Remove header prior to sort then add back in. sort done 12/22/16
grep ^# HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM.vcf > HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_header.vcf

sed '/^#/ d' HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM.vcf > HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_NOheader.vcf

sort -k2,2n HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_NOheader.vcf > HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_NOheader_sorted.vcf

cat HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_header.vcf HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_NOheader_sorted.vcf > HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_sorted.vcf

/Applications/bfx_tools/tabix-0.2.6/bgzip -c HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_sorted.vcf > HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_sorted.vcf.gz

/Applications/bfx_tools/tabix-0.2.6/tabix -f HG002_6_convGRCh38_CG_vcfBeta-GS000037263-ASM_sorted.vcf.gz

#"sorted" removed from file name on DNAnexus

##### Ion ######

## 1) need to convert CHROM number for GRCh37 vcfs
#Run 11/28/16
# Parse out vcf header
grep ^# HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385.vcf > HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_header.vcf

#create no header file
sed '/^#/ d' HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385.vcf > HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_NOheader.vcf

# replace chrom # with chr# at start of every row in no header file
sed 's/^/chr/' HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_NOheader.vcf  > HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_NOheader_CHROMfixed.vcf

# replace chrom # with chr# at start of every row in bed
sed 's/^/chr/' HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_callable.bed > HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_callable_CHROMfixed.bed

# add header back in 
cat HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_header.vcf HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_NOheader_CHROMfixed.vcf > HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_CHROMfixed.vcf

# 2) Process genome warp
#Run 12/01/16
INPUT_VCF=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/ion/HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_CHROMfixed.vcf
INPUT_BED=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/ion/HG002_GRCh37_CHROM1-X_AmpliseqExome.20141120.NA24385_callable_CHROMfixed.bed
SHARD_DIR=/Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/ion/
for c in `seq 1 22` X Y
do
  # Separate out each autosome and the X
  egrep "^(#|(chr)?${c}[[:space:]])" $INPUT_VCF > $SHARD_DIR/'GRCh37_chr'${c}.vcf
  egrep "^(chr)?${c}[[:space:]]" $INPUT_BED > $SHARD_DIR/'GRCh37_chr'${c}.bed
  #liftover (note: fasta's cannot be zipped)
  java -jar -Xmx10g /Applications/bfx_tools/genomewarp-master/target/verilylifesciences-genomewarp-1.0.0-runnable.jar --lift_over_chain_path /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/hg19ToHg38.over.chain --raw_query_vcf GRCh37_chr${c}.vcf --raw_query_bed GRCh37_chr${c}.bed --ref_query_fasta /Volumes/Boron/reference\ files/ucsc.hg19.fasta --ref_target_fasta /Volumes/Boron/reference\ files/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna --work_dir /Users/jmcdani/Documents/GiaB/GiAB_informatics_work/Liftover-GRCh37toGRCh38/HG2_b37_vcfs_to_b38/ion/genomewarp_output/ --output_variants_file 'HG002_GRCh37_CHROM_'${c}'_AmpliseqExome.20141120.NA24385_LIFTOVER_to_GRCh38.vcf' --output_regions_file 'HG002_GRCh37_CHROM_'${c}'_AmpliseqExome.20141120.NA24385_LIFTOVER_to_GRCh38.bed'
  #remove extra contigs so only left with chr${c}
  egrep "^(#|(chr)?${c}[[:space:]])" 'HG002_GRCh37_CHROM_'${c}'_AmpliseqExome.20141120.NA24385_LIFTOVER_to_GRCh38.vcf' > $SHARD_DIR/'HG002_'${c}'_convGRCh38_Ion_AmpliseqExome.20141120.NA24385.vcf'
  egrep "^(chr)?${c}[[:space:]]" 'HG002_GRCh37_CHROM_'${c}'_AmpliseqExome.20141120.NA24385_LIFTOVER_to_GRCh38.bed' > $SHARD_DIR/'HG002_'${c}'_convGRCh38_Ion_AmpliseqExome.20141120.NA24385.bed'
  # Index and zip
  /Applications/bfx_tools/tabix-0.2.6/bgzip -f 'HG002_'${c}'_convGRCh38_Ion_AmpliseqExome.20141120.NA24385.vcf'
  /Applications/bfx_tools/tabix-0.2.6/tabix -f 'HG002_'${c}'_convGRCh38_Ion_AmpliseqExome.20141120.NA24385.vcf.gz'
done
