#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs 
#
#dx-download-all-inputs --parallel
dx download "$vcfhighconfgz" -o vcfhighconfgz.vcf.gz
dx download "$vcftestgz" -o vcftestgz.vcf.gz
dx download "$vcfhighconftbi" -o vcfhighconfgz.vcf.gz.tbi
dx download "$vcftesttbi" -o vcftestgz.vcf.gz.tbi
dx download "$bedhighconf" -o vcfhighconf.bed
dx download "$bedtest" -o vcftest.bed
dx download "$svbed" -o svs.bed
dx download "$rtgsdf" -o rtgsdf.tar.gz
#
# Stream and unpack genome
#
# => RTG reference sdf files
  tar zxvf rtgsdf.tar.gz




##compare and stratify variant calls
#gunzip -c vcftestgz.vcf.gz | awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || chrom == "all" || $1 ~ /^#/) print }'


rtg vcffilter -i vcfhighconfgz.vcf.gz --include-bed vcfhighconf.bed -o vcfhighconfgz_inhighconfbed.vcf.gz
rtg vcfstats vcfhighconfgz_inhighconfbed.vcf.gz

awk '{sum+=$3-$2} END {print sum}' vcfhighconf.bed

rtg vcfeval --ref-overlap -b vcfhighconfgz.vcf.gz -c vcftestgz.vcf.gz -o comp -t rtgsdf

rtg vcfstats comp/tp.vcf.gz 
rtg vcfstats comp/fp.vcf.gz 
rtg vcfstats comp/fn.vcf.gz 

rtg vcffilter -i comp/tp.vcf.gz --include-bed vcfhighconf.bed -o comp/tp_inhighconfbed.vcf.gz

rtg vcffilter -i comp/tp.vcf.gz --include-bed vcftest.bed -o comp/tp_intestbed.vcf.gz

rtg vcffilter -i comp/tp_inhighconfbed.vcf.gz --include-bed vcftest.bed -o comp/tp_inhighconfbed_intestbed.vcf.gz

rtg vcffilter -i comp/tp.vcf.gz --exclude-bed vcfhighconf.bed -o comp/tp_notinhighconfbed.vcf.gz

rtg vcffilter -i comp/tp_notinhighconfbed.vcf.gz --exclude-bed vcftest.bed -o comp/tp_notinhighconfbed_notintestbed.vcf.gz


rtg vcffilter -i comp/tp.vcf.gz --include-bed svs.bed -o comp/tp_insvbed.vcf.gz


rtg vcffilter -i comp/fp.vcf.gz --include-bed vcfhighconf.bed -o comp/fp_inhighconfbed.vcf.gz

rtg vcffilter -i comp/fp.vcf.gz --include-bed vcftest.bed -o comp/fp_intestbed.vcf.gz

rtg vcffilter -i comp/fp_inhighconfbed.vcf.gz --include-bed vcftest.bed -o comp/fp_inhighconfbed_intestbed.vcf.gz

rtg vcffilter -i comp/fp.vcf.gz --exclude-bed vcfhighconf.bed -o comp/fp_notinhighconfbed.vcf.gz

rtg vcffilter -i comp/fp_notinhighconfbed.vcf.gz --exclude-bed vcftest.bed -o comp/fp_notinhighconfbed_notintestbed.vcf.gz

rtg vcffilter -i comp/fp.vcf.gz --include-bed svs.bed -o comp/fp_insvbed.vcf.gz


rtg vcffilter -i comp/fn.vcf.gz --include-bed vcfhighconf.bed -o comp/fn_inhighconfbed.vcf.gz

rtg vcffilter -i comp/fn.vcf.gz --include-bed vcftest.bed -o comp/fn_intestbed.vcf.gz

rtg vcffilter -i comp/fn_inhighconfbed.vcf.gz --include-bed vcftest.bed -o comp/fn_inhighconfbed_intestbed.vcf.gz

rtg vcffilter -i comp/fn.vcf.gz --exclude-bed vcfhighconf.bed -o comp/fn_notinhighconfbed.vcf.gz

rtg vcffilter -i comp/fn_notinhighconfbed.vcf.gz --exclude-bed vcftest.bed -o comp/fn_notinhighconfbed_notintestbed.vcf.gz

rtg vcffilter -i comp/fn.vcf.gz --include-bed svs.bed -o comp/fn_insvbed.vcf.gz

rtg vcffilter -i vcftestgz.vcf.gz --exclude-bed vcfhighconf.bed -o comp/vcftest_notinhighconfbed.vcf.gz

##also look at differences without requiring GT match
rtg vcfeval --ref-overlap --squash-ploidy -b vcfhighconfgz.vcf.gz -c vcftestgz.vcf.gz -o comp/compnoGT -t rtgsdf
rtg vcffilter -i comp/compnoGT/tp.vcf.gz --include-bed vcfhighconf.bed -o comp/compnoGT/tp_inhighconfbed.vcf.gz

rtg vcffilter -i comp/compnoGT/tp_inhighconfbed.vcf.gz --include-bed vcftest.bed -o comp/compnoGT/tp_inhighconfbed_intestbed.vcf.gz

rtg vcffilter -i comp/compnoGT/fp.vcf.gz --include-bed vcfhighconf.bed -o comp/compnoGT/fp_inhighconfbed.vcf.gz

rtg vcffilter -i comp/compnoGT/fp_inhighconfbed.vcf.gz --include-bed vcftest.bed -o comp/compnoGT/fp_inhighconfbed_intestbed.vcf.gz

rtg vcffilter -i comp/compnoGT/fn.vcf.gz --include-bed vcfhighconf.bed -o comp/compnoGT/fn_inhighconfbed.vcf.gz

rtg vcffilter -i comp/compnoGT/fn_inhighconfbed.vcf.gz --include-bed vcftest.bed -o comp/compnoGT/fn_inhighconfbed_intestbed.vcf.gz


tar -zcvf comparison.tar.gz comp


#
# Upload results
#
file_idcomp=`dx upload comparison.tar.gz -o "$prefix".tar.gz --brief`
dx-jobutil-add-output "compfolder" "$file_idcomp"

file_idtpinbedsvcf=`dx upload comp/tp_inhighconfbed_intestbed.vcf.gz -o "$prefix"_tpinbeds.vcf.gz --brief`
dx-jobutil-add-output "tpinbedsvcf" "$file_idtpinbedsvcf"

file_idtpinbedsvcftbi=`dx upload comp/tp_inhighconfbed_intestbed.vcf.gz.tbi -o "$prefix"_tpinbeds.vcf.gz.tbi --brief`
dx-jobutil-add-output "tpinbedsvcftbi" "$file_idtpinbedsvcftbi"

file_idfpinbedsvcf=`dx upload comp/fp_inhighconfbed_intestbed.vcf.gz -o "$prefix"_fpinbeds.vcf.gz --brief`
dx-jobutil-add-output "fpinbedsvcf" "$file_idfpinbedsvcf"

file_idfpinbedsvcftbi=`dx upload comp/fp_inhighconfbed_intestbed.vcf.gz.tbi -o "$prefix"_fpinbeds.vcf.gz.tbi --brief`
dx-jobutil-add-output "fpinbedsvcftbi" "$file_idfpinbedsvcftbi"

file_idfninbedsvcf=`dx upload comp/fn_inhighconfbed_intestbed.vcf.gz -o "$prefix"_fninbeds.vcf.gz --brief`
dx-jobutil-add-output "fninbedsvcf" "$file_idfninbedsvcf"

file_idfninbedsvcftbi=`dx upload comp/fn_inhighconfbed_intestbed.vcf.gz.tbi -o "$prefix"_fninbeds.vcf.gz.tbi --brief`
dx-jobutil-add-output "fninbedsvcftbi" "$file_idfninbedsvcftbi"

file_idfninhighconfbedvcf=`dx upload comp/fn_inhighconfbed.vcf.gz -o "$prefix"_fninhighconfbed.vcf.gz --brief`
dx-jobutil-add-output "fninhighconfbedvcf" "$file_idfninhighconfbedvcf"

file_idfninhighconfbedvcftbi=`dx upload comp/fn_inhighconfbed.vcf.gz.tbi -o "$prefix"_fninhighconfbed.vcf.gz.tbi --brief`
dx-jobutil-add-output "fninhighconfbedvcftbi" "$file_idfninhighconfbedvcftbi"

