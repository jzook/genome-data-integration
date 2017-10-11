#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch inputs (~/in/vcfs/0/*, ~/in/vcfs/1/*, etc...)
#
dx-download-all-inputs --parallel
dx download "$callsettable" -o callsettable.txt
dx download "$refn" -o refn.bed
#dx download "$svbed" -o svs.bed

##Move all files from file arrays to the same data folder
mkdir data
mv ~/in/vcfs/*/* data
mv ~/in/beds/*/* data
mv ~/in/annotations/*/* data
mv ~/in/filtbeds/*/* .
mv ~/in/ref/* .
mv ~/in/rtgsdf/* .
#unzip ref and rtgsdf files
for v in *.tar.gz; do
  tar zxvf "$v"
done
mv /usr/bin/*.pl .

#make genome and bed files for bedtools
grep '^@SQ' genome.dict  | awk 'NR < 26' | sed 's/.*SN:\(.*\)\tLN:\(.*\)\tM5:.*/\1\t\2/' > human.genome
grep '^@SQ' genome.dict  | awk 'NR < 26' | sed 's/.*SN:\(.*\)\tLN:\(.*\)\tM5:.*/\1\t1\t\2/' > human.genome.bed

##Process vcf and bed files (adding sample names, removing chr, removing homref sites)
perl preprocess_combine_vcfs.pl -cstable callsettable.txt -p data/

##Create environment variable with names of processed vcf files and make union vcf using vcflib vcfcombine
##The union vcf has a column corresponding to the genotype call from each vcf
vars=""
for v in data/*_nohomref_samplename.vcf.gz; do
  vars="$vars $v"
done
vcfcombine $vars > union.1.vcf

##Annotate union vcf with each of the callable bed files in the "callable" INFO field
i=1
for v in data/*_callable_processed.bed; do
  ((j=i+1))
  vcfannotate -b "$v" -k callable union."$i".vcf > union."$j".vcf
  rm union."$i".vcf
  awk '{ sum+=$3; sum-=$2 } END { print sum }' "$v"
  ((i=i+1))
done


##Annotate union vcf with difficult region bed files in the "difficultregion" INFO field (currently only used for information purposes)
for v in *_diffbedannotcol.bed; do
  ((j=i+1))
  vcfannotate -b "$v" -k difficultregion union."$i".vcf > union."$j".vcf
  rm union."$i".vcf
  awk '{ sum+=$3; sum-=$2 } END { print sum }' "$v"
  ((i=i+1))
done
mv union."$j".vcf union_callableannotate.vcf

##Run first pass of integration to find training variants found in 2 technologies
perl VcfClassifyUsingFilters_v3.3.pl union_callableannotate.vcf callsettable.txt union_callableannotate

#find number of lines in vcfs for debugging
#grep -v ^## union_callableannotate_ClassifyUsingFilters_2platforms.vcf | head | awk '{print $0}'
wc -l union_callableannotate_ClassifyUsingFilters_2platforms.vcf
wc -l union_callableannotate_ClassifyUsingFilters_arbitrated.vcf

#bgzip and index integrated vcfs
bgzip union_callableannotate_ClassifyUsingFilters_2platforms.vcf
tabix -p vcf union_callableannotate_ClassifyUsingFilters_2platforms.vcf.gz
bgzip union_callableannotate_ClassifyUsingFilters_allcalls.vcf
tabix -p vcf union_callableannotate_ClassifyUsingFilters_allcalls.vcf.gz

##intersect calls found in 2 technologies with calls from each callset to give training variants to one class filter (not needed for vgraph method)
for v in data/*_nohomref_samplename.vcf.gz; do
  vcfintersect --intersect-vcf union_callableannotate_ClassifyUsingFilters_2platforms.vcf.gz -r genome.fa "$v" > "$v"_indivintersect2platforms.vcf
  wc -l "$v"_indivintersect2platforms.vcf
done

##Run perl script that runs another perl script that performs one-class filtering for each callset and generates filtering bed files for each callset (also generates filtered vcf files for vgraph method)
perl RunOneClassFilter.pl -cstable callsettable.txt -p data/

##Annotate union vcf with each of the filtering bed files in the "oneclassfilt" INFO field
cp union_callableannotate.vcf union_callableannotate.1.vcf
i=1
for v in data/*_samplename_oneclassfilter_filtered.bed; do
  ((j=i+1))
  vcfannotate -b "$v" -k oneclassfilt union_callableannotate."$i".vcf > union_callableannotate."$j".vcf
  rm union_callableannotate."$i".vcf
  ((i=i+1))
done
mv union_callableannotate."$j".vcf union_callableannotate_filterannotate.vcf

##Run final pass of integration to find high-confidence variants using not-callable bed and filtering bed for arbitartion
perl VcfClassifyUsingFilters_v3.3.pl union_callableannotate_filterannotate.vcf callsettable.txt union_callableannotate_filterannotate

bgzip union_callableannotate_filterannotate.vcf
tabix -p vcf union_callableannotate_filterannotate.vcf.gz

bgzip union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf
tabix -p vcf union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz

bgzip union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf
tabix -p vcf union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf.gz

bgzip union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf
tabix -p vcf union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf.gz

##Calculate some stats about integrated variants
rtg vcfstats union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz

rtg vcfstats union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf.gz

rtg vcfstats union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf.gz

zgrep GQlessthan70 union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l
zgrep allfilteredanddisagree union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l
zgrep allfilteredbutagree union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l
zgrep discordantunfiltered union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l
zgrep discordanthet union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l
zgrep questionableindel union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l
zgrep cgonly union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l
zgrep alleleimbalance union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l

##Create high-confidence bed file
##first, find how many datasets cover each region
awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || chrom == "allbutY" || chrom == "all" || $1 ~ /^#/) print }' human.genome > human.chrom.genome
awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || chrom == "allbutY" || chrom == "all" || $1 ~ /^#/) print }' human.genome.bed > human.chrom.genome.bed

vars=""
for v in data/*_callable_processed.bed; do
  count=`wc -l $v | awk '{print $1}'`
  if [ ${count} -gt 1 ]; then
  	vars="$vars ${v}"
  fi
done

multiIntersectBed -i $vars | awk -v chrom="$chrom" 'BEGIN {FS = OFS = "\t"} {if($1 == chrom || chrom == "allbutY" || chrom == "all" || $1 ~ /^#/) print }' > callablemultinter.bed

##Create bed files with regions that are callable in at least 1 or at least 2 datasets
awk '{if($4>1) print $0}' callablemultinter.bed | mergeBed -i stdin > callablemultinter_gt1.bed
awk '{if($4>0) print $0}' callablemultinter.bed | mergeBed -i stdin > callablemultinter_gt0.bed
awk '{ sum+=$3; sum-=$2 } END { print sum }' callablemultinter_gt1.bed
awk '{ sum+=$3; sum-=$2 } END { print sum }' callablemultinter_gt0.bed

##Create bed file with every potential variant that is not high-confidence +-50bp
gunzip -c union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | awk '{OFS="\t"; if (!($7 ~ /^PASS/ || $7 ~ /^\./ || $1 ~ /^#/)) {print $1,$2-1,$2,$4"/"$5}}' | subtractBed -a human.chrom.genome.bed -b stdin | complementBed -i stdin -g human.chrom.genome | slopBed -i stdin -g human.chrom.genome -b 50 > union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed

rtg vcfeval -b union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf.gz -c union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf.gz -o selfcompare -t rtgsdf

zgrep -v ^# selfcompare/fp.vcf.gz | awk 'BEGIN {FS = OFS = "\t"} {print $1,$2-50,$2+length($4)+50}' | sed 's/^X/23/;s/^Y/24/;s/^M/25/' | sort -k1,1n -k2,2n | sed 's/^23/X/;s/^24/Y/;s/^25/MT/' | mergeBed -i stdin -d 50 > excludeoverlappingvars.bed

##Subtract regions around filtered variant, and overlapping variants; before v3.3, also excluded tandem repeats >200bp in length, homopolymers>10bp in length, segmental duplications>10kb in size, segmental duplications from Eichler group, and potential SVs to generate final high-confidence regions
subtractBed -a callablemultinter_gt0.bed -b union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed | subtractBed -a stdin -b excludeoverlappingvars.bed | subtractBed -a stdin -b refn.bed  > callablemultinter_gt0_nofilt_nooverlapvar.bed
subtractBed -a callablemultinter_gt1.bed -b union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed | subtractBed -a stdin -b excludeoverlappingvars.bed | subtractBed -a stdin -b refn.bed  > callablemultinter_gt1_nofilt_nooverlapvar.bed
awk '{ sum+=$3; sum-=$2 } END { print sum }' callablemultinter_gt1_nofilt_nooverlapvar.bed
awk '{ sum+=$3; sum-=$2 } END { print sum }' callablemultinter_gt0_nofilt_nooverlapvar.bed

##Calculate some stats about integrated variants

rtg vcffilter -i selfcompare/tp.vcf.gz --include-bed callablemultinter_gt0_nofilt_nooverlapvar.bed -o selfcompare/tp_inhighconfbed.vcf.gz

rtg vcfstats selfcompare/tp_inhighconfbed.vcf.gz

rtg vcffilter -i selfcompare/tp.vcf.gz --exclude-bed callablemultinter_gt0.bed -o selfcompare/tp_notinuniongt0bed.vcf.gz

rtg vcffilter -i selfcompare/tp.vcf.gz --include-bed union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed -o selfcompare/tp_infiltslop50bed.vcf.gz

rtg vcffilter -i selfcompare/tp.vcf.gz --include-bed excludeoverlappingvars.bed -o selfcompare/tp_inexcludeoverlappingvarsbed.vcf.gz

tar -zcf indivfilteredbeds.tar.gz data/*_callable_processed.bed

#
# Upload results
#
file_idvcfann=`dx upload union_callableannotate_filterannotate.vcf.gz -o "$prefix"_annotated.vcf.gz --brief`
dx-jobutil-add-output "vcfanngz" "$file_idvcfann"
file_idvcfanntbi=`dx upload union_callableannotate_filterannotate.vcf.gz.tbi -o "$prefix"_annotated.vcf.gz.tbi --brief`
dx-jobutil-add-output "vcfanntbi" "$file_idvcfanntbi"

file_idvcfall=`dx upload union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz -o "$prefix"_all.vcf.gz --brief`
dx-jobutil-add-output "vcfallgz" "$file_idvcfall"
file_idvcfalltbi=`dx upload union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz.tbi -o "$prefix"_all.vcf.gz.tbi --brief`
dx-jobutil-add-output "vcfalltbi" "$file_idvcfalltbi"

file_idvcfhighconf=`dx upload selfcompare/tp.vcf.gz -o "$prefix"_highconf.vcf.gz --brief`
dx-jobutil-add-output "vcfhighconfgz" "$file_idvcfhighconf"
file_idvcfhighconftbi=`dx upload selfcompare/tp.vcf.gz.tbi -o "$prefix"_highconf.vcf.gz.tbi --brief`
dx-jobutil-add-output "vcfhighconftbi" "$file_idvcfhighconftbi"

file_idvcfhighconf2tech=`dx upload union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf.gz -o "$prefix"_highconf_2tech.vcf.gz --brief`
dx-jobutil-add-output "vcfhighconf2techgz" "$file_idvcfhighconf2tech"
file_idvcfhighconf2techtbi=`dx upload union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf.gz.tbi -o "$prefix"_highconf_2tech.vcf.gz.tbi --brief`
dx-jobutil-add-output "vcfhighconf2techtbi" "$file_idvcfhighconf2techtbi"

file_idvcfallnofilt=`dx upload union_callableannotate_ClassifyUsingFilters_allcalls.vcf.gz -o "$prefix"_nofilt_all.vcf.gz --brief`
dx-jobutil-add-output "vcfallnofiltgz" "$file_idvcfallnofilt"
file_idvcfallnofilttbi=`dx upload union_callableannotate_ClassifyUsingFilters_allcalls.vcf.gz.tbi -o "$prefix"_nofilt_all.vcf.gz.tbi --brief`
dx-jobutil-add-output "vcfallnofilttbi" "$file_idvcfallnofilttbi"

file_idbed1=`dx upload callablemultinter_gt0_nofilt_nooverlapvar.bed -o "$prefix"_highconf.bed --brief`
dx-jobutil-add-output "bed1" "$file_idbed1"
file_idbed2=`dx upload callablemultinter_gt1_nofilt_nooverlapvar.bed -o "$prefix"_highconf_gt1.bed --brief`
dx-jobutil-add-output "bed2" "$file_idbed2"

file_idbedcallablegt0=`dx upload callablemultinter_gt0.bed -o "$prefix"_callablemultinter_gt0.bed --brief`
dx-jobutil-add-output "bedcallablegt0" "$file_idbedcallablegt0"
file_idbedcallableall=`dx upload callablemultinter.bed -o "$prefix"_callablemultinter_all.bed --brief`
dx-jobutil-add-output "bedcallableall" "$file_idbedcallableall"

file_idbedfilteredsites=`dx upload union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed -o "$prefix"_filteredsites.bed --brief`
dx-jobutil-add-output "bedfilteredsites" "$file_idbedfilteredsites"

file_idindivbedfilteredsites=`dx upload indivfilteredbeds.tar.gz -o "$prefix"_indivfilteredbeds.tar.gz --brief`
dx-jobutil-add-output "bedindivfilteredsites" "$file_idindivbedfilteredsites"
