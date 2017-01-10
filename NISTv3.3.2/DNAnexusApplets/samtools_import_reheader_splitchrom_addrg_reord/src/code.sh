#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x
java="java -Xmx1200m"

#
# Fetch mappings
#
aria2c "$urlbam" -o input.bam -x6 -s6 -j6 --check-certificate=false --file-allocation=none
aria2c "$urlbai" -o input.bai -x6 -s6 -j6 --check-certificate=false --file-allocation=none

opts="VALIDATION_STRINGENCY=$validation_stringency"
if [ "$advanced_options" != "" ]; then
  opts="$advanced_options"
fi

ls -l 
#
# Run samtools view to select reads from each chromosome, then reheader, add read groups, and output as new file
#    Note that chrM is excluded because sometimes chrM implies a different ref from MT and causes collisions downstream

samtools view -H input.bam | sed 's/chr//' > input_v37.sam
cat input_v37.sam

#samtools reheader input_v37.sam input.bam > input_v37.bam
#rm input.bam
#samtools view -H input_v37.bam
#samtools index input_v37.bam

samtools view -buh input.bam chr1 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam1.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam1.bam bam1.bai &
ls -l
samtools view -buh input.bam chr2 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam2.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam2.bam bam2.bai &
samtools view -buh input.bam chr3 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam3.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam3.bam bam3.bai &
samtools view -buh input.bam chr4 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam4.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam4.bam bam4.bai &
samtools view -buh input.bam chr5 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam5.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam5.bam bam5.bai &
samtools view -buh input.bam chr6 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam6.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam6.bam bam6.bai &
samtools view -buh input.bam chr7 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam7.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam7.bam bam7.bai &
samtools view -buh input.bam chr8 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam8.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam8.bam bam8.bai &
samtools view -buh input.bam chr9 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam9.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam9.bam bam9.bai &
samtools view -buh input.bam chr10 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam10.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam10.bam bam10.bai &
samtools view -buh input.bam chr11 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam11.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam11.bam bam11.bai &
samtools view -buh input.bam chr12 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam12.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam12.bam bam12.bai &
samtools view -buh input.bam chr13 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam13.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam13.bam bam13.bai &
samtools view -buh input.bam chr14 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam14.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam14.bam bam14.bai &
samtools view -buh input.bam chr15 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam15.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam15.bam bam15.bai &
samtools view -buh input.bam chr16 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam16.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam16.bam bam16.bai &
samtools view -buh input.bam chr17 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam17.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam17.bam bam17.bai &
samtools view -buh input.bam chr18 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam18.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam18.bam bam18.bai &
samtools view -buh input.bam chr19 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam19.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam19.bam bam19.bai &
samtools view -buh input.bam chr20 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam20.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam20.bam bam20.bai &
samtools view -buh input.bam chr21 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam21.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam21.bam bam21.bai &
samtools view -buh input.bam chr22 | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam22.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam22.bam bam22.bai &
samtools view -buh input.bam chrX | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bamX.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bamX.bam bamX.bai &
samtools view -buh input.bam chrY | samtools reheader input_v37.sam - | $java -jar /ReorderSam.jar I=/dev/stdin O=/dev/stdout S=true R=/human_g1k_v37_1to22XYonly.fasta | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bamY.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bamY.bam bamY.bai 

#
# Upload results
#
file_id1=`dx upload bam1.bam -o "$prefix.1.bam" --brief`
file_id1i=`dx upload bam1.bai -o "$prefix.1.bai" --brief`
dx-jobutil-add-output "bam1" "$file_id1"
dx-jobutil-add-output "bai1" "$file_id1i"
file_id2=`dx upload bam2.bam -o "$prefix.2.bam" --brief`
file_id2i=`dx upload bam2.bai -o "$prefix.2.bai" --brief`
dx-jobutil-add-output "bam2" "$file_id2"
dx-jobutil-add-output "bai2" "$file_id2i"
file_id3=`dx upload bam3.bam -o "$prefix.3.bam" --brief`
file_id3i=`dx upload bam3.bai -o "$prefix.3.bai" --brief`
dx-jobutil-add-output "bam3" "$file_id3"
dx-jobutil-add-output "bai3" "$file_id3i"
file_id4=`dx upload bam4.bam -o "$prefix.4.bam" --brief`
file_id4i=`dx upload bam4.bai -o "$prefix.4.bai" --brief`
dx-jobutil-add-output "bam4" "$file_id4"
dx-jobutil-add-output "bai4" "$file_id4i"
file_id5=`dx upload bam5.bam -o "$prefix.5.bam" --brief`
file_id5i=`dx upload bam5.bai -o "$prefix.5.bai" --brief`
dx-jobutil-add-output "bam5" "$file_id5"
dx-jobutil-add-output "bai5" "$file_id5i"
file_id6=`dx upload bam6.bam -o "$prefix.6.bam" --brief`
file_id6i=`dx upload bam6.bai -o "$prefix.6.bai" --brief`
dx-jobutil-add-output "bam6" "$file_id6"
dx-jobutil-add-output "bai6" "$file_id6i"
file_id7=`dx upload bam7.bam -o "$prefix.7.bam" --brief`
file_id7i=`dx upload bam7.bai -o "$prefix.7.bai" --brief`
dx-jobutil-add-output "bam7" "$file_id7"
dx-jobutil-add-output "bai7" "$file_id7i"
file_id8=`dx upload bam8.bam -o "$prefix.8.bam" --brief`
file_id8i=`dx upload bam8.bai -o "$prefix.8.bai" --brief`
dx-jobutil-add-output "bam8" "$file_id8"
dx-jobutil-add-output "bai8" "$file_id8i"
file_id9=`dx upload bam9.bam -o "$prefix.9.bam" --brief`
file_id9i=`dx upload bam9.bai -o "$prefix.9.bai" --brief`
dx-jobutil-add-output "bam9" "$file_id9"
dx-jobutil-add-output "bai9" "$file_id9i"
file_id10=`dx upload bam10.bam -o "$prefix.10.bam" --brief`
file_id10i=`dx upload bam10.bai -o "$prefix.10.bai" --brief`
dx-jobutil-add-output "bam10" "$file_id10"
dx-jobutil-add-output "bai10" "$file_id10i"
file_id11=`dx upload bam11.bam -o "$prefix.11.bam" --brief`
file_id11i=`dx upload bam11.bai -o "$prefix.11.bai" --brief`
dx-jobutil-add-output "bam11" "$file_id11"
dx-jobutil-add-output "bai11" "$file_id11i"
file_id12=`dx upload bam12.bam -o "$prefix.12.bam" --brief`
file_id12i=`dx upload bam12.bai -o "$prefix.12.bai" --brief`
dx-jobutil-add-output "bam12" "$file_id12"
dx-jobutil-add-output "bai12" "$file_id12i"
file_id13=`dx upload bam13.bam -o "$prefix.13.bam" --brief`
file_id13i=`dx upload bam13.bai -o "$prefix.13.bai" --brief`
dx-jobutil-add-output "bam13" "$file_id13"
dx-jobutil-add-output "bai13" "$file_id13i"
file_id14=`dx upload bam14.bam -o "$prefix.14.bam" --brief`
file_id14i=`dx upload bam14.bai -o "$prefix.14.bai" --brief`
dx-jobutil-add-output "bam14" "$file_id14"
dx-jobutil-add-output "bai14" "$file_id14i"
file_id15=`dx upload bam15.bam -o "$prefix.15.bam" --brief`
file_id15i=`dx upload bam15.bai -o "$prefix.15.bai" --brief`
dx-jobutil-add-output "bam15" "$file_id15"
dx-jobutil-add-output "bai15" "$file_id15i"
file_id16=`dx upload bam16.bam -o "$prefix.16.bam" --brief`
file_id16i=`dx upload bam16.bai -o "$prefix.16.bai" --brief`
dx-jobutil-add-output "bam16" "$file_id16"
dx-jobutil-add-output "bai16" "$file_id16i"
file_id17=`dx upload bam17.bam -o "$prefix.17.bam" --brief`
file_id17i=`dx upload bam17.bai -o "$prefix.17.bai" --brief`
dx-jobutil-add-output "bam17" "$file_id17"
dx-jobutil-add-output "bai17" "$file_id17i"
file_id18=`dx upload bam18.bam -o "$prefix.18.bam" --brief`
file_id18i=`dx upload bam18.bai -o "$prefix.18.bai" --brief`
dx-jobutil-add-output "bam18" "$file_id18"
dx-jobutil-add-output "bai18" "$file_id18i"
file_id19=`dx upload bam19.bam -o "$prefix.19.bam" --brief`
file_id19i=`dx upload bam19.bai -o "$prefix.19.bai" --brief`
dx-jobutil-add-output "bam19" "$file_id19"
dx-jobutil-add-output "bai19" "$file_id19i"
file_id20=`dx upload bam20.bam -o "$prefix.20.bam" --brief`
file_id20i=`dx upload bam20.bai -o "$prefix.20.bai" --brief`
dx-jobutil-add-output "bam20" "$file_id20"
dx-jobutil-add-output "bai20" "$file_id20i"
file_id21=`dx upload bam21.bam -o "$prefix.21.bam" --brief`
file_id21i=`dx upload bam21.bai -o "$prefix.21.bai" --brief`
dx-jobutil-add-output "bam21" "$file_id21"
dx-jobutil-add-output "bai21" "$file_id21i"
file_id22=`dx upload bam22.bam -o "$prefix.22.bam" --brief`
file_id22i=`dx upload bam22.bai -o "$prefix.22.bai" --brief`
dx-jobutil-add-output "bam22" "$file_id22"
dx-jobutil-add-output "bai22" "$file_id22i"
file_idX=`dx upload bamX.bam -o "$prefix.X.bam" --brief`
file_idXi=`dx upload bamX.bai -o "$prefix.X.bai" --brief`
dx-jobutil-add-output "bamX" "$file_idX"
dx-jobutil-add-output "baiX" "$file_idXi"
file_idY=`dx upload bamY.bam -o "$prefix.Y.bam" --brief`
file_idYi=`dx upload bamY.bai -o "$prefix.Y.bai" --brief`
dx-jobutil-add-output "bamY" "$file_idY"
dx-jobutil-add-output "baiY" "$file_idYi"
