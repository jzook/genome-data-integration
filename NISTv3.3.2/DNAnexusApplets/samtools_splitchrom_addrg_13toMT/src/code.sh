#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x
java="java -Xmx1200m"

#
# Fetch mappings
#
#dx-download-all-inputs --parallel
dx download "$sorted_bam" -o input.bam
dx download "$index_bai" -o input.bai
#mv ~/in/index_bai/* ~/in/sorted_bam/

opts="VALIDATION_STRINGENCY=$validation_stringency"
if [ "$advanced_options" != "" ]; then
  opts="$advanced_options"
fi

ls -l 
#
# Run samtools view to select reads from each chromosome and output as new file
#
samtools view -buh input.bam 13 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam13.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam13.bam bam13.bai &
samtools view -buh input.bam 14 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam14.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam14.bam bam14.bai &
samtools view -buh input.bam 15 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam15.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam15.bam bam15.bai &
samtools view -buh input.bam 16 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam16.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam16.bam bam16.bai &
samtools view -buh input.bam 17 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam17.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam17.bam bam17.bai &
samtools view -buh input.bam 18 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam18.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam18.bam bam18.bai &
samtools view -buh input.bam 19 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam19.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam19.bam bam19.bai &
samtools view -buh input.bam 20 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam20.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam20.bam bam20.bai &
samtools view -buh input.bam 21 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam21.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam21.bam bam21.bai &
samtools view -buh input.bam 22 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam22.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam22.bam bam22.bai &
samtools view -buh input.bam X | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bamX.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bamX.bam bamX.bai &
samtools view -buh input.bam Y | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bamY.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bamY.bam bamY.bai 
samtools view -buh input.bam MT | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bamMT.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bamMT.bam bamMT.bai 

ls -l

#
# Upload results
#
name=`dx describe "$sorted_bam" --name`
name="${name%.bam}"
file_id13=`dx upload bam13.bam -o "$name.13.bam" --brief`
file_id13i=`dx upload bam13.bai -o "$name.13.bai" --brief`
dx-jobutil-add-output "bam13" "$file_id13"
dx-jobutil-add-output "bai13" "$file_id13i"
file_id14=`dx upload bam14.bam -o "$name.14.bam" --brief`
file_id14i=`dx upload bam14.bai -o "$name.14.bai" --brief`
dx-jobutil-add-output "bam14" "$file_id14"
dx-jobutil-add-output "bai14" "$file_id14i"
file_id15=`dx upload bam15.bam -o "$name.15.bam" --brief`
file_id15i=`dx upload bam15.bai -o "$name.15.bai" --brief`
dx-jobutil-add-output "bam15" "$file_id15"
dx-jobutil-add-output "bai15" "$file_id15i"
file_id16=`dx upload bam16.bam -o "$name.16.bam" --brief`
file_id16i=`dx upload bam16.bai -o "$name.16.bai" --brief`
dx-jobutil-add-output "bam16" "$file_id16"
dx-jobutil-add-output "bai16" "$file_id16i"
file_id17=`dx upload bam17.bam -o "$name.17.bam" --brief`
file_id17i=`dx upload bam17.bai -o "$name.17.bai" --brief`
dx-jobutil-add-output "bam17" "$file_id17"
dx-jobutil-add-output "bai17" "$file_id17i"
file_id18=`dx upload bam18.bam -o "$name.18.bam" --brief`
file_id18i=`dx upload bam18.bai -o "$name.18.bai" --brief`
dx-jobutil-add-output "bam18" "$file_id18"
dx-jobutil-add-output "bai18" "$file_id18i"
file_id19=`dx upload bam19.bam -o "$name.19.bam" --brief`
file_id19i=`dx upload bam19.bai -o "$name.19.bai" --brief`
dx-jobutil-add-output "bam19" "$file_id19"
dx-jobutil-add-output "bai19" "$file_id19i"
file_id20=`dx upload bam20.bam -o "$name.20.bam" --brief`
file_id20i=`dx upload bam20.bai -o "$name.20.bai" --brief`
dx-jobutil-add-output "bam20" "$file_id20"
dx-jobutil-add-output "bai20" "$file_id20i"
file_id21=`dx upload bam21.bam -o "$name.21.bam" --brief`
file_id21i=`dx upload bam21.bai -o "$name.21.bai" --brief`
dx-jobutil-add-output "bam21" "$file_id21"
dx-jobutil-add-output "bai21" "$file_id21i"
file_id22=`dx upload bam22.bam -o "$name.22.bam" --brief`
file_id22i=`dx upload bam22.bai -o "$name.22.bai" --brief`
dx-jobutil-add-output "bam22" "$file_id22"
dx-jobutil-add-output "bai22" "$file_id22i"
file_idX=`dx upload bamX.bam -o "$name.X.bam" --brief`
file_idXi=`dx upload bamX.bai -o "$name.X.bai" --brief`
dx-jobutil-add-output "bamX" "$file_idX"
dx-jobutil-add-output "baiX" "$file_idXi"
file_idY=`dx upload bamY.bam -o "$name.Y.bam" --brief`
file_idYi=`dx upload bamY.bai -o "$name.Y.bai" --brief`
dx-jobutil-add-output "bamY" "$file_idY"
dx-jobutil-add-output "baiY" "$file_idYi"
file_idMT=`dx upload bamMT.bam -o "$name.MT.bam" --brief`
file_idMTi=`dx upload bamMT.bai -o "$name.MT.bai" --brief`
dx-jobutil-add-output "bamMT" "$file_idMT"
dx-jobutil-add-output "baiMT" "$file_idMTi"
