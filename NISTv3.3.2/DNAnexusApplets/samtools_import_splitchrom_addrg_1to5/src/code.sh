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
# Run samtools view to select reads from each chromosome and output as new file
#
samtools view -buh input.bam 1 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam1.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam1.bam bam1.bai &
samtools view -buh input.bam 2 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam2.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam2.bam bam2.bai &
samtools view -buh input.bam 3 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam3.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam3.bam bam3.bai &
samtools view -buh input.bam 4 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam4.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam4.bam bam4.bai &
samtools view -buh input.bam 5 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam5.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam5.bam bam5.bai 

ls -l

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
