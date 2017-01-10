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
samtools view -buh input.bam 6 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam6.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam6.bam bam6.bai &
samtools view -buh input.bam 7 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam7.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam7.bam bam7.bai &
samtools view -buh input.bam 8 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam8.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam8.bam bam8.bai &
samtools view -buh input.bam 9 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam9.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam9.bam bam9.bai &
samtools view -buh input.bam 10 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam10.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam10.bam bam10.bai &
samtools view -buh input.bam 11 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam11.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam11.bam bam11.bai &
samtools view -buh input.bam 12 | $java -jar /AddOrReplaceReadGroups.jar I=/dev/stdin O=bam12.bam RGID=$rgid RGLB=$rglb RGPL=$rgpl RGPU=$rgpu RGSM=$rgsm $opts
samtools index bam12.bam bam12.bai 

ls -l

#
# Upload results
#
name=`dx describe "$sorted_bam" --name`
name="${name%.bam}"
file_id6=`dx upload bam6.bam -o "$name.6.bam" --brief`
file_id6i=`dx upload bam6.bai -o "$name.6.bai" --brief`
dx-jobutil-add-output "bam6" "$file_id6"
dx-jobutil-add-output "bai6" "$file_id6i"
file_id7=`dx upload bam7.bam -o "$name.7.bam" --brief`
file_id7i=`dx upload bam7.bai -o "$name.7.bai" --brief`
dx-jobutil-add-output "bam7" "$file_id7"
dx-jobutil-add-output "bai7" "$file_id7i"
file_id8=`dx upload bam8.bam -o "$name.8.bam" --brief`
file_id8i=`dx upload bam8.bai -o "$name.8.bai" --brief`
dx-jobutil-add-output "bam8" "$file_id8"
dx-jobutil-add-output "bai8" "$file_id8i"
file_id9=`dx upload bam9.bam -o "$name.9.bam" --brief`
file_id9i=`dx upload bam9.bai -o "$name.9.bai" --brief`
dx-jobutil-add-output "bam9" "$file_id9"
dx-jobutil-add-output "bai9" "$file_id9i"
file_id10=`dx upload bam10.bam -o "$name.10.bam" --brief`
file_id10i=`dx upload bam10.bai -o "$name.10.bai" --brief`
dx-jobutil-add-output "bam10" "$file_id10"
dx-jobutil-add-output "bai10" "$file_id10i"
file_id11=`dx upload bam11.bam -o "$name.11.bam" --brief`
file_id11i=`dx upload bam11.bai -o "$name.11.bai" --brief`
dx-jobutil-add-output "bam11" "$file_id11"
dx-jobutil-add-output "bai11" "$file_id11i"
file_id12=`dx upload bam12.bam -o "$name.12.bam" --brief`
file_id12i=`dx upload bam12.bai -o "$name.12.bai" --brief`
dx-jobutil-add-output "bam12" "$file_id12"
dx-jobutil-add-output "bai12" "$file_id12i"
