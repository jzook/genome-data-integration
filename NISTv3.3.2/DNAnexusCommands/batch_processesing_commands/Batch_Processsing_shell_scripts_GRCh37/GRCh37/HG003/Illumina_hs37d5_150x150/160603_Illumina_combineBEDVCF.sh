#Combine Illumina bed/vcf files for all chromosomes into single bed file

#All ran in less than 10 min

#HG003 run 6/3/16
#Adjustments to options needed for each run:
#directory where bed files are located
#set name of output file, -iprefix
#set output directory, --destination

bed_inputs=""
for l in `dx ls "GIAB:/HG003/Illumina/Illumina_hs37d5_150x150/CallableLoci_output/*.bed"`; do bed_inputs="-ibeds=/HG003/Illumina/Illumina_hs37d5_150x150/CallableLoci_output/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG003_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_callableloci --destination=/HG003/Illumina/Illumina_hs37d5_150x150/CallableLoci_output/

#6/3/16 downloaded HG003_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_callableloci and parse out "callable" regions. Onlycallable file was used in comparisons with high conf calls

#HG003 run 6/3/16
#Combine individual vcf files from GATK into one vcf
#Adjustments to options needed for each run:
#directory where vcf files are located
#set name of output file, -iprefix
#set output directory, --destination

vcf_inputs=""
for l in `dx ls "GIAB:/HG003/Illumina/Illumina_hs37d5_150x150/GATKHC_output/*GATKHC.vcf.gz"`; do vcf_inputs="-ivcfs=/HG003/Illumina/Illumina_hs37d5_150x150/GATKHC_output/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG003_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_GATKHC  --destination=/HG003/Illumina/Illumina_hs37d5_150x150/GATKHC_output/

#HG003 run 6/3/16
#Combine individual vcf files from FreeBayes into one vcf
#Adjustments to options needed for each run:
#directory where vcf files are located
#set name of output file, -iprefix
#set output directory, --destination

vcf_inputs=""
for l in `dx ls "GIAB:/HG003/Illumina/Illumina_hs37d5_150x150/FreeBayes_output/*FB.vcf.gz"`; do vcf_inputs="-ivcfs=/HG003/Illumina/Illumina_hs37d5_150x150/FreeBayes_output/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG003_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_FB  --destination=/HG003/Illumina/Illumina_hs37d5_150x150/FreeBayes_output/