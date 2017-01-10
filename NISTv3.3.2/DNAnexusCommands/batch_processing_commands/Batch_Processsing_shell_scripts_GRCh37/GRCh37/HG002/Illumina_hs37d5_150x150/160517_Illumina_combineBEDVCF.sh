#Combine Illumina bed/vcf files for all chromosomes into single bed file

#All ran in less than 10 min

#HG002 run 5/18/16
#Adjustments to options needed for each run:
#directory where bed files are located
#set name of output file, -iprefix
#set output directory, --destination

bed_inputs=""
for l in `dx ls "GIAB:/HG002/Illumina/CallableLoci_output/*.bed"`; do bed_inputs="-ibeds=/HG002/Illumina/CallableLoci_output/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG002_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_callableloci --destination=/HG002/Illumina/CallableLoci_output/

#5/23/16 downloaded HG002_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_callableloci and parse out "callable" regions. Combined onlycallable using GUI.  Onlycallable file was used in comparisons.

#HG002 run 5/17/16
#Combine individual vcf files from GATK into one vcf
#Adjustments to options needed for each run:
#directory where vcf files are located
#set name of output file, -iprefix
#set output directory, --destination

vcf_inputs=""
for l in `dx ls "GIAB:/HG002/Illumina/GATKHC_output/*GATKHC.vcf.gz"`; do vcf_inputs="-ivcfs=/HG002/Illumina/GATKHC_output/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG002_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_GATKHC  --destination=/HG002/Illumina/GATKHC_output/

#HG002 run 5/17/16
#Combine individual vcf files from FreeBayes into one vcf
#Adjustments to options needed for each run:
#directory where vcf files are located
#set name of output file, -iprefix
#set output directory, --destination

vcf_inputs=""
for l in `dx ls "GIAB:/HG002/Illumina/FreeBayes_output/*FB.vcf.gz"`; do vcf_inputs="-ivcfs=/HG002/Illumina/FreeBayes_output/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG002_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_FB  --destination=/HG002/Illumina/FreeBayes_output/