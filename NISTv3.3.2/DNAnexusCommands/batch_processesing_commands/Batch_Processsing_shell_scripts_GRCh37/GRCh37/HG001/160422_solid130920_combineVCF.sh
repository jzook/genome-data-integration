#Combine vcf files for each chromosome into a single vcf file for solid run 130920

vcf_inputs=""
for l in `dx ls "GIAB:/NA12878/SOLID/130920_PE50x50bp/GATKHC_output/*GATKHC.vcf.gz"`; do vcf_inputs="-ivcfs=/NA12878/SOLID/130920_PE50x50bp/GATKHC_output/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG001_ALLCHROM_hg19_solid5500_PE50x50bp_GATKHC  --destination=/NA12878/SOLID/130920_PE50x50bp/GATKHC_output/