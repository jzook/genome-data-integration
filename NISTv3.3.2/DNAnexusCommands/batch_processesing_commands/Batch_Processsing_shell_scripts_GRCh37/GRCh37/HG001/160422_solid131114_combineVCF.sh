#Combine vcf files for each chromosome into a single vcf file for solid run 131114
bed_inputs=""
for l in `dx ls "GIAB:/NA12878/SOLID/131114_SE75bp/GATKHC_output/*GATKHC.vcf.gz"`; do vcf_inputs="-ivcfs=/NA12878/SOLID/131114_SE75bp/GATKHC_output/$l $vcf_inputs"; done
dx run GIAB:/Workflow/vcf-combineallchrom $vcf_inputs -iprefix=HG001_ALLCHROM_hg19_solid5500_SE75bp_GATKHC --destination=/NA12878/SOLID/131114_SE75bp/GATKHC_output