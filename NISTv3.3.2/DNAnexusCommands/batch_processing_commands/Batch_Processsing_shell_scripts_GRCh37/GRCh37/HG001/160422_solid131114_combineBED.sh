#Combine bed files for each chromosome into a single bed file for solid run 131114
bed_inputs=""
for l in `dx ls "GIAB:/NA12878/SOLID/131114_SE75bp/callableLoci_output/*.bed"`; do bed_inputs="-ibeds=/NA12878/SOLID/131114_SE75bp/callableLoci_output/$l $bed_inputs"; done
dx run GIAB:/Workflow/bed-combineallchrom $bed_inputs -iprefix=HG001_ALLCHROM_hg19_solid5500_SE75bp_callableloci --destination=/NA12878/SOLID/131114_SE75bp/callableLoci_output/