# Prep 10X HP1/HP2 bam mapped to GRCh38
# Run on 9/21/16, ran 9 hr 48 min

dx run -y GIAB:/Workflow/samtools_splitchrom_addrg_withchr -isorted_bam=/NA12878/GRCh38/10XGenomics/NA12878_GRCh38_HP1.bam -iindex_bai=/NA12878/GRCh38/10XGenomics/NA12878_GRCh38_HP1.bam.bai -iprefix=HG001_GRCh38_10X_HP1_ -irgid=HP1 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=NA12878  --destination=/NA12878/GRCh38/10XGenomics/ --instance-type=mem2_hdd2_x1
dx run -y GIAB:/Workflow/samtools_splitchrom_addrg_withchr -isorted_bam=/NA12878/GRCh38/10XGenomics/NA12878_GRCh38_HP2.bam -iindex_bai=/NA12878/GRCh38/10XGenomics/NA12878_GRCh38_HP2.bam.bai -iprefix=HG001_GRCh38_10X_HP2_ -irgid=HP2 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=NA12878  --destination=/NA12878/GRCh38/10XGenomics/ --instance-type=mem2_hdd2_x1
