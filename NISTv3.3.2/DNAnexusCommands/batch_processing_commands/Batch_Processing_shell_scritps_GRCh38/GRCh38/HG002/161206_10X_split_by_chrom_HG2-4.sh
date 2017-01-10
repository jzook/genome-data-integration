#Split bam by chromosome
#All Run 12/6/16

#Read group ID (rgid): HP1 or HP2
#Read group LB (rglb): 10X
#Read group PL (rgpl): illumina
#Read group PU (rgpu): all
#Read group SM (rgsm): NA24385

#HG002
dx run -y GIAB:/Workflow/samtools_splitchrom_addrg_withchr -isorted_bam=/HG002/GRCh38/10XGenomics/NA24385_GRCh38_phased_possorted_bam_HP1.bam -iindex_bai=/HG002/GRCh38/10XGenomics/NA24385_GRCh38_phased_possorted_bam_HP1.bam.bai -iprefix=HG002_GRCh38_10X_HP1_ -irgid=HP1 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=NA24385 --destination=/HG002/GRCh38/10XGenomics/ --instance-type=mem2_hdd2_x1
dx run -y GIAB:/Workflow/samtools_splitchrom_addrg_withchr -isorted_bam=/HG002/GRCh38/10XGenomics/NA24385_GRCh38_phased_possorted_bam_HP2.bam -iindex_bai=/HG002/GRCh38/10XGenomics/NA24385_GRCh38_phased_possorted_bam_HP2.bam.bai -iprefix=HG002_GRCh38_10X_HP2_ -irgid=HP2 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=NA24385 --destination=/HG002/GRCh38/10XGenomics/ --instance-type=mem2_hdd2_x1

#HG003 
#reran HP1 12/21/16  accidentally used wrong index
dx run -y GIAB:/Workflow/samtools_splitchrom_addrg_withchr -isorted_bam=/HG003/GRCh38/10XGenomics/NA24149_GRCh38_phased_possorted_bam_HP1.bam -iindex_bai=/HG003/GRCh38/10XGenomics/NA24149_GRCh38_phased_possorted_bam_HP1.bam.bai -iprefix=HG003_GRCh38_10X_HP1_ -irgid=HP1 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=NA24149 --destination=/HG003/GRCh38/10XGenomics/ --instance-type=mem2_hdd2_x1
dx run -y GIAB:/Workflow/samtools_splitchrom_addrg_withchr -isorted_bam=/HG003/GRCh38/10XGenomics/NA24149_GRCh38_phased_possorted_bam_HP2.bam -iindex_bai=/HG003/GRCh38/10XGenomics/NA24149_GRCh38_phased_possorted_bam_HP2.bam.bai -iprefix=HG003_GRCh38_10X_HP2_ -irgid=HP2 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=NA24149 --destination=/HG003/GRCh38/10XGenomics/ --instance-type=mem2_hdd2_x1

#HG004
dx run -y GIAB:/Workflow/samtools_splitchrom_addrg_withchr -isorted_bam=/HG004/GRCh38/10XGenomics/NA24143_GRCh38_phased_possorted_bam_HP1.bam -iindex_bai=/HG004/GRCh38/10XGenomics/NA24143_GRCh38_phased_possorted_bam_HP1.bam.bai -iprefix=HG004_GRCh38_10X_HP1_ -irgid=HP1 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=NA24143 --destination=/HG004/GRCh38/10XGenomics/ --instance-type=mem2_hdd2_x1
dx run -y GIAB:/Workflow/samtools_splitchrom_addrg_withchr -isorted_bam=/HG004/GRCh38/10XGenomics/NA24143_GRCh38_phased_possorted_bam_HP2.bam -iindex_bai=/HG004/GRCh38/10XGenomics/NA24143_GRCh38_phased_possorted_bam_HP2.bam.bai -iprefix=HG004_GRCh38_10X_HP2_ -irgid=HP2 -irglb=10X -irgpl=illumina -irgpu=all -irgsm=NA24143 --destination=/HG004/GRCh38/10XGenomics/ --instance-type=mem2_hdd2_x1

