#HG004 combine
#All combine run -y 1/3/17

#*_callablemultinter_gtO.bed -- BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/*_callablemultinter_gt0.bed"`; do bed_inputs="-ibeds=GIAB:/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/$l $bed_inputs"; done
dx run -y GIAB:/Workflow/bed-combineallchromwithchr $bed_inputs -iprefix=HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_callablemultinter_gt0 --destination=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/

#*_filteredsites.bed  --BED file with 50bp regions around not-high-confidence sites
bed_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/*_filteredsites.bed"`; do bed_inputs="-ibeds=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/$l $bed_inputs"; done
dx run -y GIAB:/Workflow/bed-combineallchromwithchr $bed_inputs -iprefix=HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_filteredsites --destination=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/

#*_highconf.bed  --High-confidence BED file with at least 1 technology callable
bed_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/*_highconf.bed"`; do bed_inputs="-ibeds=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/$l $bed_inputs"; done
dx run -y GIAB:/Workflow/bed-combineallchromwithchr $bed_inputs -iprefix=HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf --destination=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/

#*_all.vcf.gz  -- All integrated variants, including uncertain variants as filtered rows.
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/*v.3.3.2_all.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/$l $vcf_inputs"; done
dx run -y GIAB:/Workflow/vcf-combineallchromwithchr $vcf_inputs -iprefix=HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_all  --destination=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/ --instance-type=mem2_hdd2_x2

#*_highconf_vcf.gz -- High-confidence integrated variants.
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/*highconf.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/$l $vcf_inputs"; done
dx run -y GIAB:/Workflow/vcf-combineallchromwithchr $vcf_inputs -iprefix=HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf --destination=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/ --instance-type=mem2_hdd2_x2

#*_annotated.vcf.gz 
vcf_inputs=""
for l in `dx ls "GIAB:/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/*annotated.vcf.gz"`; do vcf_inputs="-ivcfs=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/$l $vcf_inputs"; done
dx run -y GIAB:/Workflow/vcf-combineallchromwithchr $vcf_inputs -iprefix=HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_annotated --destination=/HG004/GRCh38/Integration_v.3.3.2_output/170103_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_v.3.3.2/ --instance-type=mem2_hdd2_x2