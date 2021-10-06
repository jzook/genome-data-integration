#!/usr/bin/bash

for i in *vcf.gz; do
	gunzip -k ${i} 
	vcftools --site-depth --vcf ${i%.gz} --out ${i}
done
