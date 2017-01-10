# Python script to prepare gatk-genotype-gvcfs-v3.5-anyref files

#Adjust file name for desired output name
out = open("161003_gatk-genotype-gvcfs-v3.5-anyref_ILMN_script.sh", "w") 

#Set all of these values as appropriate 
HG = "HG001"		#genome
plat= "Illumina"    #platform
abv = "Ilmn"        #abreviated platform name
cc = "150bp300x"    #chemistry and/or coverage
map = "novoalign"   #mapper


chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5-anyref " +
	 	
		"-ivcfs=/"+HG+"/GRCh38/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + cc + "_GATKHC_gvcf.vcf.gz " + 
		"-ivcfs=/"+HG+"/GRCh38/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + cc + "_GATKHC_gvcf.vcf.gz.tbi " + 
		"-iprefix=" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + abv + cc + "_" + "GATKHC "  +
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"--destination=/" + HG + "/" + "GRCh38/" + plat + "/GATKHC_output/" +	
		"\n") 
	
	chr = chr+1 

out.close()	

