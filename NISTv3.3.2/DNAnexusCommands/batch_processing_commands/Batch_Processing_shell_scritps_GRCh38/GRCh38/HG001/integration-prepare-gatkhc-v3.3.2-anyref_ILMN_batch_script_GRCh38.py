# Python script to prepare integration-prepare-gatkhc-v3.3.1-anyref files

#Adjust file name for desired output name
out = open("161107_integration-prepare-gatkhc-v3.3.2-anyref_ILMN_script.sh", "w") 

#Set all of these values as appropriate 
HG = "HG001"		#genome
plat= "Illumina"    #platform
abv = "Ilmn"        #abreviated platform name
cc = "150bp300x"    #chemistry and/or coverage
map = "novoalign"   #mapper
ref = "/assets/GRCh38hs38d1noalt.fasta-index.tar.gz"

chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/integration-prepare-gatkhc-v3.3.2-anyref " +
	 	
		"-igvcf=/"+HG+"/GRCh38/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + cc + "_GATKHC_gvcf.vcf.gz " + 
		"-igvcftbi=/"+HG+"/GRCh38/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + cc + "_GATKHC_gvcf.vcf.gz.tbi " + 
		#"-iprefix=" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + abv + cc + "_" + "GATKHC "  +
		"-iref=" + ref + " " +
		"-ichrom=chr" + str(chr)+ " " +
		"--destination=/" + HG + "/" + "GRCh38/" + plat + "/Integration_prepare_GATKHC_v.3.3.2/" +	
		"\n") 
	
	chr = chr+1 

out.close()	
