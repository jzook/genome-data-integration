# Python script to prepare integration-prepare-gatkhc-v3.3.1-anyref files

#Adjust file name for desired output name
out = open("161005_integration-prepare-gatkhc-v3.3.1-anyref_ILMN_script.sh", "w") 


#Set all of these values as appropriate 
HG = "HG001"		#genome
plat= "Illumina"    #platform
abv = "Ilmn"        #abreviated platform name
cc = "150bp300x"    #chemistry and/or coverage
map = "novoalign"   #mapper
ref = "/assets/hs37d5.fasta-index.tar.gz"

chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/integration-prepare-gatkhc-v3.3.1-anyref " +
	 	
		"-igvcf=/"+HG+"/GRCh37/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_hs37d5_" + map + "_" + "Ilmn150bp300X_GATKHC_gvcf.vcf.gz " + 
		"-igvcftbi=/"+HG+"/GRCh37/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_hs37d5_" + map + "_" + "Ilmn150bp300X_GATKHC_gvcf.vcf.gz " + 
		#"-iprefix=" + HG + "_" + str(chr) + "_GRCh37_" + map + "_" + abv + cc + "_" + "GATKHC "  +
		"-iref=" + ref + " " +
		"-ichrom=" + str(chr)+ " " +
		"--destination=/" + HG + "/GRCh37/" + plat + "/Integration_prepare_GATKHC_v.3.3.1/" +	
		"\n") 
	
	chr = chr+1 

out.close()	
