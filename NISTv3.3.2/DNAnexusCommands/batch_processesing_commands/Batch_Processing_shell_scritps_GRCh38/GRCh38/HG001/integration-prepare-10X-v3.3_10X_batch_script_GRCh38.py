# Python script to prepare integration-prepare-10X-v3.3-anyref files

#Adjust file name for desired output name
out = open("161003_integration-prepare-10X-v3.3-anyref_10X_script.sh", "w") 

#Set all of these values as appropriate 
HG = "HG001"		#genome
plat= "10XGenomics"    #platform
abv = "10X"        #abreviated platform name
ref = "/assets/GRCh38hs38d1noalt.fasta-index.tar.gz"

chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/integration-prepare-10X-v3.3-anyref " +
	 	
		"-igvcf1=/"+HG+"/GRCh38/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_" + abv + "_HP1_GATKHC_gvcf.vcf.gz " + 
		"-igvcftbi1=/"+HG+"/GRCh38/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_" + abv + "_HP1_GATKHC_gvcf.vcf.gz.tbi " + 
		"-igvcf2=/"+HG+"/GRCh38/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_" + abv + "_HP2_GATKHC_gvcf.vcf.gz " + 
		"-igvcftbi2=/"+HG+"/GRCh38/"+ plat +"/GATKHC_output/" + HG + "_" + str(chr) + "_" + abv + "_HP2_GATKHC_gvcf.vcf.gz.tbi " + 
		"-ibed1=/"+HG+"/GRCh38/"+ plat +"/CallableLoci_output/" + HG + "_" + str(chr) + "_" + abv + "_HP1_callableloci.bed "
		"-ibed2=/"+HG+"/GRCh38/"+ plat +"/CallableLoci_output/" + HG + "_" + str(chr) + "_" + abv + "_HP2_callableloci.bed "
		"-iprefix=" + HG + "_" + str(chr) + "_" + abv + "_GATKHCbyhaplo "  +
		"-iref=" + ref + " " +
		"-ichrom=" + "chr" + str(chr)+ " " +
		"-imaxcov=36 " +
		"--destination=/" + HG + "/" + "GRCh38/" + plat + "/Integration_prepare_10X_output_v3.3/" +	
		"\n") 
	
	chr = chr+1 

out.close()	

