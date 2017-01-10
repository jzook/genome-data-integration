# CALLABLE LOCI

# Python script to prepare gatk-callableloci-v3.5-anyref files


#Adjust file name for desired output name
out = open("160922_gatk-callableloci-v3.5-anyref_10X_script.sh", "w") 

#Set all of these values as appropriate 
HG = "HG001"           #genome
plat= "10XGenomics"    #platform
max= "36"              #maxDepth determined by 2x coverage (DP) in vcf for chr 20
abv = "10X"	           #abreviated platform name


chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/GATK_V3.5/gatk-callableloci-v3.5-anyref " +
	 	
		"-iinput_bam=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_HP1_"  + str(chr) + ".bam " +
		"-iinput_bai=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_HP1_"  + str(chr) + ".bai " +
		"-ioutput_prefix=" + HG + "_" + str(chr) + "_" + abv + "_HP1_" + "callableloci " +
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"-iextra_options=\"-L chr" + str(chr) + " -minDepth 20 -mmq 20 -maxDepth "+ max + "\" " +
		"--destination=/" + HG+ "/" + "GRCh38/" + plat + "/CallableLoci_output/" +	
		"\n") 
	
	chr = chr+1 

	

chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/GATK_V3.5/gatk-callableloci-v3.5-anyref " +
	 	
		"-iinput_bam=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_HP2_"  + str(chr) + ".bam " +
		"-iinput_bai=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_HP2_"  + str(chr) + ".bai " +
		"-ioutput_prefix=" + HG + "_" + str(chr) + "_" + abv + "_HP2_" + "callableloci " +
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"-iextra_options=\"-L chr" + str(chr) + " -minDepth 20 -mmq 20 -maxDepth "+ max + "\" " +
		"--destination=/" + HG+ "/" + "GRCh38/" + plat + "/CallableLoci_output/" +	
		"\n") 
	
	chr = chr+1 

out.close()