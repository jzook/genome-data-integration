# CALLABLE LOCI

# Python script to prepare gatk-callableloci-v3.5-anyref files


#Adjust file name for desired output name
out = open("160922_gatk-callableloci-v3.5-anyref_ILMN_script.sh", "w") 

#Set all of these values as appropriate 
HG = "HG001"		#genome
plat= "Illumina"    #platform
max= "560"          #maxDepth determined by 2x coverage (DP) in vcf for chr 20
abv= "Ilmn"         #platform name abbreviated
cc= "150bp300x"		#chemistry and coverage
map = "novoalign"   #mapper

chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/GATK_V3.5/gatk-callableloci-v3.5-anyref " +
	 	
		"-iinput_bam=/" + HG + "/GRCh38/" + plat + "/" + HG + "_GRCh38_" + abv + "_" + cc + "_" + map + "_" + str(chr) + ".bam " +
		"-iinput_bai=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_" + cc + "_" + map + "_" + str(chr) + ".bai " +
		"-ioutput_prefix=" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + cc + "_" + "callableloci " +
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"-iextra_options=\"-L chr" + str(chr) + " -minDepth 20 -mmq 20 -maxDepth "+ max + "\" " +
		"--destination=/" +HG+ "/" + "GRCh38/" + plat + "/CallableLoci_output/" +	
		"\n") 
	
	chr = chr+1 

out.close()	
