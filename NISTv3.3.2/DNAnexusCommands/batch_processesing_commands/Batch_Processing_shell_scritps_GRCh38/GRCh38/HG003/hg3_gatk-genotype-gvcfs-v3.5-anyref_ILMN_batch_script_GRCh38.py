# Python script to prepare gatk-genotype-gvcfs-v3.5-anyref files

#Adjust file name for desired output name
out = open("161219_gatk-genotype-gvcfs-v3.5-anyref_ILMN150_hg.sh", "w") 

#Set all of these values as appropriate 
ref= "GRCh38"
HG = "HG003"		#genome
platform= "Illumina"
calls= "Illumina_GRCh38_150x150"
cc = "Ilmn150bp300X"    #chemistry and/or coverage
map = "novoalign"   #mapper
vc= "sentieonHC"


chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/GATK_V3.5/gatk-genotype-gvcfs-v3.5-anyref " +
	 	
		"-ivcfs=" + HG + "/"+ ref + "/" + platform + "/" + calls + "/Sentieon_output/" + HG + "_" + str(chr) + "_GRCh38_"+ map + "_" + cc + "_" + vc + "_gvcf.vcf.gz " + 
		"-ivcfs=" + HG + "/"+ ref + "/" + platform + "/" + calls + "/Sentieon_output/" + HG + "_" + str(chr) + "_GRCh38_"+ map + "_" + cc + "_" + vc + "_gvcf.vcf.gz.tbi " + 
		"-iprefix=" + HG + "_" + str(chr) + "_GRCh38_"+ map + "_" + cc + "_" + vc + " " +
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"--destination=/" + HG + "/"+ ref + "/" + platform + "/" + calls + "/Sentieon_output/" +	
		"\n") 
	
	chr = chr+1 

out.close()	

