# Python script to prepare FreeBayes files

#Adjust file name for desired output name
out = open("161003_FreeBayes_Ilmn_script_test.sh", "w") 

#Set all of these values as appropriate 
HG = "HG001"		#genome
plat= "Illumina"    #platform
abv= "Ilmn"			#abbreviated platform
cc= "150bp300x"	    #chemistry and coverage
map="novoalign"	    #mapper
ref="/assets/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa.gz"  #reference

chr = 1 
for i in range(22):  
	
	out.write("dx run -y freebayes " +
	 	
		"-isorted_bams=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_" + cc + "_"+ map + "_" + str(chr) + ".bam " +
		"-ioutput_prefix=" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + abv + cc + "_FB " +
		"-igenotype_qualities=TRUE " +
		"-istandard_filters=FALSE " +
		"-iadvanced_options=\"-F 0.05 -m 0\" " +
		"-igenome_fastagz=" ref + " "	
		"--destination=/" + HG + "/GRCh38/"+ plat +"/FreeBayes_output/" +
		"\n") 
	
	chr = chr+1 

out.close()	

