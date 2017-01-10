# Python script to prepare gatk-genotype-gvcfs-v3.5-anyref files

#Adjust file name for desired output name
out = open("160922_gatk-haplotypecaller-v35-anyref_ILMN_script.sh", "w") 

#Set all of these values as appropriate 
HG = "HG001"		#genome
plat= "Illumina"    #platform
abv = "Ilmn"        #abreviated platform name
cc = "150bp300x"    #chemistry and/or coverage
map = "novoalign"   #mapper


chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/GATK_V3.5/gatk-haplotypecaller-v35-anyref " +
	 	
		"-isorted_bam=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_" + cc + "_" + map + "_" + str(chr) + ".bam " +
		"-isorted_bai=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_" + cc + "_" + map + "_" + str(chr) + ".bai " +
		"-ioutput_prefix=" + HG + "_" + str(chr) + "_GRCh38_" + map + "_" + cc + "_" + "GATKHC_gvcf "  +
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"-iextra_options=\"-L chr" + str(chr) + " -stand_call_conf 2 -stand_emit_conf 2 -A BaseQualityRankSumTest -A ClippingRankSumTest -A Coverage -A FisherStrand -A LowMQ -A RMSMappingQuality -A ReadPosRankSumTest -A StrandOddsRatio -A HomopolymerRun -A TandemRepeatAnnotator\" " +
		"--destination=/" + HG + "/" + "GRCh38/" + plat + "/GATKHC_output/" +	
		"\n") 
	
	chr = chr+1 

out.close()	
