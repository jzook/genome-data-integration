#Before running fill in...
#out file name
#intdate
#cstable path




# Python script to prepare integration commands for chr 1-22
# This script was based on integration v3.3.2 calls and callset tables
# To run :  $ python Integration_batch_script.sh
# Integration script file output to directory script was run in

#Adjust file name for desired output name
out = open("170103_nist-integration-v3.3.2-anyref_hg4_b38.sh", "w") 

#Set all of these values as appropriate for current integration
int="v.3.3.2"                              #Integration version
HG="HG004" 								#genome number
intdate= "170103" 						#date integration is being run, yymmdd
CG="GS000037262" 						#identifier for Complete Genomics dataset
ion="20141120.NA24143" 					#identifier for Ion dataset
callsets= "CG-Illfb-IllsentieonHC-Ion-10XsentieonHC"  #callsets used for integration
RM= "RM8392"        						  #RM#
datasets = "CG_Illumina_Ion_10X"    #data set for callset tables

#confirm these annotation files are being used for current integration, if not add or delete files and adjust in iannot section in BOTH loops
iannot1= "-iannotations=/Annotation_files/GATK_Annotations_160509.txt"
iannot2= "-iannotations=/Annotation_files/Freebayes_Annotations_160408.txt"
iannot3= "-iannotations=/Annotation_files/CG_Annotations_160408.txt"
iannot4= "-iannotations=/Annotation_files/Ion_Annotations_160408.txt"

#set path to desired callset table
cstable= "170103_HG004_Callset_tables_GC_Ilmn_10X_Ion_for_GRCh38_v3.3.2/"

#confirm bed files for filtering are being used for current integration, if not add or delete files and adjust in ifiltbed section in BOTH loops
ifiltbed1= "-ifiltbeds="+HG+"/GRCh38/remapped_HG002_HG003_HG004_allsvs_merged.bed"
ifiltbed2= "-ifiltbeds=filtbeds/GRCh38/AllRepeats_gt200bp_gt95identity_merged.bed.gz"
ifiltbed3= "-ifiltbeds=filtbeds/GRCh38/hg38_self_chain_withalts_gt10k.bed.gz"
ifiltbed4= "-ifiltbeds=filtbeds/GRCh38/SimpleRepeat_imperfecthomopolgt10_slop5.bed"
ifiltbed5= "-ifiltbeds=filtbeds/GRCh38/remapped_superdupsmerged_all_sort.bed"
ifiltbed6= "-ifiltbeds=filtbeds/GRCh38/AllRepeats_lt51bp_gt95identity_merged.bed.gz"
ifiltbed7= "-ifiltbeds=filtbeds/GRCh38/AllRepeats_51to200bp_gt95identity_merged.bed.gz"

#Confirm the following callsets are being used for current integration: 
	 # CG
	 # Ilmn 150x150 300X GATKHC
	 # Ilmn 150x150 FB
	 # IonExome
	 # 10X	 		
     # Ilmn 250x250 GATKHC
	 # Ilmn 250x250 FB
	 # Ilmn MatePair GATKHC
	 # Ilmn MatePair FB
	 # SOLiD PE50x50bp GATKHC
	 # SOLiD SE75bp GATKHC
	 # If callsets vary delete and or add vcf/beds to BOTH loops
	 
# First for loop for chromosome 1-15 and has call for additional memory, adjust as neccessary


chr = 1 
for i in range(15):  
	
	out.write("dx run -y GIAB:/Workflow/nist-integration-v3.3.2-anyref "+
	 	
		"-ivcfs=/"+HG+"/GRCh38/Complete_Genomics/converted_Integration_prepare_cg_output/"+ HG + "_" + str(chr) + "_convGRCh38_CG_vcfBeta-" + CG + "-ASM.vcf.gz " + 
		
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_150x150/Sentieon_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn150bp300X_sentieonHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_150x150/FreeBayes_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn150bp300X_FB.vcf.gz " + 
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_250x250/Sentieon_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn250x250_sentieonHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_250x250/FreeBayes_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn250x250_FB.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_MatePair/Sentieon_output/"+HG+"_"+str(chr)+"_GRCh38_bwa_mem_IlmnMatePair_sentieonHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_MatePair/FreeBayes_output/"+HG+"_"+str(chr)+"_GRCh38_bwa_mem_IlmnMatePair_FB.vcf.gz " +

		"-ivcfs=/"+HG+"/GRCh38/IonExome/converted_Integration_prepare_ion_and_callableloci_output/" + HG + "_" + str(chr) + "_convGRCh38_Ion_AmpliseqExome."+ ion + ".vcf.gz " + 
		
		"-ivcfs=/"+HG+"/GRCh38/10XGenomics/Integration_prepare_10X_output_v3.3/"+HG+"_"+str(chr)+"_GRCh38_10X_sentieonHCbyhaplo.vcf.gz "+
				
		
		"-ibeds=/"+HG+"/GRCh38/Complete_Genomics/converted_Integration_prepare_cg_output/"+ HG + "_" + str(chr) + "_convGRCh38_CG_vcfBeta-" + CG + "-ASM.bed " +
		
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_150x150/Integration_prepare_sentieon_v.3.3.2/" + HG + "_" + str(chr) + "_GRCh38_novoalign_Ilmn150bp300X_sentieonHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_150x150/CallableLoci_output/" + HG + "_" + str(chr) + "_GRCh38_novoalign_Ilmn150bp300x_callableloci.bed " +
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_250x250/Integration_prepare_sentieon_v.3.3.2/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn250x250_sentieonHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_250x250/CallableLoci_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn250x250_callableloci.bed " +
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_MatePair/Integration_prepare_sentieon_v.3.3.2/"+HG+"_"+str(chr)+"_GRCh38_bwa_mem_IlmnMatePair_sentieonHC_gvcf_callable.bed "+
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_MatePair/CallableLoci_output/"+HG+"_"+str(chr)+"_GRCh38_bwa_mem_IlmnMatePair_callableloci.bed "+
		
		"-ibeds=/"+HG+"/GRCh38/IonExome/converted_Integration_prepare_ion_and_callableloci_output/" + HG + "_" + str(chr) + "_convGRCh38_Ion_AmpliseqExome."+ ion + ".bed  " +
		
		"-ibeds=/"+HG+"/GRCh38/10XGenomics/Integration_prepare_10X_output_v3.3/" + HG + "_" + str(chr) + "_GRCh38_10X_sentieonHCbyhaplo_callable.bed " +
		 
		iannot1 + " " + iannot2 + " " + iannot3 + " " + iannot4 + " " +
		 
		"-icallsettable=/callset_tables/" + cstable + HG + "_" + RM + "_Datasets_" + datasets + "_Files_GRCh38_" + str(chr) + ".txt " +
		 
		ifiltbed1 + " " + ifiltbed2 + " " + ifiltbed3 + " " + ifiltbed4 + " " + ifiltbed5 + " " + ifiltbed6 + " " + ifiltbed7 + " " +
		
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"-irtgsdf=/assets/rtgsdf38.tar.gz " +
		"-irefn=/filtbeds/GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N.bed " + 
		"-ichrom=" + "chr" + str(chr) + " " + 
		"-iprefix="+ HG + "_GRCh38_GIAB_highconf_"+ callsets + "_" + str(chr) + "_" + int + " " + 
		"--destination=/" + HG + "/GRCh38/Integration_"+ int + "_output/" + intdate + "_"+ callsets + "_" + int + "/ " +
		"--instance-type=mem3_hdd2_x2"+ 
		"\n") 
	
	chr = chr+1 

# Second for loop for smaller chromosomes 16-22, default memory used
chr = 16
for i in range(7):  
	
	out.write("dx run -y GIAB:/Workflow/nist-integration-v3.3.2-anyref "+
	 	
		"-ivcfs=/"+HG+"/GRCh38/Complete_Genomics/converted_Integration_prepare_cg_output/"+ HG + "_" + str(chr) + "_convGRCh38_CG_vcfBeta-" + CG + "-ASM.vcf.gz " + 
		
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_150x150/Sentieon_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn150bp300X_sentieonHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_150x150/FreeBayes_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn150bp300X_FB.vcf.gz " + 
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_250x250/Sentieon_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn250x250_sentieonHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_250x250/FreeBayes_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn250x250_FB.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_MatePair/Sentieon_output/"+HG+"_"+str(chr)+"_GRCh38_bwa_mem_IlmnMatePair_sentieonHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_MatePair/FreeBayes_output/"+HG+"_"+str(chr)+"_GRCh38_bwa_mem_IlmnMatePair_FB.vcf.gz " +

		"-ivcfs=/"+HG+"/GRCh38/IonExome/converted_Integration_prepare_ion_and_callableloci_output/" + HG + "_" + str(chr) + "_convGRCh38_Ion_AmpliseqExome."+ ion + ".vcf.gz " + 
		
		"-ivcfs=/"+HG+"/GRCh38/10XGenomics/Integration_prepare_10X_output_v3.3/"+HG+"_"+str(chr)+"_GRCh38_10X_sentieonHCbyhaplo.vcf.gz "+
				
		
		"-ibeds=/"+HG+"/GRCh38/Complete_Genomics/converted_Integration_prepare_cg_output/"+ HG + "_" + str(chr) + "_convGRCh38_CG_vcfBeta-" + CG + "-ASM.bed " +
		
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_150x150/Integration_prepare_sentieon_v.3.3.2/" + HG + "_" + str(chr) + "_GRCh38_novoalign_Ilmn150bp300X_sentieonHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_150x150/CallableLoci_output/" + HG + "_" + str(chr) + "_GRCh38_novoalign_Ilmn150bp300x_callableloci.bed " +
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_250x250/Integration_prepare_sentieon_v.3.3.2/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn250x250_sentieonHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_250x250/CallableLoci_output/"+HG+"_"+str(chr)+"_GRCh38_novoalign_Ilmn250x250_callableloci.bed " +
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_MatePair/Integration_prepare_sentieon_v.3.3.2/"+HG+"_"+str(chr)+"_GRCh38_bwa_mem_IlmnMatePair_sentieonHC_gvcf_callable.bed "+
		"-ibeds=/"+HG+"/GRCh38/Illumina/Illumina_GRCh38_MatePair/CallableLoci_output/"+HG+"_"+str(chr)+"_GRCh38_bwa_mem_IlmnMatePair_callableloci.bed "+
		
		"-ibeds=/"+HG+"/GRCh38/IonExome/converted_Integration_prepare_ion_and_callableloci_output/" + HG + "_" + str(chr) + "_convGRCh38_Ion_AmpliseqExome."+ ion + ".bed  " +
		
		"-ibeds=/"+HG+"/GRCh38/10XGenomics/Integration_prepare_10X_output_v3.3/" + HG + "_" + str(chr) + "_GRCh38_10X_sentieonHCbyhaplo_callable.bed " +
		 
		iannot1 + " " + iannot2 + " " + iannot3 + " " + iannot4 + " " +
		 
		"-icallsettable=/callset_tables/" + cstable + HG + "_" + RM + "_Datasets_" + datasets + "_Files_GRCh38_" + str(chr) + ".txt " +
		 
		ifiltbed1 + " " + ifiltbed2 + " " + ifiltbed3 + " " + ifiltbed4 + " " + ifiltbed5 + " " + ifiltbed6 + " " + ifiltbed7 + " " +
		
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"-irtgsdf=/assets/rtgsdf38.tar.gz " +
		"-irefn=/filtbeds/GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N.bed " + 
		"-ichrom=" + "chr" + str(chr) + " " + 
		"-iprefix="+ HG + "_GRCh38_GIAB_highconf_"+ callsets + "_" + str(chr) + "_" + int + " " + 
		"--destination=/" + HG + "/GRCh38/Integration_"+ int + "_output/" + intdate + "_"+ callsets + "_" + int + "/ " +
		"\n") 
		 
	chr = chr+1 	
out.close()	
