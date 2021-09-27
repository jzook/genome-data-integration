
#
# VcfClassifyUsingFilters_v3.pl - arbitrate uncertain genotype calls using filters and callable regions, determined differently for each technology
# v3.1 - separate training for SNPs and non-SNPs
# v3.1.1 - allow different annotations for SNPs and non-SNPs

import random 
import sys
import numpy
import re

class OneClassFiltering(object):
	line = None

	#create header
	header = None

	trainvcf = None
	sensvcf = None
	annotationsToFilter = None
	outputFileStart = None
	callsetName = None

	trainvcf_file = None
	sensvcf_file = None
	annotationsToFilter_file = None

	output_filtered_bed_file = None
	output_filtered_vcf_file = None

	#Read in information about each callset
	annots = None
	tails = None
	infoformat = None
	annsnp = None
	annindel = None
	fields = None
	annotation_lines = None
	annotno = None

	x = None
	vgraph = None
	vcf_end_of_header_line = None
	trainvcf_lines = None
	vcf_end_of_header_line = None

	ln_vcf_file = None
	rand_row_number = None
	annotMatrixSNP = None
	annotMatrixIndel = None

	rownosnp = None
	rownoindel = None
	skiprows = None

	annotmatrixsort = None

	cutoffleftsnp = None
	cutoffrightsnp = None
	cutoffleftindel = None
	cutoffrightindel = None

	annotmatrix1 = None
	rowno2 = None

	trainvcf_file = None
	sensvcf_file = None
	annotationsToFilter_file = None
	output_filtered_bed_file = None
	output_filtered_vcf_file = None


	def __init__(self, trainVCF, sensitiveVCF, annotationsToFilter, outputFileStart, callsetName):
		self.line = ""

		#create header
		self.header = ""

		self.trainvcf = trainVCF
		self.sensvcf = sensitiveVCF
		self.annotationsToFilter = annotationsToFilter
		self.outputFileStart = outputFileStart
		self.callsetName = callsetName

		self.trainvcf_file = None
		self.sensvcf_file = None
		self.annotationsToFilter_file = None

		self.output_filtered_bed_file = None
		self.output_filtered_vcf_file = None

		#Read in information about each callset
		self.annots = None
		self.tails = None
		self.infoformat = None
		self.annsnp = None
		self.annindel = None
		self.fields = None
		self.annotation_lines = None
		self.annotno = None

		self.x = None
		self.vgraph = None
		self.vcf_end_of_header_line = None
		self.trainvcf_lines = None
		self.vcf_end_of_header_line = None

		self.ln_vcf_file = None
		self.rand_row_number = None
		self.annotMatrixSNP = None
		self.annotMatrixIndel = None

		self.rownosnp = None
		self.rownoindel = None
		self.skiprows = None

		self.annotmatrixsort = None

		self.cutoffleftsnp = None
		self.cutoffrightsnp = None
		self.cutoffleftindel = None
		self.cutoffrightindel = None

		self.annotmatrix1 = None
		self.rowno2 = None

		try:
			self.trainvcf_file = open(self.trainvcf, "r")
		except:
			exit()

		try:
			self.sensvcf_file = open(self.sensvcf, "r")
		except:
			exit()

		try:
			self.annotationsToFilter_file = open(self.annotationsToFilter, "r")
		except:
			exit()

		try:
			self.output_filtered_bed_file = open(self.outputFileStart + "_filtered.bed", "w+")
		except:
			exit()

		try:
			self.output_filtered_vcf_file = open(self.outputFileStart + "_filtered.vcf", "w+")
		except:
			exit()

	def	parse_annotations(self):

		#Read in information about each callset
		self.annots = [""]*20
		self.tails = [""]*20
		self.infoformat = [""]*20
		self.annsnp = [""]*20
		self.annindel = [""]*20
		self.fields = []
		self.annotation_lines = self.annotationsToFilter_file.readlines() #skip header line
		self.annotno=0

		for line_idx in range(1, len(self.annotation_lines)):
			# Split up the line into an array
			self.fields = self.annotation_lines[line_idx].split("\t")
			self.annots[self.annotno]=self.fields[0]
			self.tails[self.annotno]=self.fields[1]
			self.infoformat[self.annotno]=self.fields[2]
			self.annsnp[self.annotno]=self.fields[3]
			if re.search("(.*)\n", self.fields[4]):
				self.annindel[self.annotno]=re.split("(.*)\n", self.fields[4])[1]
			else:
				self.annindel[self.annotno]=self.fields[4]
			self.annotno = self.annotno + 1

		self.x = 0
		self.vgraph = 0
		self.vcf_end_of_header_line = 0
		self.trainvcf_lines = self.trainvcf_file.readlines()
		for li in self.trainvcf_lines:
			self.vcf_end_of_header_line = self.vcf_end_of_header_line + 1
			if re.search("ID=BT,", li):
				self.vgraph = 1
			if not(re.search("^\#\#",  li)):
				self.x = 1
			if self.x == 1:
				break


		self.ln_vcf_file = len(self.trainvcf_lines)
		self.rand_row_number = float(40000.0/self.ln_vcf_file)

		#Initialize 2D array to store annotations from SNPs and indels	
		self.annotMatrixSNP = numpy.full((self.annotno,40000), 1000000.0)
		self.annotMatrixIndel = numpy.full((self.annotno,40000), 1000000.0)

	def train_classifier(self): 
		self.rownosnp = 0
		self.rownoindel = 0
		self.skiprows = 0

		for line_idx in range(self.vcf_end_of_header_line, len(self.trainvcf_lines)):
			# Split up the line into an array
			self.fields = self.trainvcf_lines[line_idx].split("\t")            
			#If length of REF or any ALT is greater than one, then treat as not a snp (i.e., an indel)
			self.nonsnp = 0
			self.ref = self.fields[3]
			self.alt = self.fields[4]
			self.alts = self.alt.split(",")
			if len(self.ref) > 1:
				self.nonsnp = 1
			
			if self.nonsnp == 0:
				for j in range(0, len(self.alts)):
					if len(self.alts[j]) > 1:
						self.nonsnp = 1
						break

			random_num = random.random()

			#randomly select 40000 rows for SNPs
			if (self.nonsnp == 0 and (self.rownosnp >= 40000 or random_num > self.rand_row_number)) or (self.nonsnp == 1 and (self.rownoindel >= 40000)):
				self.skiprows = self.skiprows + 1
				continue

			self.info = self.fields[7]
			self.format_field = self.fields[8]
			self.formats = self.format_field.split(":")
			self.chars = self.fields[9].split(":")

			#First, check if # of technologies is >1 if the input is from the vgraph multimerge
			if self.vgraph == 1:
				self.techno = -1
				for k in range(0, len(self.chars)):
					if self.formats[k] == "BT":
						self.charfld1 = self.chars[k].split(",")
						if re.search("(.*)\n", self.charfld1[0]):
							self.techno = re.split("(.*)\n", self.charfld1[0])[1]
						else:
							self.techno = self.charfld1[0]
				if (self.techno == 0 or self.techno == 1 or self.techno == -1):
					self.skiprows = self.skiprows + 1
					continue 

			#Add annotations from this line to the matrix of annotations
			self.annotval = -1

			for j in range(0,self.annotno):
				if self.infoformat[j]=="0":  #field is in INFO
					self.infoflds = self.info.split(";")
					#print("infoflds: " + str(self.infoflds))
					for k in range(0, len(self.infoflds)): # might be +1 -- k<=$#infoflds
						self.infofld = self.infoflds[k].split("=")
						if len(self.infofld)==2 and self.infofld[0] == self.annots[j]: 
							self.infofld1 = self.infofld[1].split(",") # use minimum of first 3 values if field is multivalued (e.g., multiallelic GATK sites)
							self.annotval = self.infofld1[0]
							if(len(self.infofld1) > 1 and float(self.annotval) > float(self.infofld1[1])):
								self.annotval = self.infofld1[1]
							if(len(self.infofld1) > 2 and float(self.annotval) > float(self.infofld1[2])):
								self.annotval = self.infofld1[2]
				else: #field is in FORMAT	
					for k in range(0, len(self.chars)): # might be +1 -- ; $k<=$#chars
						if self.formats[k] == self.annots[j]:
							self.charfld1=self.chars[k].split(",") #use minimum of first 3 values if field is multivalued (e.g., CGA_CEHQ or multiallelic GATK sites)
							if re.search("(.*)\n", self.charfld1[0]):
								self.annotval = re.split("(.*)\n", self.charfld1[0])[1]
							else: 
								self.annotval = self.charfld1[0]
								if len(self.charfld1) > 1:
									if re.search("(.*)\n", self.charfld1[1]):
										if float(self.annotval) > float(re.split("(.*)\n", self.charfld1[1])[1]): 
											self.annotval = re.split("(.*)\n", self.charfld1[1])[1]
									else: 
										if float(self.annotval) > float(self.charfld1[1]):
											self.annotval = self.charfld1[1]
									
									if len(self.charfld1)>2:
										if re.search("(.*)\n",self.charfld1[2]):
											if float(self.annotval) > float(re.split("(.*)\n", self.charfld1[2])[1]):
												self.annotval = re.split("(.*)\n", self.charfld1[2])[1]
										else: 
											if float(self.annotval) > float(self.charfld1[2]):
												self.annotval = self.charfld1[2]
														
				#Put annotation value in SNP or indel matrix
				if self.nonsnp==0: 
					self.annotMatrixSNP[j,self.rownosnp]=float(self.annotval)
				else:
					self.annotMatrixIndel[j,self.rownoindel]=float(self.annotval)
			
			if self.nonsnp==0:
				self.rownosnp = self.rownosnp + 1
			else:
				self.rownoindel = self.rownoindel + 1
			
			
		self.trainvcf_file.close()

		#Step 2: For each annotation, sort the values from smallest to largest and find the appropriate cut-off(s)
		self.annotmatrixsort = None

		self.cutoffleftsnp = [-1000000.0]*self.annotno
		self.cutoffrightsnp = [1000000.0]*self.annotno
		self.cutoffleftindel = [-1000000.0]*self.annotno
		self.cutoffrightindel = [1000000.0]*self.annotno
		for j in range(0,self.annotno):
			#First, find cut-offs for SNPs
			self.annotmatrix1 = [1000000.0]*self.rownosnp
			self.rowno2 = 0
			while self.rowno2 < self.rownosnp: #first make new 1D array with column of interest
				self.annotmatrix1[self.rowno2] = self.annotMatrixSNP[j,self.rowno2]
				self.rowno2 = self.rowno2 + 1
			self.annotmatrix1.sort()
			self.annotmatrixsort = self.annotmatrix1 #sort

			#Finding the number of rows that had a value for this annotation
			self.rowno2 = 0
			while self.rowno2 < self.rownosnp:
				if self.annotmatrixsort[self.rowno2] > 999999:
					break
				self.rowno2 = self.rowno2 + 1

			#Find cut-off values for this annotation (Exclude extreme values less than or greater than 5%/(number of annotations)/(number of tails excluded))
			if self.tails[j]=="0": #filter both tails
				self.cutoffleftsnp[j] = self.annotmatrixsort[int(0.5+self.rowno2*(0.05/self.annotno/2))] 
				self.cutoffrightsnp[j] = self.annotmatrixsort[int(0.5+self.rowno2*(1.0-0.05/self.annotno/2))] 
			elif self.tails[j]=="-1": #filter left tail only
				self.cutoffleftsnp[j] = self.annotmatrixsort[int(0.5+self.rowno2*(0.05/self.annotno))] 
			elif self.tails[j]=="1": #filter right tail only
				self.cutoffrightsnp[j] = self.annotmatrixsort[int(0.5+self.rowno2*(1.0-0.05/self.annotno))] 
					
			print("Number of records with values and left and right cutoffs for " + str(self.annots[j]) + " for SNPs: " + str(self.rowno2) + ", " + str(self.cutoffleftsnp[j]) + ", " + str(self.cutoffrightsnp[j]) + "\n")
			#Now, find cut-offs for indels
			self.annotmatrix1 = [1000000.0]*self.rownoindel
			self.rowno2 = 0
			while self.rowno2 < self.rownoindel: ## first make new 1D array with column of interest
				self.annotmatrix1[self.rowno2] = self.annotMatrixIndel[j,self.rowno2]
				self.rowno2 = self.rowno2 + 1
			self.annotmatrix1.sort()
			self.annotmatrixsort = self.annotmatrix1

			self.rowno2 = 0
			while self.rowno2 < self.rownoindel:
				if self.annotmatrixsort[self.rowno2] > 999999:
					break
				self.rowno2  = self.rowno2 + 1 
			#Only find cut-offs if there are >100 values to use
			if self.tails[j] == "0" and self.rowno2 > 100: #filter both tails
				self.cutoffleftindel[j] = self.annotmatrixsort[int(0.499+self.rowno2*(0.05/self.annotno/2))-1]; 
				self.cutoffrightindel[j] = self.annotmatrixsort[int(0.499+self.rowno2*(1.0-(0.05/self.annotno/2)))-1]; 
			elif self.tails[j]=="-1" and self.rowno2 > 100: #filter left tail only
				self.cutoffleftindel[j] = self.annotmatrixsort[int(0.499+self.rowno2*(0.05/self.annotno))-1]
			elif self.tails[j]=="1" and self.rowno2 > 100: #filter right tail only
				self.cutoffrightindel[j] = self.annotmatrixsort[int(0.499+self.rowno2*(1.0-(0.05/self.annotno)))-1] 
			print("Number of records with values and left and right cutoffs for " + str(self.annots[j]) + " for indels: " + str(self.rowno2) +  ", " + str(self.cutoffleftindel[j]) + ", " + str(self.cutoffrightindel[j]) + "\n")

	##finally, create a bed file using the vcf records that fall outside the cut-offs or are already filtered

	def run_classifier_sensitive_vcf(self):
		#skip header lines in vcf file containing all sites
		self.x = 0
		self.addfilt = 1
		self.sensvcf_file_lines = self.sensvcf_file.readlines() 
		self.sensvcf_file_idx = 0
		for line in self.sensvcf_file_lines:
			self.sensvcf_file_idx = self.sensvcf_file_idx + 1
			if not(re.search("^\#\#", line)):
				self.x = 1 #skip header lines, but including line starting with #
			if self.addfilt == 1 and re.search("^\#\#INFO", line): #add new filter in header right before INFO header lines
				self.output_filtered_vcf_file.write("##FILTER=<ID=outlier,Description=\"Filtered because one or more annotations were determined to be outliers using VcfOneClassFiltering_v3.1.pl\">\n")
				self.addfilt = 0
			self.output_filtered_vcf_file.write(line)
			if self.x == 1:
				break

		self.bedchrom = ""
		self.bedstart = -1
		self.bedend = -1
		self.filtcnt = 0
		self.totalcnt = 0
		self.filtindivcntsnp = [0]*self.annotno
		self.filtindivcntindel = [0]*self.annotno

		for line_idx in range(self.sensvcf_file_idx, len(self.sensvcf_file_lines)):    
			line = self.sensvcf_file_lines[line_idx]
			# Split up the line into an array
			self.fields = line.split("\t")

			self.nonsnp=0
			self.nogeno=0 #don't output records with no genotype or partial genotypes
			self.ref = self.fields[3]
			self.alt = self.fields[4]
			self.alts = self.alt.split(",")
			
			#is this a SNP or not
			if len(self.ref) > 1:
				self.nonsnp = 1
			if self.nonsnp==0:
				for j in range(0,len(self.alts)):
					if len(self.alts[j]) > 1: 
						self.nonsnp = 1
						break
			

			self.filt = self.fields[6]
			self.info = self.fields[7]
			self.format_field = self.fields[8]
			self.formats = self.format_field.split(":")
			self.chars = self.fields[9].split(":")

			self.filtline=0 #flag to set to 1 if line should be filtered
			if not(self.filt == "PASS" or self.filt == "."):
				self.filtline = 1
			else:
				for j in range(0,self.annotno):
					if self.infoformat[j]=="0": #field is in INFO
						self.infoflds = self.info.split(";")
						for k in range(0,len(self.infoflds)):
							self.infofld = self.infoflds[k].split("=")
							if len(self.infofld)==2 and self.infofld[0] == self.annots[j]: 
								self.infofld1 = self.infofld[1].split(",") #only use first value if field is multivalued (e.g., multiallelic GATK sites)
								self.annotval = self.infofld1[0]
								if(len(self.infofld1) > 1 and float(self.annotval) > float(self.infofld1[1])):
									self.annotval = self.infofld1[1]
								if(len(self.infofld1) > 2 and float(self.annotval) > float(self.infofld1[2])):
									self.annotval = self.infofld1[2]
								if (self.nonsnp==0 and self.annsnp[j]=="1" and (float(self.annotval) < float(self.cutoffleftsnp[j]) or float(self.annotval) > float(self.cutoffrightsnp[j]))): 
									self.filtline = 1 
									self.filtindivcntsnp[j] = self.filtindivcntsnp[j] + 1
								if (self.nonsnp==1 and self.annindel[j]=="1" and (float(self.annotval) < float(self.cutoffleftindel[j]) or float(self.annotval) > float(self.cutoffrightindel[j]))): 
									self.filtline = 1 
									self.filtindivcntindel[j] = self.filtindivcntindel[j] + 1
					else: #field is in FORMAT
						for k in range(0,len(self.chars)):
							if self.formats[k] == self.annots[j]:
								self.charfld1 = self.chars[k].split(",") #use minimum of first 3 values if field is multivalued (e.g., CGA_CEHQ or multiallelic GATK sites)
								self.charval = 0
								if re.search("(.*)\n", self.charfld1[0]):
									self.charval = re.split("(.*)\n", self.charfld1[0])[1]
								else:
									self.charval = self.charfld1[0]
									if len(self.charfld1) > 1:
										if re.search("(.*)\n", self.charfld1[1]): 
											if float(self.charval) > float(re.split("(.*)\n", self.charfld1[1])[1]):
												self.charval = re.split("(.*)\n", self.charfld1[1])[1]
										else: 
											if float(self.charval) > float(self.charfld1[1]):
												self.charval = self.charfld1[1]
										
										if len(self.charfld1) > 2:
											if re.search("(.*)\n", self.charfld1[2]):
												if float(self.charval) > float(re.split("(.*)\n", self.charfld1[2])[1]): 
													self.charval = re.split("(.*)\n", self.charfld1[2])[1]
											else: 
												if float(self.charval) > float(self.charfld1[2]):
													self.charval = self.charfld1[2]
								##print "insidetest:$annotmatrixsnp[$j][$rowno];$rowno\n";
								if (self.nonsnp==0 and self.annsnp[j]=="1" and (not(self.charval == ".") and (float(self.charval) < float(self.cutoffleftsnp[j]) or float(self.charval) > float(self.cutoffrightsnp[j])))): 
									self.filtline = 1 
									self.filtindivcntsnp[j] = self.filtindivcntsnp[j] + 1
								if (self.nonsnp==1 and self.annindel[j]=="1" and (not(self.charval == ".") and (float(self.charval) < float(self.cutoffleftindel[j]) or float(self.charval) > float(self.cutoffrightindel[j])))): 
									self.filtline = 1
									self.filtindivcntindel[j] = self.filtindivcntindel[j] + 1
							
				#also check if GQ<20
				for k in range(0,len(self.chars)):
					if self.formats[k] == "GQ":
						self.charval = 0
						self.charfld1 = self.chars[k].split(",") #only use first value if field is multivalued (e.g., multiallelic GATK sites)
						if re.search("(.*)\n", self.charfld1[0]): 
							self.charval = re.split("(.*)\n", self.charfld1[0])[1]
						else: 
							self.charval = self.charfld1[0]
						##print "insidetest:$annotmatrixsnp[$j][$rowno];$rowno\n";
						if (not(self.charval == ".") and (float(self.charval)<20)):
							self.filtline = 1
					if self.formats[k] == "GT":
						if re.search("\.", self.chars[k]): 
							self.nogeno = 1
			#print("fields: " + str(self.fields))
			if self.filtline == 1: 
				if self.bedend == -1: 
					self.bedchrom = self.fields[0]
					self.bedstart = int(self.fields[1]) - 1 - 50 #subtract 1 to create 0-based bed and add 50bp padding
					self.bedend = int(self.fields[1]) + len(self.fields[3]) + 50 #add number of reference bases to uncertain region + 50bp padding			
				elif ((int(self.fields[1])-1-50) > self.bedend) or not(self.fields[0] == self.bedchrom): ##print previous record if far from current row 
					self.output_filtered_bed_file.write(str(self.bedchrom) + "\t" + str(self.bedstart) + "\t" + str(self.bedend) + "\t" + "CS_" + str(self.callsetName) + "_filt\n")
					self.bedchrom = self.fields[0]
					self.bedstart = int(self.fields[1]) - 1 - 50 #subtract 1 to create 0-based bed and add 50bp padding
					self.bedend = int(self.fields[1]) + len(self.fields[3]) + 50 #add number of reference bases to uncertain region + 50bp padding
				else: #extend previous record
					if (self.bedend < (int(self.fields[1]) + len(self.fields[3]) + 50)): 
						self.bedend = int(self.fields[1]) + len(self.fields[3]) + 50 #add number of reference bases to uncertain region + 50bp padding
				self.filtcnt = self.filtcnt + 1
				self.filt="outlier"
			
			if (self.nogeno==0):
				self.output_filtered_vcf_file.write(str(self.fields[0]) + "\t" + str(self.fields[1]) + "\t" + str(self.fields[2]) + "\t" + str(self.fields[3]) + "\t" + str(self.fields[4]) + "\t" + str(self.fields[5]) + "\t" + str(self.filt) + "\t" + str(self.fields[7]) + "\t" + str(self.fields[8]) + "\t" + str(self.fields[9]) + "\n")
				self.totalcnt = self.totalcnt + 1

		self.output_filtered_bed_file.write(str(self.bedchrom) + "\t" + str(self.bedstart) + "\t" + str(self.bedend) + "\t" + "CS_" + str(self.callsetName) + "_filt\n") ##print last record

		print("Number of vcf records filtered: " + str(self.filtcnt) + "out of " + str(self.totalcnt) + "\n")
		for j in range(0, self.annotno):
			print("Number of snps filtered due to outlier " + str(self.annots[j]) + ": " + str(self.filtindivcntsnp[j]) + "\n")
			print("Number of non-snps filtered due to outlier " + str(self.annots[j]) + ": " + str(self.filtindivcntindel[j]) + "\n")

		self.sensvcf_file.close()
		self.output_filtered_bed_file.close()
		self.output_filtered_vcf_file.close()
		self.annotationsToFilter_file.close()

	def one_class_filtering(self):
		self.parse_annotations()
		self.train_classifier()
		self.run_classifier_sensitive_vcf()
