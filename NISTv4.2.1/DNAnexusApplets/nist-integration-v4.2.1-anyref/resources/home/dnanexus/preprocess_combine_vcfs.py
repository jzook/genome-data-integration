import argparse
import subprocess
import os
import glob
import re

class PreProcessCombine(object):
	parser = None
	args = None
	path = None
	cstable_file = None
	cstable_file_lines = None
	cstable_file_header = None
	output_file_lines = None
	m = None
	vcfs = None
	cstable_line_head = None
	vcf_no_string = None
	line = None
	line_data = None
	vcf_prefix = None
	vcfin_path = None
	vcfin_file = None
			
	outindiv_path = None
	outindiv_file = None
	
	outdupbed_path = None
	outdupbed_file = None
	
	vcfin_lines = None
	duplicate_position = None
	anyduplicate = None
	vcf_previous_line = None
	line_data_previous = None
	split_vcf_line = None

	duplicate_position = None
	anyduplicate = None
	bed1 = None
	bed2 = None
	temp_callable_loci_bed_file = None
	temp_callable_loci_cols = None
	bedname = None

	def __init__(self, cstable, p):
		# self.parser = argparse.ArgumentParser(description='Process and combine VCF files.')

		# self.parser.add_argument('-cstable', metavar='txt', help='The input Tab-delimited text file containing columns Platform, Dataset, Callset, vcfAll.vcf.gz, callableBed.bed, and annotationsFile')
		# self.parser.add_argument('-p', metavar='FILEPATH', help='Path where vcf and bed files are located')

		
		self.path = p
		self.cstable_file = open(cstable, "r")

		self.cstable_file_lines = self.cstable_file.readlines()
		self.cstable_file_header = None
		self.output_file_lines = None
		self.m = None
		self.vcfs = None
		self.cstable_line_head = None
		self.vcf_no_string = None
		self.vcf_prefix = None
		self.vcfin_path = None
		self.vcfin_file = None
		
		self.outindiv_path = None
		self.outindiv_file = None
				
		self.outdupbed_path = None
		self.outdupbed_file = None
				
		self.vcfin_lines = None
		self.duplicate_position = None
		self.anyduplicate = None
		self.vcf_previous_line = None
		self.line_data_previous = None
		self.split_vcf_line = None

		self.duplicate_position = None
		self.anyduplicate = None
		self.bed1 = None
		self.bed2 = None
		self.temp_callable_loci_bed_file = None
		self.temp_callable_loci_cols = None
		self.bedname = None

	def process_callset_table(self):
		self.cstable_file_header = ""
		self.output_file_lines = []
		self.m = 1
		self.vcfs = ""
		self.cstable_line_head = []
		for line in self.cstable_file_lines:
			
			self.vcf_no_string = str(self.m)
			if self.m < 10:
				self.vcf_no_string = "0" + str(self.m)
		
			line = line.rstrip("\n")
			line_data = line.split("\t")
		
			if re.search("^\#", line):
				self.cstable_file_header = line
				self.cstable_line_head = self.cstable_file_header.split("\t")
				continue	

			if re.search("(.*)\.vcf.gz", line_data[3]):
				self.vcf_prefix = re.split("(.*)\.vcf.gz", line_data[3])[1]
			else:
				#print("vcfAll does not end in .vcf.gz")
				exit(1)
			#print("gunzip -c " + self.path + line_data[3] + " | cut -f1-10 | grep -v '0/0' | grep -v '0|0' | vcfallelicprimitives -k -g | perl vcfsorter.pl genome.dict - > " + self.path + self.vcf_prefix + "_nohomref.vcf")
			subprocess.check_call("gunzip -c " + self.path + line_data[3] + " | cut -f1-10 | grep -v '0/0' | grep -v '0|0' | vcfallelicprimitives -k -g | perl vcfsorter.pl genome.dict - > " + self.path + self.vcf_prefix + "_nohomref.vcf", shell=True)
			
			self.vcfin_path = self.path + self.vcf_prefix + "_nohomref.vcf"
			#print(self.vcfin_path)
			self.vcfin_file = open(self.vcfin_path, "r")
			
			self.outindiv_path = self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf"
			#print(self.outindiv_path)
			self.outindiv_file = open(self.outindiv_path, "w+")
			
			self.outdupbed_path = self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_dup.bed"
			#print(self.outdupbed_path)
			self.outdupbed_file = open(self.outdupbed_path, "w+")
			
			self.vcfin_lines = self.vcfin_file.readlines()
			# #print(self.vcfin_lines)
			self.duplicate_position = 0
			self.anyduplicate = 0 
			self.vcf_previous_line = self.vcfin_lines[0].rstrip("\n")
			self.line_data_previous = self.vcf_previous_line.split("\t")
			
			for i in range(1,len(self.vcfin_lines)):
				self.vcfin_line = self.vcfin_lines[i].rstrip("\n")
				self.split_vcf_line = self.vcfin_line.split("\t")
				if self.vcf_previous_line.startswith('#CHROM'):
					self.outindiv_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVCF" + self.vcf_no_string + "_"+ self.split_vcf_line[2] + "\n")
				elif self.vcf_previous_line.startswith('#'):
					self.outindiv_file.write(self.vcf_previous_line + "\n")
				elif self.line_data_previous[1] == self.split_vcf_line[1]:
					self.duplicate_position = 1
					self.anyduplicate = 1
					self.bed1 = int(self.split_vcf_line[1]) - 50
					self.bed2 = int(self.split_vcf_line[1]) + 50
					#print("Duplicated vcf lines: " + self.vcf_previous_line + "\n" + self.vcfin_line + self.splint_vcf_line[0] + "\t" + str(self.bed1) + "\t" + str(self.bed2) + "\n")
					self.outdupbed_file.write(self.split_vcf_line[0]+"\t"+str(self.bed1)+"\t"+str(self.bed2)+"\n")
				elif self.duplicate_position == 1:
					self.duplicate_position = 0
				else:
					self.outindiv_file.write(self.vcf_previous_line + "\n")

				self.vcf_previous_line = self.vcfin_line
				self.line_data_previous = self.vcfin_line.split("\t")

			self.outindiv_file.write(self.vcf_previous_line + "\n")
			self.vcfin_file.close()
			self.outindiv_file.close() 
			self.outdupbed_file.close() 

			subprocess.check_call("bgzip -c " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf > " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf.gz", shell=True) 
			subprocess.check_call("tabix -p vcf " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf.gz", shell=True)

			self.vcfs = self.vcfs + " " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf.gz"	
			#print("any duplicate: " + str(self.anyduplicate))
			if not(line_data[4] == "none"):
				self.temp_callable_loci_bed_file = open(self.path+line_data[4], "r").readlines()
				self.temp_callable_loci_cols = self.temp_callable_loci_bed_file[0].split("\t")
				if len(self.temp_callable_loci_cols) > 3: 
					if self.anyduplicate == 1:
						subprocess.check_call("grep CALLABLE " + self.path + line_data[4] + " | cut -f1-3 | subtractBed -a stdin -b " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_dup.bed > " + self.path + self.vcf_prefix + "_callable_tmp_1.bed", shell=True) 
					else: 
						subprocess.check_call("grep CALLABLE " + self.path + line_data[4] + " > " + self.path + self.vcf_prefix + "_callable_tmp_1.bed", shell=True)
				else:
					if self.anyduplicate == 1:
						subprocess.check_call("cut -f1-3 " + self.path + line_data[4] + " | subtractBed -a stdin -b " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_dup.bed > " + self.path + self.vcf_prefix + "_callable_tmp_1.bed", shell=True)
					else:
						subprocess.check_call("cp " + self.path + line_data[4] + " " + self.path + self.vcf_prefix + "_callable_tmp_1.bed",shell=True)
				
				c = 6
				d = 1
				while c < len(self.cstable_line_head):
					if line_data[c] == "1":
						subprocess.check_call("subtractBed -a " + self.path + self.vcf_prefix + "_callable_tmp_" + str(d) + ".bed -b " + self.cstable_line_head[c] + " > " + self.path + self.vcf_prefix + "_callable_tmp_" + str(c) + ".bed", shell=True)
						d = c
					c = c + 1
				subprocess.check_call("awk \'BEGIN {FS = OFS = \"\t\"} {print $1 \"\t\" $2 \"\t\" $3 \"\t\" \"CS_"+ line_data[2] + "_callable\"}\' " + self.path + self.vcf_prefix + "_callable_tmp_" + str(d) + ".bed > " + self.path + self.vcf_prefix + "_callable_processed.bed", shell=True) 
			self.m = self.m + 1

		#print(self.vcfs + "\n")

		self.cstable_file.close()

	def add_difficult_annotation(self):
		c = 6
		while(c < len(self.cstable_line_head)):
			bedname = self.cstable_line_head[c]
			if re.search("(.*).bed.gz", self.cstable_line_head[c]):
				bedname = re.split("(.*).bed.gz", self.cstable_line_head[c])[1]
				#print(bedname)
				#print("gunzip -c " + self.cstable_line_head[c] + " | awk \'BEGIN {FS = OFS = \"\t\"} {#print $1 \"\t\" $2 \"\t\" $3 \"\t\" \"" + bedname + "\"}\' > " + bedname + "_diffbedannotcol.bed")
				subprocess.check_call("gunzip -c " + self.cstable_line_head[c] + " | awk \'BEGIN {FS = OFS = \"\t\"} {print $1 \"\t\" $2  \"\t\" $3 \"\t\" \"" + bedname + "\"}\' > " + bedname + "_diffbedannotcol.bed", shell=True)
			elif re.search("(.*).bed", self.cstable_line_head[c]):
				bedname = re.split("(.*).bed", self.cstable_line_head[c])[1]
				#print(bedname)
				#print("awk \'BEGIN {FS = OFS = \"\t\"} {#print $1 \"\t\" $2 \"\t\" $3 \"\t\" \"" + bedname + "\"}\' " + self.cstable_line_head[c] + " > " + bedname + "_diffbedannotcol.bed")
				subprocess.check_call("awk \'BEGIN {FS = OFS = \"\t\"} {print $1 \"\t\" $2 \"\t\" $3 \"\t\" \"" + bedname + "\"}\' " + self.cstable_line_head[c] + " > " + bedname + "_diffbedannotcol.bed", shell=True) 
			c = c + 1

	def preprocess(self):
		self.process_callset_table()
		self.add_difficult_annotation()
