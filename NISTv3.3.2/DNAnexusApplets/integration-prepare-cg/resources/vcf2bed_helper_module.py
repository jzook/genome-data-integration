#!/usr/bin/env python


import sys, re

PRINT_NO_CALLS = 'print_no-calls'
NOT_PASS_CALLS_ARE_NO_CALLS = 'not_pass_calls_are_no-calls'
USE_FT_INSTEAD_OF_FILTER = 'use_FT_instead_of_FILTER'
INCLUDE_PARTIAL_CALLS = 'include_partial_calls'

OPTIONS_LIST = [PRINT_NO_CALLS, NOT_PASS_CALLS_ARE_NO_CALLS,
			USE_FT_INSTEAD_OF_FILTER, INCLUDE_PARTIAL_CALLS]

END = "END"
PASS = "PASS"
NC_CORR = "NC-CORR"
NOCALL = "NOCALL"

FT = "FT"
GT = "GT"

EMPTY_STRING = ""

CGA_MEDEL = "CGA_MEDEL"
MATEID = "MATEID"
SVTYPE = "SVTYPE"
CGA_BNDG = "CGA_BNDG"
CGA_BNDGO = "CGA_BNDGO"
CIPOS = "CIPOS"
IMPRECISE = "IMPRECISE"
MEINFO = "MEINFO"
SVLEN = "SVLEN"

FULL_REF_CALL = 'full_ref_call'
FULL_VARIANT_CALL = 'full_variant_call'
FULL_NO_CALL = 'full_no_call'
PARTIAL_VARIANT_CALL = 'partial_variant_call'
PARTIAL_REF_CALL = 'partial_ref_call'
MIXED_REF_VARIANT_CALL = 'mixed_ref_variant_call'

## NOTE-- I was originally planning to screen out SVTYPE along with the other suspect INFO keys,
##        but it is already screened out in a separate statement (which occurs earlier in the code)
SUSPECT_INFO_KEYS = [CGA_MEDEL, MATEID, SVTYPE, CGA_BNDG, CGA_BNDGO, CIPOS, IMPRECISE, MEINFO, SVLEN]

CURR_LINE = ""

MIN_NUM_FIELDS_IN_ARGV = 4
MAX_NUM_FIELDS_IN_ARGV = MIN_NUM_FIELDS_IN_ARGV + len(OPTIONS_LIST)

##############################################################################################
## BEFORE MAIN LOOP - INITIALIZE VARIABLES, GET USER OPTIONS, AND PARSE NO-REF REGION FILE: ##
##############################################################################################

def get_user_options_dict(user_options_list):
	user_options_dict = {}
	for option in OPTIONS_LIST:
		user_options_dict[option] = False

	for user_option in user_options_list:
		if user_option[0] != '-':
			sys.exit("User options (3rd argument and on) must each begin with a '-'")
		user_option = user_option[1:]
		if not (user_option in OPTIONS_LIST):
			print "Warning: unrecognized option \"%s\"" % user_option
		else:
			user_options_dict[user_option] = True
	
	return user_options_dict

def get_chrom_to_length_dict_from_header_lines_of_vcf_file(vcf_input_filename):
	chrom_to_length_dict = {}

	vcf_input_file = open(vcf_input_filename, 'r')
	for line in vcf_input_file:
		line = line.rstrip('\n\r')
		if (len(line) != 0) and (line[0] == '#'):
			check_vcf_line_for_contig_information_and_update_chrom_length_dict(line, chrom_to_length_dict)
		elif (len(line) != 0) and (line[0] != '#'):
			break
	vcf_input_file.close()

	return chrom_to_length_dict

def get_chrom_to_no_ref_regions_dict(no_ref_regions_input_filename):
	chrom_to_no_ref_regions_dict = {}

	no_ref_regions_input_file = open(no_ref_regions_input_filename, 'r')
	for line in no_ref_regions_input_file:
		line = line.rstrip('\n\r')
		field_list = line.split('\t')
		chrom = field_list[0].replace('chr', '')
		start_pos_0bho = int(field_list[1])
		end_pos_0bho = int(field_list[2])

		## convert from 0-based half-open to 1-based closed
		start_pos_1bc = start_pos_0bho + 1
		end_pos_1bc = end_pos_0bho

		## check for zero or negative-length regions
		if start_pos_1bc > end_pos_1bc:
			sys.exit("Unexpected: zero or negative-length 'no-ref' region detected")

		## update dictionary with region
		if not (chrom in chrom_to_no_ref_regions_dict):
			chrom_to_no_ref_regions_dict[chrom] = []
		chrom_to_no_ref_regions_dict[chrom].append([start_pos_1bc, end_pos_1bc])
	no_ref_regions_input_file.close()

	return chrom_to_no_ref_regions_dict

def initialize_variables_for_main_loop():

	# dictionaries
	chrom_to_length_dict = {}
	chrom_to_highest_observed_end_position_of_a_feature_dict = {}
	call_type_to_values_of_filter_to_use_dict = {}

	# no previous feature
	prev_feature_chrom = -1
	prev_feature_start_pos = 'NA'
	prev_feature_end_pos = 'NA'
	prev_feature_is_call = 'NA'
	prev_feature_call_type = -1
	prev_feature_is_padding_base = 'NA'
	prev_feature_filter = 'NA'

	prev_feature = [prev_feature_chrom, prev_feature_start_pos, prev_feature_end_pos, prev_feature_is_call, prev_feature_call_type, prev_feature_is_padding_base, prev_feature_filter]

	# regions of interest will be either call regions or no-call regions,
	# depending on the value of the user-specified option "-print_no-calls"
	#
	# the value of this variable will be initialized later on in the code
	region_start_pos = 'NA'

	return prev_feature, chrom_to_length_dict, chrom_to_highest_observed_end_position_of_a_feature_dict, \
			call_type_to_values_of_filter_to_use_dict, region_start_pos

######################################################
## MAIN LOOP - PREPROCESS LINE AND EXTRACT ENTRIES: ##
######################################################

def check_vcf_line_for_contig_information_and_update_chrom_length_dict(line, chrom_to_length_dict):
	if line[0:8] == '##contig':
		contig_field_list = line[8:].strip('=').strip('><').split(',')	

		cf_key_to_value_dict = {}
		for contig_field in contig_field_list:
			[cf_key, cf_value] = contig_field.split("=")
			cf_key_to_value_dict[cf_key] = cf_value
		chrom_ID = cf_key_to_value_dict['ID'].strip()
		chrom_length = int(cf_key_to_value_dict['length'])
		chrom_to_length_dict[chrom_ID] = chrom_length
	return

def get_base_padding_offset(ref_string, alt_string):
	if ',' in ref_string:
		sys.exit("did not expect more than one value for the reference in a single line of the vcf file")

	if ((ref_string == "") and (alt_string == "")) or (alt_string == '<CGA_NOCALL>'):
		# leave base_padding_offset as 0 for certain special cases 
		base_padding_offset = 0
	else:
		alt_list = alt_string.split(',')

		# determine base padding offset
		base_padding_offset = 0
		break_out_of_loop = False
		while True:
			for alt_entry in alt_list:
				# plan to break if the index specified by base_padding_offset is off the end of either ref_string or alt_entry
				if (base_padding_offset == len(ref_string)) or (base_padding_offset == len(alt_entry)):
					break_out_of_loop = True
				else:
					# convert lowercase to uppercase
					refbase = ref_string[base_padding_offset].upper()
					altbase = alt_entry[base_padding_offset].upper()
					# plan to break if refbase and altbase do not match or if they are not one of 'A', 'C', 'G', or 'T'
					ACGT_array = ['A', 'C', 'G', 'T']
					if ((refbase != altbase) or (not (refbase in ACGT_array))):
						break_out_of_loop = True

				if break_out_of_loop:
					break
			if break_out_of_loop:
				break

			base_padding_offset += 1

	return base_padding_offset

def get_end_position(preliminary_start_pos, ref_string, info_key_to_value_dict):

	## if there is an end tag, then use it.
	## Otherwise, infer the end position from the start position and the length of the reference string.
	if END in info_key_to_value_dict:
		end_pos = int(info_key_to_value_dict[END])
	else:
		end_pos = preliminary_start_pos + len(ref_string) - 1 

	return end_pos

def decide_whether_feature_is_a_call_or_no_call(filter_string, format_key_to_value_dict, user_options_dict,
						call_type_to_values_of_filter_to_use_dict):

	## read GT flag to get allelic call information
	at_least_one_no_call, at_least_one_variant_call, at_least_one_reference_call = \
			read_genotype_information(format_key_to_value_dict)
	
	## classify the feature into one of six possible groups by its allelic call information
	full_ref_call, full_variant_call, full_no_call, partial_variant_call, partial_ref_call, mixed_ref_variant_call, call_type = \
			evaluate_allelic_calls(at_least_one_no_call, at_least_one_variant_call, at_least_one_reference_call)

	## Get a value that may or may not be used to decide whether a given candidate is a call or a no-call
	## (depending on the user-specified option "-not_pass_calls_are_no-calls").
	##
	## This value may come from the FILTER or the FT flag
	## (depending on the user-specified option "-use_FT_instead_of_FILTER").
	filter_to_use = get_filter_to_use(filter_string, format_key_to_value_dict, user_options_dict)

	## record all observed filter values for each type of call
	update_dictionary_of_all_observed_filter_values_for_each_call_type(call_type_to_values_of_filter_to_use_dict, call_type, filter_to_use)

	## determine whether we consider this a call or a no-call
	is_call = apply_rules_to_decide_call_or_no_call(full_ref_call, full_variant_call, full_no_call,
									partial_variant_call, partial_ref_call,
									mixed_ref_variant_call, filter_to_use, user_options_dict)

	return is_call, call_type, filter_to_use

def read_genotype_information(format_key_to_value_dict):
	## read GT (genotype) information
	GT_string = format_key_to_value_dict[GT]
	allele_list = re.split("\||/", GT_string)
	
	if len(allele_list) == 0:
		sys.exit("length of allele list is 0")

	at_least_one_no_call = False
	at_least_one_variant_call = False
	at_least_one_reference_call = False
	for allele in allele_list:
		if allele == '.':
			## NOCALL
			at_least_one_no_call = True
		elif allele == '0':
			## ref call
			at_least_one_reference_call = True
		elif re.match(r'[0-9]+', allele):
			## variant call
			at_least_one_variant_call = True
		else:
			sys.exit("Specified value of GT does not match specification: '%s'" % GT_string)

	return at_least_one_no_call, at_least_one_variant_call, at_least_one_reference_call

def evaluate_allelic_calls(at_least_one_no_call, at_least_one_variant_call, at_least_one_reference_call):

	condition_count = 0
	call_type = -1

	full_ref_call, full_variant_call, full_no_call, condition_count, call_type = \
			 test_for_full_calls_or_no_calls(at_least_one_no_call, at_least_one_variant_call,
							at_least_one_reference_call, condition_count, call_type)


	partial_variant_call, partial_ref_call, mixed_ref_variant_call, condition_count, call_type = \
			 test_for_partial_calls_and_mixed_ref_variant_calls(at_least_one_no_call, at_least_one_variant_call,
									at_least_one_reference_call, condition_count, call_type)



	## check for unexpected combinations of allele calls in GT results
	if at_least_one_reference_call and at_least_one_variant_call and at_least_one_no_call:
		sys.exit("Unexpected: saw alleles for a ref call, a variant call, and a no-call in the same line of the vcf file")

	if (not at_least_one_reference_call) and (not at_least_one_variant_call) and (not at_least_one_no_call):
		sys.exit("Unexpected: did not see any alleles for a ref call, a variant call, or a no-call in the current line of the vcf file")

	## check my own coding logic
	if condition_count != 1:
		sys.exit("Internal error: Exactly one of the expected conditions should be the case!")

	if call_type == -1:
		sys.exit("Internal error: Call type should not be -1")

	return full_ref_call, full_variant_call, full_no_call, partial_variant_call, partial_ref_call, mixed_ref_variant_call, call_type

def get_filter_to_use(filter_string, format_key_to_value_dict, user_options_dict):
	manually_setting_filter_to_use = False
	if user_options_dict[USE_FT_INSTEAD_OF_FILTER]:
		if FT in format_key_to_value_dict:
			filter_to_use = format_key_to_value_dict[FT]
		else:
			filter_to_use = EMPTY_STRING
			manually_setting_filter_to_use = True
	else:
		filter_to_use = filter_string
	if (not manually_setting_filter_to_use) and (filter_to_use == EMPTY_STRING):
		sys.exit("Unexpected: FT or FILTER was set to the empty string in the vcf file")

	return filter_to_use

def update_dictionary_of_all_observed_filter_values_for_each_call_type(call_type_to_values_of_filter_to_use_dict, call_type, filter_to_use):
	if not (call_type in call_type_to_values_of_filter_to_use_dict):
		call_type_to_values_of_filter_to_use_dict[call_type] = {}
	if not (filter_to_use in call_type_to_values_of_filter_to_use_dict[call_type]):
		call_type_to_values_of_filter_to_use_dict[call_type][filter_to_use] = 0
	call_type_to_values_of_filter_to_use_dict[call_type][filter_to_use] += 1

	return

def apply_rules_to_decide_call_or_no_call(full_ref_call, full_variant_call, full_no_call,
								partial_variant_call, partial_ref_call,
								mixed_ref_variant_call, filter_to_use, user_options_dict):

	if (full_ref_call or full_variant_call or mixed_ref_variant_call or
	    (user_options_dict[INCLUDE_PARTIAL_CALLS] and (partial_ref_call or partial_variant_call))):

		pass_strings_list = [PASS]

		if full_ref_call:
			# consider a full ref call to be high-confidence ("PASS")
			# in the case where the associated line of the vcf file has:
			#     (A) an empty entry for the filter to use [if the option -use_FT_instead_of_FILTER was specified]
			#     (B) an entry of "." for the filter to use [if the option -use_FT_instead_of_FILTER was not specified]
			if user_options_dict[USE_FT_INSTEAD_OF_FILTER]:
				pass_strings_list.append(EMPTY_STRING)
			else:
				pass_strings_list.append(".")

		if partial_ref_call:
			# consider a partial ref call to be high-confidence ("PASS")
			# in the case where the associated line of the vcf file has:
			#     (A) an empty entry for the filter to use [if the option -use_FT_instead_of_FILTER was specified]
			#     (B) an entry of "NC-CORR" or "NOCALL" for the filter to use [if the option -use_FT_instead_of_FILTER was not specified]
			if user_options_dict[USE_FT_INSTEAD_OF_FILTER]:
				pass_strings_list.append(EMPTY_STRING)
			else:
				pass_strings_list.append(NC_CORR)
				pass_strings_list.append(NOCALL)

		if (filter_to_use in pass_strings_list) or (not user_options_dict[NOT_PASS_CALLS_ARE_NO_CALLS]):
			## passed filter or filter doesn't matter
			is_call = True
		else:
			## did not pass filter and filter matters
			is_call = False
	else:
		## consider a full no-call to always be a no-call
		## no matter the value of its FT/FILTER flag
		## and even if it does not have an FT flag
		##
		## consider the same to be true of partial calls
		## (which are assumed to include mixed ref/variant calls),
		## in the case where the user has not supplied
		## the option "-include_partial_calls")
		is_call = False
	
	if not (full_ref_call or full_variant_call or full_no_call or partial_ref_call or partial_variant_call or mixed_ref_variant_call):
		sys.exit("Internal error: expected at least one condition to be set to True")

	return is_call

def test_for_full_calls_or_no_calls(at_least_one_no_call, at_least_one_variant_call,
					at_least_one_reference_call, condition_count, call_type):

	if at_least_one_reference_call and (not at_least_one_variant_call) and (not at_least_one_no_call):
		call_type = FULL_REF_CALL
		full_ref_call = True
		condition_count += 1
	else:
		full_ref_call = False

	if at_least_one_variant_call and (not at_least_one_reference_call) and (not at_least_one_no_call):
		call_type = FULL_VARIANT_CALL
		full_variant_call = True
		condition_count += 1
	else:
		full_variant_call = False

	if at_least_one_no_call and (not at_least_one_reference_call) and (not at_least_one_variant_call):
		call_type = FULL_NO_CALL
		full_no_call = True
		condition_count += 1
	else:
		full_no_call = False

	return full_ref_call, full_variant_call, full_no_call, condition_count, call_type

def test_for_partial_calls_and_mixed_ref_variant_calls(at_least_one_no_call, at_least_one_variant_call,
					at_least_one_reference_call, condition_count, call_type):

	## this excludes variant/ref calls:
	if at_least_one_variant_call and at_least_one_no_call and (not at_least_one_reference_call):
		call_type = PARTIAL_VARIANT_CALL
		partial_variant_call = True
		condition_count += 1
	else:
		partial_variant_call = False

	## this excludes variant/ref calls:
	if at_least_one_reference_call and at_least_one_no_call and (not at_least_one_variant_call):
		call_type = PARTIAL_REF_CALL
		partial_ref_call = True
		condition_count += 1
	else:
		partial_ref_call = False

	## heterozygous ref/variant call
	if at_least_one_reference_call and at_least_one_variant_call and (not at_least_one_no_call):
		call_type = MIXED_REF_VARIANT_CALL
		mixed_ref_variant_call = True
		condition_count += 1
	else:
		mixed_ref_variant_call = False

	return partial_variant_call, partial_ref_call, mixed_ref_variant_call, condition_count, call_type

def append_padding_base_feature_to_features_to_consider_list(chrom, start_pos, base_padding_offset, features_to_consider_list, filter_to_use):
	padding_base_start_pos = start_pos - base_padding_offset
	padding_base_end_pos = start_pos - 1
	
	PB_is_padding_base = True
	# always assume padding base is high quality reference
	padding_base_is_call = True
	PB_call_type = FULL_REF_CALL
	
	# append padding base vcf feature:
	#     [chromosome, start position, end position, "is it a call?", call type, "is it a padding base?"]
	features_to_consider_list.append([chrom, padding_base_start_pos, padding_base_end_pos, padding_base_is_call,
						PB_call_type, PB_is_padding_base, filter_to_use])

	return

def append_main_feature_to_features_to_consider_list(chrom, start_pos, end_pos, is_call, call_type, features_to_consider_list, filter_to_use):
	# append main vcf feature:
	#     [chromosome, start position, end position, "is it a call?", call type, "is it a padding base?"]
	MF_is_padding_base = False
	features_to_consider_list.append([chrom, start_pos, end_pos, is_call,
						call_type, MF_is_padding_base, filter_to_use])
	
	return

##################################################################################################
## MAIN LOOP - CHECK ASSUMPTIONS ABOUT CURRENT FEATURE AND ITS RELATIONSHIP TO PREVIOUS FEATURE ##
##################################################################################################

def check_assumptions_about_current_and_previous_feature(prev_feature_chrom, prev_feature_end_pos, chrom, start_pos, end_pos):
	## check three assumptions:

	# if previous feature was on the same chromosome
	if (prev_feature_end_pos != 'NA') and (chrom == prev_feature_chrom):
		# check that the current feature's start is after the previous feature's end
		if start_pos <= prev_feature_end_pos:
			print CURR_LINE
			sys.exit("Unexpected: current feature start pos (%d) <= prev feature end pos (%d)" % (start_pos, prev_feature_end_pos))
		# check that the current feature's end is at or after the previous feature's end
		if end_pos < prev_feature_end_pos:
			print CURR_LINE
			sys.exit("Unexpected: current feature end pos < prev feature end pos")
	# check that the current feature's end is not before the position before its start
	if end_pos < (start_pos - 1):
		print CURR_LINE
		sys.exit("Unexpected: current feature end pos < one less than current feature start pos")
	
	return

######################################################
## MAIN LOOP - UPDATE AND PRINT REGIONS OF INTEREST ##
######################################################

def update_and_print_regions_of_interest(prev_feature_is_call, prev_feature_chrom, prev_feature_end_pos,
					is_call, chrom, start_pos, end_pos, 
					region_start_pos, user_options_dict, chrom_to_length_dict,
					chrom_to_highest_observed_end_position_of_a_feature_dict,
					chrom_to_no_ref_regions_dict, bed_output_file):

	## get convenient terms to use in upcoming conditional statements
	is_adjacent, is_no_call, prev_feature_is_no_call = \
			get_convenient_terms_for_print_regions_logic(prev_feature_is_call, prev_feature_chrom, prev_feature_end_pos,
									is_call, chrom, start_pos)

							

	if user_options_dict[PRINT_NO_CALLS]:
		region_start_pos = update_and_print_no_call_regions(prev_feature_is_call, prev_feature_is_no_call,
										prev_feature_chrom, prev_feature_end_pos,
										is_call, is_no_call, start_pos,
										is_adjacent, region_start_pos,
										chrom_to_no_ref_regions_dict, bed_output_file)

	else:
		## only recording and printing call regions
		if prev_feature_chrom != chrom:
			region_start_pos = update_and_print_call_regions__new_chromosome(
										prev_feature_is_no_call, prev_feature_chrom,
										prev_feature_end_pos,
										chrom, is_no_call, start_pos,
										region_start_pos, chrom_to_length_dict,
										chrom_to_no_ref_regions_dict, bed_output_file)

		else:
			region_start_pos = update_and_print_call_regions__same_chromosome(
							prev_feature_is_call, prev_feature_is_no_call, prev_feature_end_pos,
							chrom, is_call, is_no_call, start_pos, is_adjacent, region_start_pos, 
							chrom_to_no_ref_regions_dict, bed_output_file)

	
		

	## update end position of most recent feature for this chromosome
	if ((not (chrom in chrom_to_highest_observed_end_position_of_a_feature_dict)) or
	    (end_pos >= chrom_to_highest_observed_end_position_of_a_feature_dict[chrom])):
		chrom_to_highest_observed_end_position_of_a_feature_dict[chrom] = end_pos
	else:
		sys.exit('Internal error: chrom %s: end pos of most recent feature (%d) is < highest observed end pos (%d)' % \
			 (chrom, end_pos, chrom_to_highest_observed_end_position_of_a_feature_dict[chrom]))


	return region_start_pos

def get_convenient_terms_for_print_regions_logic(prev_feature_is_call, prev_feature_chrom, prev_feature_end_pos,
						is_call, chrom, start_pos):
								
	# current feature is adjacent to previous feature
	if (prev_feature_chrom == chrom) and (prev_feature_end_pos == (start_pos - 1)):
		is_adjacent = True
	else:
		is_adjacent = False

	# current feature is no-call
	is_no_call = (not is_call)

	# prev feature is no-call
	prev_feature_is_no_call = (not prev_feature_is_call)

	return is_adjacent, is_no_call, prev_feature_is_no_call

def update_and_print_no_call_regions(prev_feature_is_call, prev_feature_is_no_call, prev_feature_chrom, prev_feature_end_pos,
						is_call, is_no_call, start_pos, 
						is_adjacent, region_start_pos, chrom_to_no_ref_regions_dict, bed_output_file):
	## only recording and printing no-call regions
	if (prev_feature_is_no_call and
	    (is_call or
	     (is_no_call and (not is_adjacent)))):
		## report old region
		convert_idx_and_write_region_to_bed_file(prev_feature_chrom, region_start_pos, prev_feature_end_pos,
							chrom_to_no_ref_regions_dict, bed_output_file)
	if (is_no_call and 
	    (prev_feature_is_call or
	     (prev_feature_is_no_call and (not is_adjacent)))):
		## start new region
		region_start_pos = start_pos

	return region_start_pos

def update_and_print_call_regions__new_chromosome(prev_feature_is_no_call, prev_feature_chrom,
								prev_feature_end_pos,
								chrom, is_no_call, start_pos,
								region_start_pos, chrom_to_length_dict,
								chrom_to_no_ref_regions_dict, bed_output_file):

	region_start_pos = update_and_print_call_regions_for_end_of_old_chromosome(prev_feature_is_no_call, prev_feature_chrom, prev_feature_end_pos,
										region_start_pos, chrom_to_length_dict,
										chrom_to_no_ref_regions_dict, bed_output_file)

	region_start_pos = update_and_print_call_regions_for_start_of_new_chromosome(chrom, is_no_call, start_pos, region_start_pos, 
											chrom_to_no_ref_regions_dict, bed_output_file)

	return region_start_pos

def update_and_print_call_regions_for_end_of_old_chromosome(prev_feature_is_no_call, prev_feature_chrom, prev_feature_end_pos,
								region_start_pos, chrom_to_length_dict,
								chrom_to_no_ref_regions_dict, bed_output_file):
	## if there is an old chromosome
	if prev_feature_chrom != -1:
		## handle end of old chromosome
		if prev_feature_is_no_call:
			# create new region
			region_start_pos = prev_feature_end_pos + 1
		if region_start_pos <= chrom_to_length_dict[prev_feature_chrom]:
			# report new region
			convert_idx_and_write_region_to_bed_file(prev_feature_chrom, region_start_pos,
								chrom_to_length_dict[prev_feature_chrom],
								chrom_to_no_ref_regions_dict, bed_output_file)

	return region_start_pos

def update_and_print_call_regions_for_start_of_new_chromosome(chrom, is_no_call, start_pos, region_start_pos, 
								chrom_to_no_ref_regions_dict, bed_output_file):
	## handle beginning of new chromosome
	# create new region
	region_start_pos = 1
	if is_no_call:
		if region_start_pos < start_pos:
			# report new region
			convert_idx_and_write_region_to_bed_file(chrom, region_start_pos, (start_pos - 1),
								chrom_to_no_ref_regions_dict, bed_output_file)

	return region_start_pos

def update_and_print_call_regions__same_chromosome(prev_feature_is_call, prev_feature_is_no_call, prev_feature_end_pos,
								chrom, is_call, is_no_call, start_pos,
								is_adjacent, region_start_pos,
								chrom_to_no_ref_regions_dict, bed_output_file):

	## if this feature and the previous feature are on the same chromosome
	if is_no_call and prev_feature_is_call:
		# report old region
		if region_start_pos > start_pos:
			sys.exit("Unexpected: region to print start pos > current feature start pos")
		convert_idx_and_write_region_to_bed_file(chrom, region_start_pos, (start_pos - 1),
							chrom_to_no_ref_regions_dict, bed_output_file)

	if is_no_call and prev_feature_is_no_call and (not is_adjacent):
		# create and report gap region
		region_start_pos = (prev_feature_end_pos + 1)
		if region_start_pos >= start_pos:
			sys.exit("Internal error: region to print start pos >= current feature start pos")
		convert_idx_and_write_region_to_bed_file(chrom, region_start_pos, (start_pos - 1),
							chrom_to_no_ref_regions_dict, bed_output_file)

	if is_call and prev_feature_is_no_call:
		## start new region
		region_start_pos = prev_feature_end_pos + 1

	return region_start_pos

def convert_idx_and_write_region_to_bed_file(chrom, original_start_pos_1bc, original_end_pos_1bc, chrom_to_no_ref_regions_dict, bed_output_file):
	# no ref regions are currently stored in vcf format
	region_to_print_list = trim_and_split_region_to_print_by_no_ref_regions(chrom, original_start_pos_1bc, original_end_pos_1bc,
										chrom_to_no_ref_regions_dict)

	for region_to_print in region_to_print_list:
		[start_pos_1bc, end_pos_1bc] = region_to_print

		# convert from 1-based closed (vcf format) to 0-based half-open (bed format)
		start_pos_0bho = start_pos_1bc - 1
		end_pos_0bho = end_pos_1bc
        
		# report region-to-print in bed file
		bed_output_file.write('%s\t%s\t%s\n' % (chrom, start_pos_0bho, end_pos_0bho))
	return

## both no-ref regions and start_pos/end_pos should be stored in vcf format at this point
def trim_and_split_region_to_print_by_no_ref_regions(chrom, original_start_pos, original_end_pos, chrom_to_no_ref_regions_dict):
	if not (chrom in chrom_to_no_ref_regions_dict):
		no_ref_regions_for_chrom = []
	else:
		no_ref_regions_for_chrom = chrom_to_no_ref_regions_dict[chrom]

	region_to_print_list = [[original_start_pos, original_end_pos]]

	## go over each no-ref region, one at a time
	for no_ref_region in no_ref_regions_for_chrom:
		[nrr_start_pos, nrr_end_pos] = no_ref_region

		## go over each region to print, one at a time
		temp_region_to_print_list = []
		for region in region_to_print_list:
			[start_pos, end_pos] = region

			if (end_pos < nrr_start_pos) or (start_pos > nrr_end_pos):
				# region to print does not overlap with no-ref region,
				# so keep region to print as-is
				temp_region_to_print_list.append([start_pos, end_pos])
			else:
				# there is an overlap, so consider the type of overlap:
				if end_pos > nrr_end_pos:
					# part of region to print is after no-ref region
					# so keep this part
					temp_region_to_print_list.append([nrr_end_pos + 1, end_pos])
				if start_pos < nrr_start_pos:
					# part of region to print is before no-ref region
					# so keep this part
					temp_region_to_print_list.append([start_pos, nrr_start_pos - 1])
		region_to_print_list = temp_region_to_print_list

	## sort regions to print
	region_to_print_list.sort()

	return region_to_print_list

def update_and_print_regions_of_interest_for_end_of_chromosome_for_last_feature(prev_feature_is_call, prev_feature_chrom, prev_feature_end_pos,
								region_start_pos, user_options_dict, chrom_to_length_dict,
								chrom_to_highest_observed_end_position_of_a_feature_dict,
								chrom_to_no_ref_regions_dict, bed_output_file):

	## get convenient term to use in upcoming conditional statement
	prev_feature_is_no_call = (not prev_feature_is_call)

	if not user_options_dict[PRINT_NO_CALLS]:
		region_start_pos = update_and_print_call_regions_for_end_of_old_chromosome(
										prev_feature_is_no_call, prev_feature_chrom, prev_feature_end_pos,
										region_start_pos, chrom_to_length_dict,
										chrom_to_no_ref_regions_dict, bed_output_file)

	return region_start_pos


