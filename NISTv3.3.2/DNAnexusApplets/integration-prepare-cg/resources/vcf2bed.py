#!/usr/bin/env python

import sys, re
from vcf2bed_helper_module import *

if not (len(sys.argv) in range(MIN_NUM_FIELDS_IN_ARGV, MAX_NUM_FIELDS_IN_ARGV + 1)):
	usage_string = "Usage: ./vcf2bed.py vcf_input_file no_ref_regions_input_file bed_output_file [ options ]"
	usage_string += '\n\nOPTIONS:'
	usage_string += '\n    -' + '\n    -'.join(OPTIONS_LIST)
	sys.exit(usage_string)

def main(vcf_input_filename, no_ref_regions_input_filename, bed_output_filename, user_options_list):

	global CURR_LINE

	# parse through user-specified options
	user_options_dict = get_user_options_dict(user_options_list)
	# read in no-ref regions
	chrom_to_no_ref_regions_dict = get_chrom_to_no_ref_regions_dict(no_ref_regions_input_filename)
	# initialize variables to use in main loop
	prev_feature, chrom_to_length_dict, chrom_to_highest_observed_end_position_of_a_feature_dict, \
			call_type_to_values_of_filter_to_use_dict, region_start_pos = \
					initialize_variables_for_main_loop()
	# open input and output files and write header line to bed output file
	vcf_input_file = open(vcf_input_filename, "r")
	bed_output_file = open(bed_output_filename, "w")
	bed_output_file.write('%s\t%s\t%s\n' % ("chrom", "chromStart", "chromEnd"))


	for line in vcf_input_file:
		line = line.rstrip('\n\r')
		CURR_LINE = line

		features_to_consider_list = preprocess_line_and_extract_features_to_consider(line, user_options_dict, chrom_to_length_dict, 
										call_type_to_values_of_filter_to_use_dict)

       
		for feature in features_to_consider_list:
			#print feature

			[chrom, start_pos, end_pos, is_call, call_type, is_padding_base, filter] = feature
			[prev_feature_chrom, prev_feature_start_pos, prev_feature_end_pos, prev_feature_is_call, \
					prev_feature_call_type, prev_feature_is_padding_base, prev_feature_filter] = prev_feature

			## assert that we do not have any features that are just a padding region
			if is_padding_base:
				sys.exit("We should not have any features that are just a padding region")

			# check assumptions about the length of the current feature
			# and its spatial relationship to the previous feature
			check_assumptions_about_current_and_previous_feature(prev_feature_chrom, prev_feature_end_pos, chrom, start_pos, end_pos)
                        
			# update and print regions of interest, based on current feature --
			# regions of interest will be either call regions or no-call regions,
			# depending on the value of the user-specified option "-print_no-calls"
			region_start_pos = \
				update_and_print_regions_of_interest(prev_feature_is_call, prev_feature_chrom, prev_feature_end_pos,
									is_call, chrom, start_pos, end_pos, 
									region_start_pos, user_options_dict, chrom_to_length_dict,
									chrom_to_highest_observed_end_position_of_a_feature_dict,
									chrom_to_no_ref_regions_dict, bed_output_file)

                                
                        
			## keep track of certain peices of information about most recent feature
			prev_feature = feature


	# assign variables providing information about previous vcf feature
	[prev_feature_chrom, prev_feature_start_pos, prev_feature_end_pos, prev_feature_is_call, \
			prev_feature_call_type, prev_feature_is_padding_base, prev_feature_filter] = prev_feature

	# update and print region of interest for end of chromosome for last feature
	region_start_pos = \
		update_and_print_regions_of_interest_for_end_of_chromosome_for_last_feature(
							prev_feature_is_call, prev_feature_chrom, prev_feature_end_pos,
							region_start_pos, user_options_dict, chrom_to_length_dict,
							chrom_to_highest_observed_end_position_of_a_feature_dict,
							chrom_to_no_ref_regions_dict, bed_output_file)


	## close input and output files
	vcf_input_file.close()
	bed_output_file.close()

	return

######################################################
## MAIN LOOP - PREPROCESS LINE AND EXTRACT ENTRIES: ##
######################################################

def preprocess_line_and_extract_features_to_consider(line, user_options_dict, chrom_to_length_dict, 
							call_type_to_values_of_filter_to_use_dict):

	features_to_consider_list = []

	## extract information about contig/chromosome length from certain header lines,
	## but do not search for call/no-call information on header lines
	if (len(line) != 0) and (line[0] == '#'):
		check_vcf_line_for_contig_information_and_update_chrom_length_dict(line, chrom_to_length_dict)
	else:
		## extract fields of interest
		field_list = get_field_list_for_line(line, 10)
		[chrom, POS, ID, ref_string, alt_string, QUAL, filter_string, info_string, format_string, value_string] = field_list
		preliminary_start_pos = int(POS)
        
		## get dictionary of key:value pairs from the "INFO" field
		info_key_to_value_dict = get_info_key_to_value_dict(info_string)
        
		## get dictionary of key:value pairs from the "FORMAT" and "GENOTYPE" fields (code assumes there is exactly one listed genotype)
		format_key_to_value_dict = get_format_key_to_value_dict(format_string, value_string)
        
		## skip line for certain types of features (i.e. mobile elts, copy number analysis windows, and structural variants)
		skip_line = determine_whether_to_skip_the_feature_on_this_line(alt_string, info_key_to_value_dict, line)
		if not skip_line:

			## get base padding offset
			base_padding_offset = get_base_padding_offset(ref_string, alt_string)

			if base_padding_offset < 0:
				sys.exit("Internal error: base padding offset should be >= 0")

			## Set the start position as the preliminary start position
			## plus the base padding offset, with the exception that padding regions
			## of size greater than one should be considered as part of the main entry.
			if base_padding_offset <= 1:
				start_pos = preliminary_start_pos + base_padding_offset
			else:
				start_pos = preliminary_start_pos
        
			## determine end position
			end_pos = get_end_position(preliminary_start_pos, ref_string, info_key_to_value_dict)
			
			## decide whether feature is a call or a no-call (based on both genotypic information and user-specified options)
			is_call, call_type, filter_to_use = \
					decide_whether_feature_is_a_call_or_no_call(filter_string, format_key_to_value_dict, user_options_dict,
										call_type_to_values_of_filter_to_use_dict)

			## do not consider the padding region (if any) to be a separate feature
			append_main_feature_to_features_to_consider_list(chrom, start_pos, end_pos, is_call, call_type,
									features_to_consider_list, filter_to_use)

	return features_to_consider_list

def get_field_list_for_line(line, expected_number_of_fields):
	## get each field
	field_list = line.split('\t')
	if len(field_list) != expected_number_of_fields:
		sys.exit("expected exactly %d fields (including exactly one 'genotype value' field) on every non-header line" % \
			 expected_number_of_fields)

	## strip away whitespace from each field
	for field_idx in range(0, len(field_list)):
		field_list[field_idx] = field_list[field_idx].strip()
	return field_list

def get_info_key_to_value_dict(info_string):
	info_entry_list = info_string.split(';')

	info_key_to_value_dict = {}
	for info_entry in info_entry_list:
		ie_key_value_list = info_entry.split("=")
		if len(ie_key_value_list) == 1:
			## if there is no "=" in this INFO entry
			## then make the entire entry a key,
			## and associate that key with a placeholder value of ""
			info_key = ie_key_value_list[0]
			info_value = ""
		elif len(ie_key_value_list) == 2:
			## if there is a single '=', then split the entry around the '='
			## to get a key and value pair
			info_key = ie_key_value_list[0]
			info_value = ie_key_value_list[1]
		else:
			info_key = -1
			info_value = -1
			sys.exit("expected each entry in INFO to contain at most one '=' character")
		info_key_to_value_dict[info_key] = info_value
	return info_key_to_value_dict

def get_format_key_to_value_dict(format_string, value_string):

	format_entry_list = format_string.split(':')
	value_entry_list = value_string.split(':')
	if len(format_entry_list) != len(value_entry_list):
		sys.exit("number of format entries does not equal number of value entries")

	format_key_to_value_dict = {}
	for format_entry_idx in range(0, len(format_entry_list)):
		format_key = format_entry_list[format_entry_idx]
		format_value = value_entry_list[format_entry_idx]
		format_key_to_value_dict[format_key] = format_value

	return format_key_to_value_dict

def determine_whether_to_skip_the_feature_on_this_line(alt_string, info_key_to_value_dict, line):
	skip_line = False

	## types of features that Rebecca said I should discard
	if re.match(r".*INS:ME.*", alt_string):
		# insertion of mobile element
		skip_line = True
	if re.match(r".*CGA_CNVWIN.*", alt_string):
		# copy number analysis window
		skip_line = True
	if SVTYPE in info_key_to_value_dict:
		# structural variant
		skip_line = True

	## additional types of features that looked suspicious to me
	for suspect_info_key in SUSPECT_INFO_KEYS:
		if suspect_info_key in info_key_to_value_dict:
			# info key that suggests it is possible that the feature should be discarded
			if not skip_line:
				print "Warning: Ignoring possible structural variant not detected by earlier filters:\n\"%s\"\n" % line
			skip_line = True

	return skip_line

main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])

