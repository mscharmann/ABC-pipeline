#!/usr/local/bin/python
# Python 3
# 2020
# 


# Mathias Scharmann
# added feature: for twophase models, Nx_ancient can be specified as LARGER or SMALLER


# a common issue:
# IOError: [Errno 32] Broken pipe
# is caused by mis-specification of $msnsam_arg2!

"""

v="-t tbs[theta] -r tbs[rec_rate] tbs[nsites] -I 3[npop] tbs[nsam1] tbs[nsam2] tbs[nsam3] 0[symm_migration_matrix]
-en 0.0[present time] 1[popID] tbs[N1_recent] 
-eg 0.0[present time] 1[popID] tbs[p1_growth_rate_during_recent_phase !unfree: constant OR exponential rate given by N1_recent and N1_ancient; controlled in argfile!] 
-en 0.0[present time] 2[popID] tbs[N2_recent] 
-eg 0.0[present time] 2[popID] tbs[p2_growth_rate_during_recent_phase !unfree: constant OR exponential rate given by N2_recent and N2_ancient; controlled in argfile!]
-en 0.0[present time] 3[popID] tbs[N3_recent]
-eg 0.0[present time] 3[popID] tbs[p3_growth_rate_during_recent_phase !unfree: constant OR exponential rate given by N3_recent and N3_ancient; controlled in argfile!] 
-en tbs[T_N1_ancient] 1[popID] tbs[N1_ancient]
-en tbs[T_N2_ancient] 2[popID] tbs[N2_ancient]
-en tbs[T_N3_ancient] 3[popID] tbs[N3_ancient]
-em tbs[T_m12_start] 1[recipient_ID] 2[donor_ID] tbs[m12 !given by hyper-parameters m12_scale and m12_prop_mig!]
-em tbs[T_m12_stop] 1[recipient_ID] 2[donor_ID] 0[mig_rate]
-em tbs[T_m13_start] 1[recipient_ID] 3[donor_ID] tbs[m13 !given by hyper-parameters m13_scale and m13_prop_mig!]
-em tbs[T_m13_stop] 1[recipient_ID] 3[donor_ID] 0[mig_rate]
-em tbs[T_m21_start] 2[recipient_ID] 1[donor_ID] tbs[m21 !given by hyper-parameters m21_scale and m21_prop_mig!]
-em tbs[T_m21_stop] 2[recipient_ID] 1[donor_ID] 0[mig_rate]
-em tbs[T_m23_start] 2[recipient_ID] 3[donor_ID] tbs[m23 !given by hyper-parameters m23_scale and m23_prop_mig!]
-em tbs[T_m23_stop] 2[recipient_ID] 3[donor_ID] 0[mig_rate]
-em tbs[T_m31_start] 3[recipient_ID] 1[donor_ID] tbs[m31 !given by hyper-parameters m31_scale and m31_prop_mig!]
-em tbs[T_m31_stop] 3[recipient_ID] 1[donor_ID] 0[mig_rate]
-em tbs[T_m32_start] 3[recipient_ID] 2[donor_ID] tbs[m32 !given by hyper-parameters m32_scale and m32_prop_mig!]
-em tbs[T_m32_stop] 3[recipient_ID] 2[donor_ID] 0[mig_rate]
-ej tbs[T_merge_2_into_1] 2[popID] 1[popID]
-en tbs[T_merge_2_into_1] 1[popID] tbs[pop_size_of_anc12, constant]
-eM tbs[T_merge_2_into_1] 0[symm_migration_matrix, i.e. all migration stops at tau_merge_2_into_1]
-ej tbs[T_merge_anc12_into_3] 1[popID] 3[popID]
-en tbs[T_merge_anc12_into_3] 3[popID] tbs[pop_size_of_anc123, constant]"


"""

# example:
# pseudocode:

#### HEAD


import sys
import os
import argparse
import numpy
import subprocess

# collects command line arguments and checks plausibility
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()

	parser.add_argument("-bpfile", help="bpfile: contains the fixed locus-specific parameters", type=extant_file, required=True, metavar="FILE")	
	parser.add_argument("-argfile", help="argfile: defines n_simulations and prior ranges", type=extant_file, required=True, metavar="FILE")
	parser.add_argument("-relative_N", help="distribution of relative Ne values for contigs in the genome; locus-specific Ne will be sampled from this distribution", type=extant_file, required=True, metavar="FILE")
	
	args = parser.parse_args()

	# finish
	return args

# checks if file exists and breaks script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print ("Error: {0} does not exist".format(x))
		exit()
	x = str(x)
	return x

### BODY

def draw_from_truncated_lognormal(mode_of_lnorm, sigma, lowerbound, upperbound):
	
	sigma = float(sigma)
	mu = numpy.log( float(mode_of_lnorm) ) + sigma*sigma
	lowerbound = float(lowerbound)
	upperbound = float(upperbound)
	
	while 1:
		x = numpy.random.lognormal( mu, float(sigma) )
		if x >= lowerbound and x <= upperbound:
			break
			
	return x

	
def draw_from_maxcutoff_flat_range( lowerbound, upperbound, maxcutoff ):
	x = maxcutoff + 1
	while x > maxcutoff:
		x = numpy.random.uniform( lowerbound, upperbound )
	return x


def draw_from_mincutoff_flat_range( lowerbound, upperbound, mincutoff ):
	x = mincutoff - 1
	while x < mincutoff:
		x = numpy.random.uniform( lowerbound, upperbound )
	return x


def draw_and_pipe_priors (n_iterations, loci_properties, prior_args, relative_N_distr):
		
	migration_shape = prior_args["migration_shape"]
	
	loci_properties_unmod = loci_properties
	
	# EDIT this:
	param_log_header = ["N1_recent", "N2_recent", "N3_recent", "N1_ancient", "N2_ancient", "N3_ancient", "Nanc12", "Nanc123","T_N1_ancient", "T_N2_ancient" , "T_N3_ancient",
	"T_merge_2_into_1",
	"T_merge_anc12_into_3",
	"T_m12_start","T_m12_stop",
"T_m13_start",
"T_m13_stop",
"T_m21_start",
"T_m21_stop",
"T_m23_start",
"T_m23_stop",
"T_m31_start",
"T_m31_stop",
"T_m32_start",
"T_m32_stop",
"m12_scale",
"m13_scale",
"m21_scale",
"m23_scale",
"m31_scale",
"m32_scale",
"m12_prop_mig",
"m13_prop_mig",
"m21_prop_mig",
"m23_prop_mig",
"m31_prop_mig",
"m32_prop_mig",
"m12beta1",
"m12beta2",
"m13beta1",
"m13beta2",
"m21beta1",
"m21beta2",
"m23beta1",
"m23beta2",
"m31beta1",
"m31beta2",
"m32beta1",
"m32beta2"]
	
	with open("parameters.txt", "w") as PARAMETERS_OUT:
		PARAMETERS_OUT.write( "\t".join( param_log_header ) + "\n" )
		
		for i in range(n_iterations) : # argfile specifies how many sets!
			
			# draw independent priors (some pop sizes and scales of migration rates)
			if prior_args["N_priors_shape"] == "flat":
				N1_recent = numpy.random.uniform( prior_args["N1_recent"][0], prior_args["N1_recent"][1] )
				N2_recent = numpy.random.uniform( prior_args["N2_recent"][0], prior_args["N2_recent"][1] )
				N3_recent = numpy.random.uniform( prior_args["N3_recent"][0], prior_args["N3_recent"][1] )
				Nanc12 = numpy.random.uniform( prior_args["Nanc12"][0], prior_args["Nanc12"][1] )
				Nanc123 = numpy.random.uniform( prior_args["Nanc123"][0], prior_args["Nanc123"][1] )
			elif prior_args["N_priors_shape"] == "lognormal":
				N1_recent = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_recent"][0], prior_args["N1_recent"][1] )
				N2_recent = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N2_recent"][0], prior_args["N2_recent"][1] )
				N3_recent = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N3_recent"][0], prior_args["N3_recent"][1] )
				Nanc12 = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["Nanc12"][0], prior_args["Nanc12"][1] )
				Nanc123 = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["Nanc123"][0], prior_args["Nanc123"][1] )
				
			
			# migration scales draw			
			if not prior_args["m12_scale"][0] == 0.0:
				m12_scale = 4*100000*(10**numpy.random.uniform( prior_args["m12_scale"][0], prior_args["m12_scale"][1] ))
			else: 
				m12_scale = 0.0

			if not prior_args["m21_scale"][0] == 0.0:
				m21_scale = 4*100000*(10**numpy.random.uniform( prior_args["m21_scale"][0], prior_args["m21_scale"][1] ))
			else: 
				m21_scale = 0.0			

			if not prior_args["m13_scale"][0] == 0.0:
				m13_scale = 4*100000*(10**numpy.random.uniform( prior_args["m13_scale"][0], prior_args["m13_scale"][1] ))
			else: 
				m13_scale = 0.0			

			if not prior_args["m31_scale"][0] == 0.0:
				m31_scale = 4*100000*(10**numpy.random.uniform( prior_args["m31_scale"][0], prior_args["m31_scale"][1] ))
			else: 
				m31_scale = 0.0			

			if not prior_args["m23_scale"][0] == 0.0:
				m23_scale = 4*100000*(10**numpy.random.uniform( prior_args["m23_scale"][0], prior_args["m23_scale"][1] ))
			else: 
				m23_scale = 0.0			

			if not prior_args["m32_scale"][0] == 0.0:
				m32_scale = 4*100000*(10**numpy.random.uniform( prior_args["m32_scale"][0], prior_args["m32_scale"][1] ))
			else: 
				m32_scale = 0.0			

				
			# draw interdependent priors (Time-priors)
			T_merge_2_into_1 = numpy.random.uniform( prior_args["T_merge_2_into_1"][0],  prior_args["T_merge_2_into_1"][1])
			T_merge_anc12_into_3 = draw_from_mincutoff_flat_range( prior_args["T_merge_anc12_into_3"][0],  prior_args["T_merge_anc12_into_3"][1], T_merge_2_into_1)
				# do not simply draw T_merge_anc12_into_3 at random with T_merge_2_into_1 as the lower bound, because this would over-ride an eventual independent lower bound specificed for it
									
			try: # the special flag is caught here:
				T_m12_stop = draw_from_maxcutoff_flat_range( prior_args["T_m12_stop"][0], prior_args["T_m12_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m12_stop = T_merge_2_into_1
			T_m12_start = draw_from_maxcutoff_flat_range( prior_args["T_m12_start"][0], prior_args["T_m12_start"][1], T_m12_stop) 
			try: 
				T_m13_stop = draw_from_maxcutoff_flat_range( prior_args["T_m13_stop"][0], prior_args["T_m13_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m13_stop = T_merge_2_into_1
			T_m13_start = draw_from_maxcutoff_flat_range( prior_args["T_m13_start"][0], prior_args["T_m13_start"][1], T_m13_stop) 
			try: 
				T_m21_stop = draw_from_maxcutoff_flat_range( prior_args["T_m21_stop"][0], prior_args["T_m21_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m21_stop = T_merge_2_into_1
			T_m21_start = draw_from_maxcutoff_flat_range( prior_args["T_m21_start"][0], prior_args["T_m21_start"][1], T_m21_stop) 
			try: 
				T_m23_stop = draw_from_maxcutoff_flat_range( prior_args["T_m23_stop"][0], prior_args["T_m23_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m23_stop = T_merge_2_into_1
			T_m23_start = draw_from_maxcutoff_flat_range( prior_args["T_m23_start"][0], prior_args["T_m23_start"][1], T_m23_stop) 
			try: 
				T_m31_stop = draw_from_maxcutoff_flat_range( prior_args["T_m31_stop"][0], prior_args["T_m31_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m31_stop = T_merge_2_into_1
			T_m31_start = draw_from_maxcutoff_flat_range( prior_args["T_m31_start"][0], prior_args["T_m31_start"][1], T_m31_stop) 
			try: 
				T_m32_stop = draw_from_maxcutoff_flat_range( prior_args["T_m32_stop"][0], prior_args["T_m32_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m32_stop = T_merge_2_into_1
			T_m32_start = draw_from_maxcutoff_flat_range( prior_args["T_m32_start"][0], prior_args["T_m32_start"][1], T_m32_stop) 


			T_N1_ancient = draw_from_maxcutoff_flat_range( prior_args["T_N1_ancient"][0], prior_args["T_N1_ancient"][1], T_merge_2_into_1)
			T_N2_ancient = draw_from_maxcutoff_flat_range( prior_args["T_N2_ancient"][0], prior_args["T_N2_ancient"][1], T_merge_2_into_1) 
			T_N3_ancient = draw_from_maxcutoff_flat_range( prior_args["T_N3_ancient"][0], prior_args["T_N3_ancient"][1], T_merge_anc12_into_3) 

			
			## get the N_ancient sizes and catch special flag which will also affect growth rate!! (has to be defined as zero); constant size means growth rate zero!
			# the exponential growth rate parameter, which depends on N1_recent, N1_ancient and T_N1_ancient. This is NOT A FREE PARAMETER but necessary to supply it to ms to run!
			# end_size = start_size * exp^-growth_rate*end_time ### end_size and end_time are both BACKWARD in time!! 
			# thus growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent) ## natural logarithm
			try: # the special flag is caught here for N1_ancient:
				if prior_args["N_priors_shape"] == "flat":
					N1_ancient = numpy.random.uniform( prior_args["N1_ancient"][0], prior_args["N1_ancient"][1] )
				elif prior_args["N_priors_shape"] == "lognormal":
					N1_ancient = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_ancient"][0], prior_args["N1_ancient"][1] )
				growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent)
				growth_rate_1 = str( growth_rate_1 )
			except ValueError: # this error means that N1_ancient has no prior range but is set relative to N1 (identical, smaller or larger)
				if prior_args["N1_ancient"] == "N1_recent":# to be identical to N1_recent => makes a single-phase model!
					N1_ancient = N1_recent
					growth_rate_1 = "0.0"
				elif prior_args["N1_ancient"] == "larger_than_recent":
					N1_ancient = numpy.random.uniform( N1_recent, prior_args["N1_recent"][1] )
					growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent)
					growth_rate_1 = str( growth_rate_1 )
				elif prior_args["N1_ancient"] == "smaller_than_recent":
					N1_ancient = numpy.random.uniform( prior_args["N1_recent"][0], N1_recent )
					growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent)
					growth_rate_1 = str( growth_rate_1 )

			try: # the special flag is caught here for N2_ancient:
				if prior_args["N_priors_shape"] == "flat":
					N2_ancient = numpy.random.uniform( prior_args["N2_ancient"][0], prior_args["N2_ancient"][1] )
				elif prior_args["N_priors_shape"] == "lognormal":
					N2_ancient = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N2_ancient"][0], prior_args["N2_ancient"][1] )
				growth_rate_2 = -(1.0/T_N2_ancient)*numpy.log(N2_ancient/N2_recent)
				growth_rate_2 = str( growth_rate_2 )	
			except ValueError:  # this error means that N2_ancient has no prior range but is set relative to N2 (identical, smaller or larger)
				if prior_args["N2_ancient"] == "N2_recent":# to be identical to N2_recent => makes a single-phase model!
					N2_ancient = N2_recent
					growth_rate_2 = "0.0"
				elif prior_args["N2_ancient"] == "larger_than_recent":
					N2_ancient = numpy.random.uniform( N2_recent, prior_args["N2_recent"][1] )
					growth_rate_2 = -(1.0/T_N2_ancient)*numpy.log(N2_ancient/N2_recent)
					growth_rate_2 = str( growth_rate_2 )	
				elif prior_args["N2_ancient"] == "smaller_than_recent":
					N2_ancient = numpy.random.uniform( prior_args["N2_recent"][0], N2_recent )
					growth_rate_2 = -(1.0/T_N2_ancient)*numpy.log(N2_ancient/N2_recent)
					growth_rate_2 = str( growth_rate_2 )	

			try: # the special flag is caught here for N3_ancient:
				if prior_args["N_priors_shape"] == "flat":
					N3_ancient = numpy.random.uniform( prior_args["N3_ancient"][0], prior_args["N3_ancient"][1] )
				elif prior_args["N_priors_shape"] == "lognormal":
					N3_ancient = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N3_ancient"][0], prior_args["N3_ancient"][1] )
				growth_rate_3 = -(1.0/T_N3_ancient)*numpy.log(N3_ancient/N3_recent)
				growth_rate_3 = str( growth_rate_3 )	
			except ValueError:  # this error means that N3_ancient has no prior range but is set relative to N3 (identical, smaller or larger)
				if prior_args["N3_ancient"] == "N3_recent":# to be identical to N3_recent => makes a single-phase model!
					N3_ancient = N3_recent
					growth_rate_3 = "0.0"
				elif prior_args["N3_ancient"] == "larger_than_recent":
					N3_ancient = numpy.random.uniform( N3_recent, prior_args["N3_recent"][1] )
					growth_rate_3 = -(1.0/T_N3_ancient)*numpy.log(N3_ancient/N3_recent)
					growth_rate_3 = str( growth_rate_3 )	
				elif prior_args["N3_ancient"] == "smaller_than_recent":
					N3_ancient = numpy.random.uniform( prior_args["N3_recent"][0], N3_recent )
					growth_rate_3 = -(1.0/T_N3_ancient)*numpy.log(N3_ancient/N3_recent)
					growth_rate_3 = str( growth_rate_3 )	


			# now the locus-specific stuff:
			# get the migration rates first, this may be clever for speed if all the draws are done before the stdout loop, but have to test that
			nloci = len( loci_properties )			
			
			if migration_shape == "homogenous":
				m12_list = numpy.array( [1]*nloci ) * m12_scale
				m13_list = numpy.array( [1]*nloci ) * m13_scale
				m21_list = numpy.array( [1]*nloci ) * m21_scale
				m23_list = numpy.array( [1]*nloci ) * m23_scale
				m31_list = numpy.array( [1]*nloci ) * m31_scale
				m32_list = numpy.array( [1]*nloci ) * m32_scale
				
				m12beta1 = m12beta2 = m13beta1 = m13beta2 = m21beta1 = m21beta2 = m23beta1 = m23beta2 =m31beta1 = m31beta2 = m32beta1 = m32beta2 = "NA"
				m12_prop_mig = m13_prop_mig = m21_prop_mig = m23_prop_mig = m31_prop_mig = m32_prop_mig = "NA"
				
				
			elif migration_shape == "heterogenous_beta":
				
				# hyperpriors:
				m12beta1 = numpy.random.uniform( prior_args["m12beta1"][0], prior_args["m12beta1"][1] )
				m12beta2 = numpy.random.uniform( prior_args["m12beta2"][0], prior_args["m12beta2"][1] )
				m13beta1 = numpy.random.uniform( prior_args["m13beta1"][0], prior_args["m13beta1"][1] )
				m13beta2 = numpy.random.uniform( prior_args["m13beta2"][0], prior_args["m13beta2"][1] )
				m21beta1 = numpy.random.uniform( prior_args["m21beta1"][0], prior_args["m21beta1"][1] )
				m21beta2 = numpy.random.uniform( prior_args["m21beta2"][0], prior_args["m21beta2"][1] )
				m23beta1 = numpy.random.uniform( prior_args["m23beta1"][0], prior_args["m23beta1"][1] )
				m23beta2 = numpy.random.uniform( prior_args["m23beta2"][0], prior_args["m23beta2"][1] )
				m31beta1 = numpy.random.uniform( prior_args["m31beta1"][0], prior_args["m31beta1"][1] )
				m31beta2 = numpy.random.uniform( prior_args["m31beta2"][0], prior_args["m31beta2"][1] )
				m32beta1 = numpy.random.uniform( prior_args["m32beta1"][0], prior_args["m32beta1"][1] )
				m32beta2 = numpy.random.uniform( prior_args["m32beta2"][0], prior_args["m32beta2"][1] )
			
				m12_list = numpy.random.beta( m12beta1 , m12beta2, nloci) * m12_scale		
				m13_list = numpy.random.beta( m13beta1 , m13beta2, nloci) * m13_scale		
				m21_list = numpy.random.beta( m21beta1 , m21beta2, nloci) * m21_scale
				m23_list = numpy.random.beta( m23beta1 , m23beta2, nloci) * m23_scale
				m31_list = numpy.random.beta( m31beta1 , m31beta2, nloci) * m31_scale
				m32_list = numpy.random.beta( m32beta1 , m32beta2, nloci) * m32_scale
				
				m12_prop_mig = m13_prop_mig = m21_prop_mig = m23_prop_mig = m31_prop_mig = m32_prop_mig = "NA"
				
			elif migration_shape == "heterogenous_binomial":
				
				m12_prop_mig = numpy.random.uniform( prior_args["m12_prop_mig"][0], prior_args["m12_prop_mig"][1] )
				m13_prop_mig = numpy.random.uniform( prior_args["m13_prop_mig"][0], prior_args["m13_prop_mig"][1] )
				m21_prop_mig = numpy.random.uniform( prior_args["m21_prop_mig"][0], prior_args["m21_prop_mig"][1] )
				m23_prop_mig = numpy.random.uniform( prior_args["m23_prop_mig"][0], prior_args["m23_prop_mig"][1] )
				m31_prop_mig = numpy.random.uniform( prior_args["m31_prop_mig"][0], prior_args["m31_prop_mig"][1] )
				m32_prop_mig = numpy.random.uniform( prior_args["m32_prop_mig"][0], prior_args["m32_prop_mig"][1] )
				
				m12_list = numpy.random.binomial( 1 , m12_prop_mig, nloci) * m12_scale		
				m13_list = numpy.random.binomial( 1 , m13_prop_mig, nloci) * m13_scale		
				m21_list = numpy.random.binomial( 1 , m21_prop_mig, nloci) * m21_scale
				m23_list = numpy.random.binomial( 1 , m23_prop_mig, nloci) * m23_scale
				m31_list = numpy.random.binomial( 1 , m31_prop_mig, nloci) * m31_scale
				m32_list = numpy.random.binomial( 1 , m32_prop_mig, nloci) * m32_scale
				
				m12beta1 = m12beta2 = m13beta1 = m13beta2 = m21beta1 = m21beta2 = m23beta1 = m23beta2 =m31beta1 = m31beta2 = m32beta1 = m32beta2 = "NA"

			## now N_hetero
			relative_Ns = numpy.random.choice(relative_N_distr, size=nloci, replace=True, p=None)
#			print (numpy.mean(relative_Ns))
			hetero_Ns = []
			for idx in range( nloci ):
				loc_Ns = [x*relative_Ns[idx] for x in [N1_recent, N2_recent, N3_recent, N1_ancient, N2_ancient, N3_ancient, Nanc12, Nanc123]]
#				print(idx,relative_Ns[idx],[N1_recent,N1_ancient,N2_recent,N2_ancient,Nanc12],loc_Ns)
				hetero_Ns.append( [str(round_except_TypeError(x, 5)) for x in loc_Ns ] )
#			print (hetero_Ns)

			parameters = [str(round_except_TypeError(x, 5)) for x in [		
			N1_recent, N2_recent, N3_recent, N1_ancient, N2_ancient, N3_ancient, Nanc12, Nanc123, T_N1_ancient, T_N2_ancient , T_N3_ancient, T_merge_2_into_1,	T_merge_anc12_into_3, T_m12_start,T_m12_stop,T_m13_start,T_m13_stop,T_m21_start,T_m21_stop,T_m23_start,T_m23_stop,T_m31_start,T_m31_stop,T_m32_start, T_m32_stop,
			
			 log10_except_zero(m12_scale/(4*100000)), log10_except_zero(m13_scale/(4*100000)), log10_except_zero(m21_scale/(4*100000)), log10_except_zero(m23_scale/(4*100000)), log10_except_zero(m31_scale/(4*100000)), log10_except_zero(m32_scale/(4*100000)),
			
			 m12_prop_mig, m13_prop_mig, m21_prop_mig, m23_prop_mig, m31_prop_mig, m32_prop_mig, m12beta1, m12beta2, m13beta1, m13beta2, m21beta1, m21beta2, m23beta1, m23beta2, m31beta1, m31beta2, m32beta1, m32beta2	] ]
			
			[N1_recent, N2_recent, N3_recent, N1_ancient, N2_ancient, N3_ancient, Nanc12, Nanc123, T_N1_ancient, T_N2_ancient , T_N3_ancient, T_merge_2_into_1,	T_merge_anc12_into_3, T_m12_start,T_m12_stop,T_m13_start,T_m13_stop,T_m21_start,T_m21_stop,T_m23_start,T_m23_stop,T_m31_start,T_m31_stop,T_m32_start, T_m32_stop, m12_scale, m13_scale, m21_scale, m23_scale, m31_scale, m32_scale, m12_prop_mig, m13_prop_mig, m21_prop_mig, m23_prop_mig, m31_prop_mig, m32_prop_mig, m12beta1, m12beta2, m13beta1, m13beta2, m21beta1, m21beta2, m23beta1, m23beta2, m31beta1, m31beta2, m32beta1, m32beta2] = parameters 			

			# now the locus-specific stuff:
			migrates = []			
			for idx in range( nloci ):

				migrates.append( [str(round_except_TypeError(x, 5)) for x in [ m12_list[idx], m13_list[idx] , m21_list[idx] , m23_list[idx] , m31_list[idx] , m32_list[idx] ]] )
		
			for idx, locus in enumerate(loci_properties):

				[locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, nsam_3, total_nsam_per_locus], [m12, m13, m21, m23, m31, m32], [N1_recent, N2_recent, N3_recent, N1_ancient, N2_ancient, N3_ancient, Nanc12, Nanc123] = loci_properties[idx], migrates[idx], hetero_Ns[idx]	

				# writes directly to stdout -> output can be piped to msnsam 
				sys.stdout.write(
				" ".join( [total_nsam_per_locus, locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2,  nsam_3,
				N1_recent, growth_rate_1, 
				N2_recent, growth_rate_2, 
				N3_recent, growth_rate_3, 
				T_N1_ancient, N1_ancient, 
				T_N2_ancient, N2_ancient, 
				T_N3_ancient, N3_ancient,
				T_m12_start, m12, 
				T_m12_stop,
				T_m13_start, m13 ,
				T_m13_stop,
				T_m21_start, m21 ,
				T_m21_stop,
				T_m23_start, m23 ,
				T_m23_stop,
				T_m31_start, m31 ,
				T_m31_stop,
				T_m32_start, m32 ,
				T_m32_stop,
				T_merge_2_into_1,
				T_merge_2_into_1, Nanc12,
				T_merge_2_into_1,
				T_merge_anc12_into_3,
				T_merge_anc12_into_3, Nanc123
				]) + "\n")
	
			
			# write prior realisations to parameters.txt log file
			PARAMETERS_OUT.write( "\t".join( parameters  ) + "\n" )

				


###
def log10_except_zero (infloat):
	
	if not infloat == 0.0:
		return(numpy.log10(infloat))
	else:
		return 0.0


def round_except_TypeError(invalue, n_digits):
	
	try:
		outvalue = round(invalue, n_digits)
	except TypeError:
		outvalue = invalue
	return outvalue
	

def read_bpfile(bpfile):
	
	## for 3 populations!
	
	with open(bpfile, "r") as INFILE:
		header = INFILE.readline() # remove first line
		locilengthline = INFILE.readline().strip("\n").split("\t")
		nsam1_line = INFILE.readline().strip("\n").split("\t")
		nsam2_line = INFILE.readline().strip("\n").split("\t")
		nsam3_line = INFILE.readline().strip("\n").split("\t")
		thetaline = INFILE.readline().strip("\n").split("\t")
		rholine = INFILE.readline().strip("\n").split("\t")
	
	# [locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, total_nsam_per_locus]

	loci_properties = []
	for i in range(len(locilengthline)):
		locus_theta = thetaline[i]		
		locus_rho = rholine[i]
		locus_nsites = locilengthline[i]
		nsam_1 = nsam1_line[i]
		nsam_2 = nsam2_line[i]
		nsam_3 = nsam3_line[i]
		loci_properties.append( [locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, nsam_3, str(int(nsam_1)+int(nsam_2)+int(nsam_3)) ] )
	
	return loci_properties	

def read_argfile(argfile):
	
	# EDIT this function
	
	nreps = None
	
	prior_args = {}		
	with open(argfile, "r") as INFILE:
		for line in INFILE:
			if not line.startswith("#"): # allow uncommented lines 			
				try:			
					parameter, entry = line.strip("\n").split(" = ")
				except ValueError:
					continue # allows for empty lines in argfile
					
				if parameter == "nreps":
					nreps = int(entry)
				elif parameter == "migration_shape":
					if entry in ["homogenous", "heterogenous_beta", "heterogenous_binomial"]:
						prior_args[parameter] = entry
					else:
						print ("migration_shape must be one of 'homogenous', 'heterogenous_beta' or 'heterogenous_binomial'")
						exit()
				elif parameter == "N_priors_shape":
					if entry in ["flat", "lognormal"]:
						prior_args[parameter] = entry
					else:
						print ("N_priors_shape must be one of 'flat' or 'lognormal'")
						exit()

				else:	
					try:		
						prior_args[parameter] = [ float(x) for x in entry.split(",") ]
					except ValueError:
						if parameter == "T_m12_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						if parameter == "T_m13_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						if parameter == "T_m21_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						if parameter == "T_m23_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						if parameter == "T_m31_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						if parameter == "T_m32_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						elif parameter == "N1_ancient" and entry == "N1_recent":
							prior_args[parameter] = entry
						elif parameter == "N1_ancient" and entry == "larger_than_recent":
							prior_args[parameter] = entry
						elif parameter == "N1_ancient" and entry == "smaller_than_recent":
							prior_args[parameter] = entry
						elif parameter == "N2_ancient" and entry == "N2_recent":
							prior_args[parameter] = entry
						elif parameter == "N2_ancient" and entry == "larger_than_recent":
							prior_args[parameter] = entry
						elif parameter == "N2_ancient" and entry == "smaller_than_recent":
							prior_args[parameter] = entry
						elif parameter == "N3_ancient" and entry == "N3_recent":
							prior_args[parameter] = entry
						elif parameter == "N3_ancient" and entry == "larger_than_recent":
							prior_args[parameter] = entry
						elif parameter == "N3_ancient" and entry == "smaller_than_recent":
							prior_args[parameter] = entry
						else:							
							print ("argfile format is offended by ", parameter)
							exit()

	# check that migration scale is meaningful
	assert 0 <= 10**prior_args["m12_scale"][0] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m12_scale"][1] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m21_scale"][0] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m21_scale"][1] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m13_scale"][0] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m13_scale"][1] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m31_scale"][0] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m31_scale"][1] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m23_scale"][0] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m23_scale"][1] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m32_scale"][0] <= 1.0, "fraction of immigrants not within [0,1]"
	assert 0 <= 10**prior_args["m32_scale"][1] <= 1.0, "fraction of immigrants not within [0,1]"
	
	if nreps == None:
		print ("argfile is missing specification of nreps")
		exit()
	
	# EDIT this:		
	params_to_be_read = [ "N1_recent", "N2_recent", "N3_recent", "N1_ancient", "N2_ancient", "N3_ancient", "Nanc12", "Nanc123","T_N1_ancient", "T_N2_ancient" , "T_N3_ancient",
	"T_merge_2_into_1",
	"T_merge_anc12_into_3",
	"T_m12_start","T_m12_stop",
"T_m13_start",
"T_m13_stop",
"T_m21_start",
"T_m21_stop",
"T_m23_start",
"T_m23_stop",
"T_m31_start",
"T_m31_stop",
"T_m32_start",
"T_m32_stop",
"m12_scale",
"m13_scale",
"m21_scale",
"m23_scale",
"m31_scale",
"m32_scale",
"migration_shape",
"m12_prop_mig",
"m13_prop_mig",
"m21_prop_mig",
"m23_prop_mig",
"m31_prop_mig",
"m32_prop_mig",
"m12beta1",
"m12beta2",
"m13beta1",
"m13beta2",
"m21beta1",
"m21beta2",
"m23beta1",
"m23beta2",
"m31beta1",
"m31beta2",
"m32beta1",
"m32beta2",
"N_priors_shape","N_mode_of_lognormal","N_sigma"]

	for param in params_to_be_read:
		if param not in prior_args.keys():
			print ("argfile is missing specification of ", param)
			exit()
			
	return nreps, prior_args

def modify_spinput_txt(n_iterations):
	
	# appends correct number of iterations and filename to spinput.txt for mscalc
	
	with open("spinput.txt", "r") as FILE:
		num_lines = sum(1 for line in FILE)
	
	with open("spinput.txt", "r") as FILE:
				
		outspinput = ""
		linecnt = 0
		while linecnt < (num_lines - 2) :
			inline = FILE.readline()
			outspinput += inline
			linecnt += 1
		
		outspinput += str(n_iterations) + "\n"
		outspinput += "myfifo\n"
	
	with open("spinput.txt", "w") as FILE:
		FILE.write(outspinput)
	
def read_relative_N (infile):

#	print("readling relative N distribution from file")
	
	hetero_N = []
	with open(infile,"r") as I:
		for line in I:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				hetero_N.append( float(fields[0]) )
	return(numpy.array(hetero_N))

	
####################			
### MAIN

# take the verbose of msnsam_model_cmd and remove all the ms flags etc until only the parameters to be generated by this python script are left:
all_free_parameters_in_order="""[N1_recent] 
[growth_rate_1 !unfree: constant OR exponential rate given by N1_recent and N1_ancient; controlled in argfile!] 
[N2_recent] 
[growth_rate_2 !unfree: constant OR exponential rate given by N2_recent and N2_ancient; controlled in argfile!]
[N3_recent]
[growth_rate_3 !unfree: constant OR exponential rate given by N3_recent and N3_ancient; controlled in argfile!] 
[T_N1_ancient] [N1_ancient]
[T_N2_ancient] [N2_ancient]
[T_N3_ancient] [N3_ancient]
[T_m12_start] [m12 !given by hyper-parameters m12_scale and m12_prop_mig!]
[T_m12_stop]
[T_m13_start] [m13 !given by hyper-parameters m13_scale and m13_prop_mig!]
[T_m13_stop]
[T_m21_start] [m21 !given by hyper-parameters m21_scale and m21_prop_mig!]
[T_m21_stop]
[T_m23_start] [m23 !given by hyper-parameters m23_scale and m23_prop_mig!]
[T_m23_stop]
[T_m31_start] [m31 !given by hyper-parameters m31_scale and m31_prop_mig!]
[T_m31_stop]
[T_m32_start] [m32 !given by hyper-parameters m32_scale and m32_prop_mig!]
[T_m32_stop]
[T_merge_2_into_1]
[T_merge_2_into_1] [Nanc12 !constant!]
[T_merge_2_into_1]
[T_merge_anc12_into_3]
[T_merge_anc12_into_3] [Nanc123 !constant!]
"""

# edit all_free_parameters_in_order until you have a clean python list of the parameters:
# TextWrangler: \] to , and \[ to nothing and \!.*?\! to nothing
STDOUT_structure = """N1_recent, growth_rate_1, N2_recent, growth_rate_2, N3_recent, growth_rate_3, 
T_N1_ancient, N1_ancient, T_N2_ancient, N2_ancient, T_N3_ancient, N3_ancient,
T_m12_start, m12, 
T_m12_stop,
T_m13_start, m13 ,
T_m13_stop,
T_m21_start, m21 ,
T_m21_stop,
T_m23_start, m23 ,
T_m23_stop,
T_m31_start, m31 ,
T_m31_stop,
T_m32_start, m32 ,
T_m32_stop,
T_merge_2_into_1,
T_merge_2_into_1, Nanc12,
T_merge_2_into_1,
T_merge_anc12_into_3,
T_merge_anc12_into_3, Nanc123"""

# now use this python command to get a python list: [x.strip() for x in STDOUT_structure.split(",")]
STDOUT_structure = ['N1_recent', 'growth_rate_1', 'N2_recent', 'growth_rate_2', 'N3_recent', 'growth_rate_3', 'T_N1_ancient', 'N1_ancient', 'T_N2_ancient', 'N2_ancient', 'T_N3_ancient', 'N3_ancient', 'T_m12_start', 'm12', 'T_m12_stop', 'T_m13_start', 'm13', 'T_m13_stop', 'T_m21_start', 'm21', 'T_m21_stop', 'T_m23_start', 'm23', 'T_m23_stop', 'T_m31_start', 'm31', 'T_m31_stop', 'T_m32_start', 'm32', 'T_m32_stop', 'T_merge_2_into_1', 'T_merge_2_into_1', 'Nanc12', 'T_merge_2_into_1', 'T_merge_anc12_into_3', 'T_merge_anc12_into_3', 'Nanc123']

parameters_log_structure = sorted( list(set(STDOUT_structure)) )


args = 	get_commandline_arguments ()

loci_properties = read_bpfile(args.bpfile)
#print loci_properties

relative_N_distr = read_relative_N(args.relative_N)

n_iterations, prior_args = read_argfile(args.argfile)
#print n_iterations, prior_args

draw_and_pipe_priors (n_iterations, loci_properties, prior_args, relative_N_distr)







