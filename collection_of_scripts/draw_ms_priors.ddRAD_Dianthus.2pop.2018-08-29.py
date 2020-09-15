#!/usr/local/bin/python
# Python 2.7.6
# 2018
# 


# Mathias Scharmann
# added feature: for twophase models, Nx_ancient can be specified as LARGER or SMALLER
# added feature: linked selection! => heterongenous Ne along the genome. Modelled as in Roux et al 2016 PLOS Biology.

# a common issue:
# IOError: [Errno 32] Broken pipe
# is caused by mis-specification of $msnsam_arg2!

"""

v="-t tbs[theta] -r tbs[rec_rate] tbs[nsites] -I 2[npop] tbs[nsam1] tbs[nsam2] 0[symm_migration_matrix]
-en 0.0[present time] 1[popID] tbs[N1_recent] 
-eg 0.0[present time] 1[popID] tbs[p1_growth_rate_during_recent_phase !unfree: constant OR exponential rate given by N1_recent and N1_ancient; controlled in argfile!] 
-en 0.0[present time] 2[popID] tbs[N2_recent] 
-eg 0.0[present time] 2[popID] tbs[p2_growth_rate_during_recent_phase !unfree: constant OR exponential rate given by N2_recent and N2_ancient; controlled in argfile!]
-en tbs[T_N1_ancient] 1[popID] tbs[N1_ancient]
-en tbs[T_N2_ancient] 2[popID] tbs[N2_ancient]
-em tbs[T_m12_start] 1[recipient_ID] 2[donor_ID] tbs[m12 !given by hyper-parameters m12_scale and m12_prop_mig!]
-em tbs[T_m12_stop] 1[recipient_ID] 2[donor_ID] 0[mig_rate]
-em tbs[T_m21_start] 2[recipient_ID] 1[donor_ID] tbs[m21 !given by hyper-parameters m21_scale and m21_prop_mig!]
-em tbs[T_m21_stop] 2[recipient_ID] 1[donor_ID] 0[mig_rate]
-ej tbs[T_merge_2_into_1] 2[popID] 1[popID]
-en tbs[T_merge_2_into_1] 1[popID] tbs[pop_size_of_anc12, constant]"


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
	
	args = parser.parse_args()

	# finish
	return args

# checks if file exists and breaks script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
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
		x = numpy.random.lognormal( mu, float(sigma), 1 )
		if x >= lowerbound and x <= upperbound:
			break
			
	return x

	
def draw_from_maxcutoff_flat_range( lowerbound, upperbound, maxcutoff ):
	x = maxcutoff + 1
	while x > maxcutoff:
		x = numpy.random.uniform( lowerbound, upperbound, 1 )
	return x


def draw_from_mincutoff_flat_range( lowerbound, upperbound, mincutoff ):
	x = mincutoff - 1
	while x < mincutoff:
		x = numpy.random.uniform( lowerbound, upperbound, 1 )
	return x


def draw_and_pipe_priors (n_iterations, loci_properties, prior_args):
		
	migration_shape = prior_args["migration_shape"]
	
	loci_properties_unmod = loci_properties
	
	# EDIT this:
	param_log_header = ["N1_recent", "N2_recent", "N1_ancient", "N2_ancient", "Nanc12", "T_N1_ancient", "T_N2_ancient" , 
	"T_merge_2_into_1",
	"T_m12_start","T_m12_stop",
"T_m21_start",
"T_m21_stop",
"m12_scale",
"m21_scale",
"m12_prop_mig",
"m21_prop_mig",
"m12beta1",
"m12beta2",
"m21beta1",
"m21beta2",
"N1_prop_hetero",
"N2_prop_hetero",
"Nanc12_prop_hetero",
"N_hetero_beta1",
"N_hetero_beta2"
]
	
	with open("parameters.txt", "w") as PARAMETERS_OUT:
		PARAMETERS_OUT.write( "\t".join( param_log_header ) + "\n" )
		
		for i in range(n_iterations) : # argfile specifies how many sets!
			
			# draw independent priors (some pop sizes, scales of migration rates, N_prop_hetero, shapes of N_hetero distribution)
			if prior_args["N_priors_shape"] == "flat":
				N1_recent = numpy.random.uniform( prior_args["N1_recent"][0], prior_args["N1_recent"][1], 1 )
				N2_recent = numpy.random.uniform( prior_args["N2_recent"][0], prior_args["N2_recent"][1], 1 )
				Nanc12 = numpy.random.uniform( prior_args["Nanc12"][0], prior_args["Nanc12"][1], 1 )
			elif prior_args["N_priors_shape"] == "lognormal":
				N1_recent = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_recent"][0], prior_args["N1_recent"][1] )
				N2_recent = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N2_recent"][0], prior_args["N2_recent"][1] )

				Nanc12 = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["Nanc12"][0], prior_args["Nanc12"][1] )
				
			
			m12_scale = numpy.random.uniform( prior_args["m12_scale"][0], prior_args["m12_scale"][1], 1 )
			m21_scale = numpy.random.uniform( prior_args["m21_scale"][0], prior_args["m21_scale"][1], 1 )
			
			N1_prop_hetero = numpy.random.uniform( prior_args["N_prop_hetero"][0], prior_args["N_prop_hetero"][1] )
			N2_prop_hetero = numpy.random.uniform( prior_args["N_prop_hetero"][0], prior_args["N_prop_hetero"][1] )
			Nanc12_prop_hetero = numpy.random.uniform( prior_args["N_prop_hetero"][0], prior_args["N_prop_hetero"][1] )
			N_hetero_beta1 = numpy.random.uniform( prior_args["N_hetero_beta1"][0], prior_args["N_hetero_beta1"][1] )
			N_hetero_beta2 = numpy.random.uniform( prior_args["N_hetero_beta2"][0], prior_args["N_hetero_beta2"][1] )

			# draw interdependent priors (Time-priors)
			T_merge_2_into_1 = numpy.random.uniform( prior_args["T_merge_2_into_1"][0],  prior_args["T_merge_2_into_1"][1], 1)
									
			try: # the special flag is caught here:
				T_m12_stop = draw_from_maxcutoff_flat_range( prior_args["T_m12_stop"][0], prior_args["T_m12_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m12_stop = T_merge_2_into_1
			T_m12_start = draw_from_maxcutoff_flat_range( prior_args["T_m12_start"][0], prior_args["T_m12_start"][1], T_m12_stop) 
			try: 
				T_m21_stop = draw_from_maxcutoff_flat_range( prior_args["T_m21_stop"][0], prior_args["T_m21_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m21_stop = T_merge_2_into_1
			T_m21_start = draw_from_maxcutoff_flat_range( prior_args["T_m21_start"][0], prior_args["T_m21_start"][1], T_m21_stop) 

			T_N1_ancient = draw_from_maxcutoff_flat_range( prior_args["T_N1_ancient"][0], prior_args["T_N1_ancient"][1], T_merge_2_into_1)
			T_N2_ancient = draw_from_maxcutoff_flat_range( prior_args["T_N2_ancient"][0], prior_args["T_N2_ancient"][1], T_merge_2_into_1) 
			
			## get the N_ancient sizes and catch special flag which will also affect growth rate!! (has to be defined as zero); constant size means growth rate zero!
			# the exponential growth rate parameter, which depends on N1_recent, N1_ancient and T_N1_ancient. This is NOT A FREE PARAMETER but necessary to supply it to ms to run!
			# end_size = start_size * exp^-growth_rate*end_time ### end_size and end_time are both BACKWARD in time!! 
			# thus growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent) ## natural logarithm
			# 						=> since N_hetero is implemented such that both N_ancient and N_recent are re-scaled by the same factor 
			#							(i.e 'linked selection' is constant in the two Ne phases),  
			#							the ratio N_ancient/N_recent remains the same,  and thus the growth_rate also does not need to be locus-specific.
			try: # the special flag is caught here for N1_ancient:
				if prior_args["N_priors_shape"] == "flat":
					N1_ancient = numpy.random.uniform( prior_args["N1_ancient"][0], prior_args["N1_ancient"][1], 1 )
				elif prior_args["N_priors_shape"] == "lognormal":
					N1_ancient = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_ancient"][0], prior_args["N1_ancient"][1] )
				growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent)
				growth_rate_1 = str( growth_rate_1[0] )
			except ValueError: # this error means that N1_ancient has no prior range but is set relative to N1 (identical, smaller or larger)
				if prior_args["N1_ancient"] == "N1_recent":# to be identical to N1_recent => makes a single-phase model!
					N1_ancient = N1_recent
					growth_rate_1 = "0.0"
				elif prior_args["N1_ancient"] == "larger_than_recent":
					N1_ancient = numpy.random.uniform( N1_recent, prior_args["N1_recent"][1], 1 )[0]
					growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent)
					growth_rate_1 = str( growth_rate_1[0] )
				elif prior_args["N1_ancient"] == "smaller_than_recent":
					N1_ancient = numpy.random.uniform( prior_args["N1_recent"][0], N1_recent, 1 )[0]
					growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent)
					growth_rate_1 = str( growth_rate_1[0] )

			try: # the special flag is caught here for N2_ancient:
				if prior_args["N_priors_shape"] == "flat":
					N2_ancient = numpy.random.uniform( prior_args["N2_ancient"][0], prior_args["N2_ancient"][1], 1 )
				elif prior_args["N_priors_shape"] == "lognormal":
					N2_ancient = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N2_ancient"][0], prior_args["N2_ancient"][1] )
				growth_rate_2 = -(1.0/T_N2_ancient)*numpy.log(N2_ancient/N2_recent)
				growth_rate_2 = str( growth_rate_2[0] )	
			except ValueError:  # this error means that N2_ancient has no prior range but is set relative to N2 (identical, smaller or larger)
				if prior_args["N2_ancient"] == "N2_recent":# to be identical to N2_recent => makes a single-phase model!
					N2_ancient = N2_recent
					growth_rate_2 = "0.0"
				elif prior_args["N2_ancient"] == "larger_than_recent":
					N2_ancient = numpy.random.uniform( N2_recent, prior_args["N2_recent"][1], 1 )[0]
					growth_rate_2 = -(1.0/T_N2_ancient)*numpy.log(N2_ancient/N2_recent)
					growth_rate_2 = str( growth_rate_2[0] )	
				elif prior_args["N2_ancient"] == "smaller_than_recent":
					N2_ancient = numpy.random.uniform( prior_args["N2_recent"][0], N2_recent, 1 )[0]
					growth_rate_2 = -(1.0/T_N2_ancient)*numpy.log(N2_ancient/N2_recent)
					growth_rate_2 = str( growth_rate_2[0] )	

			# now the locus-specific stuff:
			# get the migration rates first, this may be clever for speed if all the draws are done before the stdout loop, but have to test that
			nloci = len( loci_properties )			
			
			if migration_shape == "homogenous":
				m12_list = numpy.array( [1]*nloci ) * m12_scale
				m21_list = numpy.array( [1]*nloci ) * m21_scale
			
				m12beta1 = m12beta2 = m21beta1 = m21beta2 = "NA"
				m12_prop_mig = m21_prop_mig = "NA"
				
				
			elif migration_shape == "heterogenous_beta":
				
				# hyperpriors:
				m12beta1 = numpy.random.uniform( prior_args["m12beta1"][0], prior_args["m12beta1"][1], 1 )
				m12beta2 = numpy.random.uniform( prior_args["m12beta2"][0], prior_args["m12beta2"][1], 1 )
				m21beta1 = numpy.random.uniform( prior_args["m21beta1"][0], prior_args["m21beta1"][1], 1 )
				m21beta2 = numpy.random.uniform( prior_args["m21beta2"][0], prior_args["m21beta2"][1], 1 )
			
				m12_list = numpy.random.beta( m12beta1 , m12beta2, nloci) * m12_scale		
				m21_list = numpy.random.beta( m21beta1 , m21beta2, nloci) * m21_scale

				m12_prop_mig = m21_prop_mig = "NA"
				
			elif migration_shape == "heterogenous_binomial":
				
				m12_prop_mig = numpy.random.uniform( prior_args["m12_prop_mig"][0], prior_args["m12_prop_mig"][1], 1 )
				m21_prop_mig = numpy.random.uniform( prior_args["m21_prop_mig"][0], prior_args["m21_prop_mig"][1], 1 )
				
				m12_list = numpy.random.binomial( 1 , m12_prop_mig, nloci) * m12_scale		
				m21_list = numpy.random.binomial( 1 , m21_prop_mig, nloci) * m21_scale
				
				m12beta1 = m12beta2 = m21beta1 = m21beta2 = "NA"
			
			## now N_hetero
			hetero_Ns = []
			# pre-compute the random stuff to reduce number of function calls
			rescales = numpy.random.beta( N_hetero_beta1 , N_hetero_beta2, 3*nloci) 
			resc_idx = 0
			N1_is_het = numpy.random.binomial( 1 , N1_prop_hetero, nloci) 
			N2_is_het = numpy.random.binomial( 1 , N2_prop_hetero, nloci) 
			Nanc12_is_het = numpy.random.binomial( 1 , Nanc12_prop_hetero, nloci) 
			for idx in range( nloci ):
				loc_Ns = []
				if N1_is_het[idx] == 1:
					loc_Ns += [N1_recent*rescales[resc_idx], N1_ancient*rescales[resc_idx]]
				else:
					loc_Ns += [N1_recent, N1_ancient]
				if N2_is_het[idx] == 1:
					loc_Ns += [N2_recent*rescales[resc_idx+1], N2_ancient*rescales[resc_idx+1]]
				else:
					loc_Ns += [N2_recent, N2_ancient]
				if Nanc12_is_het[idx] == 1:
					loc_Ns += [Nanc12*rescales[resc_idx+2]]
				else:
					loc_Ns += [Nanc12]
				hetero_Ns.append( [str(round_except_TypeError(x, 5)) for x in loc_Ns ] )
				resc_idx += 3
						
#			print hetero_Ns


			parameters = [str(round_except_TypeError(x, 5)) for x in [		
			N1_recent, N2_recent, N1_ancient, N2_ancient, Nanc12, T_N1_ancient, T_N2_ancient , T_merge_2_into_1, T_m12_start,T_m12_stop,T_m21_start,T_m21_stop, m12_scale,  m21_scale, m12_prop_mig, m21_prop_mig, m12beta1, m12beta2, m21beta1, m21beta2, N1_prop_hetero, N2_prop_hetero, Nanc12_prop_hetero, N_hetero_beta1, N_hetero_beta2] ]
			
			[ N1_recent, N2_recent, N1_ancient, N2_ancient, Nanc12, T_N1_ancient, T_N2_ancient , T_merge_2_into_1, T_m12_start,T_m12_stop,T_m21_start,T_m21_stop, m12_scale,  m21_scale, m12_prop_mig, m21_prop_mig, m12beta1, m12beta2, m21beta1, m21beta2, N1_prop_hetero, N2_prop_hetero, Nanc12_prop_hetero, N_hetero_beta1, N_hetero_beta2] = parameters 			

			# now the locus-specific stuff:
			migrates = []			
			for idx in range( nloci ):

				migrates.append( [str(round_except_TypeError(x, 5)) for x in [ m12_list[idx], m21_list[idx] ]] )
						
			for idx, locus in enumerate(loci_properties):

				[locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, total_nsam_per_locus], [m12, m21], [N1_recent, N1_ancient, N2_recent, N2_ancient, Nanc12] = loci_properties[idx], migrates[idx], hetero_Ns[idx]  	

				# writes directly to stdout -> output can be piped to msnsam 
				sys.stdout.write(
				" ".join( [total_nsam_per_locus, locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, 
				N1_recent, growth_rate_1, 
				N2_recent, growth_rate_2, 
				T_N1_ancient, N1_ancient, 
				T_N2_ancient, N2_ancient, 
				T_m12_start, m12, 
				T_m12_stop,
				T_m21_start, m21 ,
				T_m21_stop,
				T_merge_2_into_1,
				T_merge_2_into_1, Nanc12
				]) + "\n")
	
			
			# write prior realisations to parameters.txt log file
			PARAMETERS_OUT.write( "\t".join( parameters  ) + "\n" )

				


###
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
		loci_properties.append( [locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, str(int(nsam_1)+int(nsam_2)) ] )
	
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
						print "migration_shape must be one of 'homogenous', 'heterogenous_beta' or 'heterogenous_binomial'"
						exit()
				elif parameter == "N_priors_shape":
					if entry in ["flat", "lognormal"]:
						prior_args[parameter] = entry
					else:
						print "N_priors_shape must be one of 'flat' or 'lognormal'"
						exit()

				else:	
					try:		
						prior_args[parameter] = [ float(x) for x in entry.split(",") ]
					except ValueError:
						if parameter == "T_m12_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						elif parameter == "T_m21_stop" and entry == "T_merge_2_into_1":
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
						else:							
							print "argfile format is offended by ", parameter
							exit()

	if nreps == None:
		print "argfile is missing specification of nreps"
		exit()
	
	# EDIT this:		
	params_to_be_read = [ "N1_recent", "N2_recent", "N1_ancient", "N2_ancient", "Nanc12", "T_N1_ancient", "T_N2_ancient" ,
	"T_merge_2_into_1",
	"T_m12_start","T_m12_stop",
"T_m21_start",
"T_m21_stop",
"m12_scale",
"m21_scale",
"migration_shape",
"m12_prop_mig",
"m21_prop_mig",
"m12beta1",
"m12beta2",
"m21beta1",
"m21beta2",
"N_priors_shape","N_mode_of_lognormal","N_sigma",
"N_prop_hetero",
"N_hetero_beta1",
"N_hetero_beta2"
]

	for param in params_to_be_read:
		if param not in prior_args.keys():
			print "argfile is missing specification of ", param
			exit()
			
	return nreps, prior_args


####################			
### MAIN


args = 	get_commandline_arguments ()

loci_properties = read_bpfile(args.bpfile)
#print loci_properties

n_iterations, prior_args = read_argfile(args.argfile)
#print n_iterations, prior_args

draw_and_pipe_priors (n_iterations, loci_properties, prior_args)







