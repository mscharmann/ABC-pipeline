#!/usr/local/bin/python
# Python 3
# 2020
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
	
	loci_properties_unmod = loci_properties
	
	# EDIT this:
	param_log_header = ["N1_recent", "N1_interm", "N1_ancient","T_N1_interm", "T_N1_ancient"]
	
#	print(relative_N_distr)
	
	with open("parameters.txt", "w") as PARAMETERS_OUT:
		PARAMETERS_OUT.write( "\t".join( param_log_header ) + "\n" )
		
		for i in range(n_iterations) : # argfile specifies how many sets!
			
			# draw independent priors (some pop sizes, scales of migration rates, N_prop_hetero, shapes of N_hetero distribution)
			if prior_args["N_priors_shape"] == "flat":
				N1_recent = numpy.random.uniform( prior_args["N1_recent"][0], prior_args["N1_recent"][1] )
			elif prior_args["N_priors_shape"] == "lognormal":
				N1_recent = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_recent"][0], prior_args["N1_recent"][1] )
									
			T_N1_interm = numpy.random.uniform( prior_args["T_N1_interm"][0], prior_args["T_N1_interm"][1] )
			T_N1_ancient = draw_from_mincutoff_flat_range( prior_args["T_N1_ancient"][0], prior_args["T_N1_ancient"][1], T_N1_interm )
			
			## get the N_ancient sizes and catch special flag which will also affect growth rate!! (has to be defined as zero); constant size means growth rate zero!
			# the exponential growth rate parameter, which depends on N1_recent, N1_ancient and T_N1_ancient. This is NOT A FREE PARAMETER but necessary to supply it to ms to run!
			# end_size = start_size * exp^-growth_rate*end_time ### end_size and end_time are both BACKWARD in time!! 
			# thus growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent) ## natural logarithm
			# 						=> since N_hetero is implemented such that both N_ancient and N_recent are re-scaled by the same factor 
			#							(i.e 'linked selection' is constant in the two Ne phases),  
			#							the ratio N_ancient/N_recent remains the same,  and thus the growth_rate also does not need to be locus-specific.
			try: # the special flag is caught here for N1_interm:
				if prior_args["N_priors_shape"] == "flat":
					N1_interm = numpy.random.uniform( prior_args["N1_interm"][0], prior_args["N1_interm"][1] )
				elif prior_args["N_priors_shape"] == "lognormal":
					N1_interm = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_interm"][0], prior_args["N1_interm"][1] )
				
				growth_rate_1 = -(1.0/T_N1_interm)*numpy.log(N1_interm/N1_recent)
				growth_rate_1 = str( growth_rate_1 )
			except ValueError: # this error means that N1_ancient has no prior range but is set relative to N1 (identical, smaller or larger)
				if prior_args["N1_interm"] == "N1_recent":# to be identical to N1_recent => makes a single-phase model!
					N1_interm = N1_recent
					growth_rate_1 = "0.0"
				elif prior_args["N1_interm"] == "larger_than_recent":
					N1_interm = numpy.random.uniform( N1_recent, prior_args["N1_recent"][1] )
					growth_rate_1 = -(1.0/T_N1_interm)*numpy.log(N1_interm/N1_recent)
					growth_rate_1 = str( growth_rate_1 )
				elif prior_args["N1_interm"] == "smaller_than_recent":
					N1_interm = numpy.random.uniform( prior_args["N1_recent"][0], N1_recent )
					growth_rate_1 = -(1.0/T_N1_interm)*numpy.log(N1_interm/N1_recent)
					growth_rate_1 = str( growth_rate_1 )

			try: # the special flag is caught here for N1_ancient:
				if prior_args["N_priors_shape"] == "flat":
					N1_ancient = numpy.random.uniform( prior_args["N1_ancient"][0], prior_args["N1_ancient"][1] )
				elif prior_args["N_priors_shape"] == "lognormal":
					N1_ancient = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_ancient"][0], prior_args["N1_ancient"][1] )

			except ValueError:  # this error means that N1_ancient has no prior range but is set relative to N1_interm (identical, smaller or larger)
				if prior_args["N1_ancient"] == "N1_interm":# to be identical to N1_interm => makes a single-phase model!
					N1_ancient = N1_interm
				elif prior_args["N1_ancient"] == "larger_than_interm":
					try: # if can not get numerical value from N1_interm prior value, get it from N1_recent!
						N1_ancient = numpy.random.uniform( N1_interm, prior_args["N1_interm"][1])
					except (ValueError, TypeError) as e:
						N1_ancient = numpy.random.uniform( N1_interm, prior_args["N1_recent"][1])
				elif prior_args["N1_ancient"] == "smaller_than_interm":
					try: # if can not get numerical value from N1_interm prior value, get it from N1_recent!
						N1_ancient = numpy.random.uniform( prior_args["N1_interm"][0], N1_interm )
					except ValueError:
						N1_ancient = numpy.random.uniform(  prior_args["N1_recent"][0], N1_interm )

			# now the locus-specific stuff:
			nloci = len( loci_properties )				
			## now N_hetero
			relative_Ns = numpy.random.choice(relative_N_distr, size=nloci, replace=True, p=None)
#			print (numpy.mean(relative_Ns))
			hetero_Ns = []
			for idx in range( nloci ):
				loc_Ns = [x*relative_Ns[idx] for x in [N1_recent,N1_interm,N1_ancient]]
#				print(idx,relative_Ns[idx],[N1_recent,N1_ancient,N2_recent,N2_ancient,Nanc12],loc_Ns)
				hetero_Ns.append( [str(round_except_TypeError(x, 5)) for x in loc_Ns ] )
#			print (hetero_Ns)
						

			parameters = [str(round_except_TypeError(x, 5)) for x in [		
			N1_recent, N1_interm, N1_ancient, T_N1_interm, T_N1_ancient] ]
			
			[N1_recent, N1_interm, N1_ancient, T_N1_interm, T_N1_ancient] = parameters 			
						
			for idx, locus in enumerate(loci_properties):
#				print(loci_properties[idx], hetero_Ns[idx] )
				
				[locus_theta, locus_rho, locus_nsites, nsam_1, total_nsam_per_locus], [N1_recent, N1_interm, N1_ancient] = loci_properties[idx], hetero_Ns[idx]  	
				
				# writes directly to stdout -> output can be piped to msnsam 
				sys.stdout.write(
				" ".join( [total_nsam_per_locus, locus_theta, locus_rho, locus_nsites, nsam_1,
				N1_recent,
				growth_rate_1, 
				T_N1_interm,
				N1_interm, 
				T_N1_ancient,
				N1_ancient 
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
	
	## for 1 population! => replacing nsam1 with zero!!
	
	with open(bpfile, "r") as INFILE:
		header = INFILE.readline() # remove first line
		locilengthline = INFILE.readline().strip("\n").split("\t")
		nsam1_line = INFILE.readline().strip("\n").split("\t")
		thetaline = INFILE.readline().strip("\n").split("\t")
		rholine = INFILE.readline().strip("\n").split("\t")
	
	# [locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, total_nsam_per_locus]

	loci_properties = []
	for i in range(len(locilengthline)):
		locus_theta = thetaline[i]		
		locus_rho = rholine[i]
		locus_nsites = locilengthline[i]
		nsam_1 = nsam1_line[i]
		loci_properties.append( [locus_theta, locus_rho, locus_nsites, nsam_1, nsam_1 ] )
	
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
						if parameter == "N1_ancient" and entry == "N1_interm":
							prior_args[parameter] = entry
						elif parameter == "N1_ancient" and entry == "larger_than_interm":
							prior_args[parameter] = entry
						elif parameter == "N1_ancient" and entry == "smaller_than_interm":
							prior_args[parameter] = entry
						elif parameter == "N1_interm" and entry == "N1_recent":
							prior_args[parameter] = entry
						elif parameter == "N1_interm" and entry == "larger_than_recent":
							prior_args[parameter] = entry
						elif parameter == "N1_interm" and entry == "smaller_than_recent":
							prior_args[parameter] = entry

						else:							
							print ("argfile format is offended by ", parameter)
							exit()

	if nreps == None:
		print ("argfile is missing specification of nreps")
		exit()
	
	# EDIT this:		
	params_to_be_read = ["N1_recent", "N1_interm", "N1_ancient","T_N1_interm", "T_N1_ancient"]

	for param in params_to_be_read:
		if param not in prior_args.keys():
			print ("argfile is missing specification of ", param)
			exit()
			
	return nreps, prior_args


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


args = 	get_commandline_arguments ()

loci_properties = read_bpfile(args.bpfile)
#print loci_properties

relative_N_distr = read_relative_N(args.relative_N)

n_iterations, prior_args = read_argfile(args.argfile)
#print n_iterations, prior_args

draw_and_pipe_priors (n_iterations, loci_properties, prior_args, relative_N_distr)







