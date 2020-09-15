#!/usr/local/bin/python
# Python 2.7.6
# 12 Feb 2018
# 


# Mathias Scharmann
# added feature: for twophase models, Nx_ancient can be specified as LARGER or SMALLER
"CAREFUL: T_m_stop NOT YET READY FOR EMPIRICAL DISTRIBUTIONS!!!"	

# a common issue:
# IOError: [Errno 32] Broken pipe
# is caused by mis-specification of $msnsam_arg2!

"""
tb=$HOME/ABC_GWH/ABC_tools/

# the generalised msnsam tbs string
$tb/msnsam/msnsam tbs $msnsam_arg2 -t tbs -r tbs tbs -I 3 tbs tbs 0 0 \ # locus-specific (fixed) parameters theta rho and nsites; sampling parameters (3 pops with x,y,0 chromosomes from 1,2,ghost)  
-n 1 tbs -n 2 tbs -n 3 tbs \ 			# Ns for all pops in most recent phase
-em tbs 1 2 0 -em tbs 2 1 0 \ 			# at T_m_stop, migration rates m12 and m21 are set to zero (for mig up to present, set T_m_stop = 0)
-em tbs 1 2 tbs -em tbs 2 1 tbs \ 		# at T_m_start == T_merge_2_into_1, migration rates m12 and m21 are set to some value (for strict isolation scenario, set m = 0)
-em tbs 1 3 0 -em tbs 2 3 0 \ 			# at T_m_ghost_stop, migration rates from ghost are set to zero (for ghostmig up to present, set T_m_ghost_stop = 0)
-em tbs 1 3 tbs -em tbs 2 3 tbs \ 		# at T_m_ghost_start == T_merge_2_into_1, migration rates from ghost are set to some value (for no ghost-scenario, set m = 0)
-ej tbs 2 1 -en tbs 1 tbs -eM tbs 0 \ 		# at T_merge_2_into_1, merge 2 into 1; at same time set N1 to some value, at same time set all migration to zero	
-ej tbs 1 3 -en tbs 3 tbs			# at T_merge_1_into_ghost, merge 1 into 3 (ghost); at same time set N3 = N (common ancestor) to some value

msnsam_arg2=

$tb/msnsam/msnsam tbs $msnsam_arg2 -t tbs -r tbs tbs -I 3 tbs tbs 0 0 -n 1 tbs -n 2 tbs -n 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 1 3 0 -em tbs 2 3 0 -em tbs 1 3 tbs -em tbs 2 3 tbs -ej tbs 2 1 -en tbs 1 tbs -eM tbs 0 -ej tbs 1 3 -en tbs 3 tbs

## example:
# make asure spinput is there and correct!

rm parameters.txt ABCstat.txt myfifo

tb=$HOME/ABC_GWH/ABC_tools/
msnsam_arg2=40900
msnsam_model_cmd="-t tbs -r tbs tbs -I 3 tbs tbs 0 0 -n 1 tbs -n 2 tbs -n 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 1 3 0 -em tbs 2 3 0 -em tbs 1 3 tbs -em tbs 2 3 tbs -ej tbs 2 1 -en tbs 1 tbs -eM tbs 0 -ej tbs 1 3 -en tbs 3 tbs"

mkfifo myfifo
$tb/mscalc/mscalc <myfifo &

python draw_ms_priors.py -bpfile final_filtered.allsites.pres0.9.minDP10.vcf.genepop.txt.bpfile.txt -argfile argfile_ghostswitch.txt  | $tb/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd >myfifo

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
	parser.add_argument("-emp_distr", help="distributions file: contains empricial distributions of parameter values", type=extant_file, required=True, metavar="FILE")
	
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



def draw_and_pipe_priors_and_emp_distr_ghost_off (n_iterations, loci_properties, prior_args, empirical_distributions):
	
	"""
	the empirical distributions have to override the priors, if they exist. If not, then sample from prior!!
	numpy.random.choice( emp_dr )
	"""
	for empdr in empirical_distributions.keys():
		if empdr not in prior_args.keys():
			print empdr, " is not a recognised parameter empirical distribution, dying"
			exit()
		else:
			prior_args[empdr] = "empirical"
		
	migration_shape = prior_args["migration_shape"]
	
	loci_properties_unmod = loci_properties
	
	with open("parameters.txt", "w") as PARAMETERS_OUT:
		PARAMETERS_OUT.write( "\t".join( ["N1_recent", "N2_recent", "N1_ancient", "N2_ancient", "Nghost", "Nanc12", "Nanc12ghost", "T_merge_2_into_1", "T_merge_1_into_ghost", "T_N1_ancient", "T_N2_ancient", "T_m_start", "T_m_stop", "T_m_ghost_start", "T_m_ghost_stop", "m12_scale", "m21_scale", "m1ghost_scale", "m2ghost_scale", "m12beta1", "m12beta2", "m21beta1", "m21beta2", "m1ghostbeta1", "m1ghostbeta2", "m2ghostbeta1", "m2ghostbeta2","m12_prop_mig", "m21_prop_mig", "m1ghost_prop_mig", "m2ghost_prop_mig"] ) + "\n" )
		
		for i in range(n_iterations) : # argfile specifies how many sets!
			
			# draw independent priors:
			if prior_args["N_priors_shape"] == "flat":
				
				if not prior_args["N1_recent"] == "empirical":
					N1_recent = numpy.random.uniform( prior_args["N1_recent"][0], prior_args["N1_recent"][1], 1 )
				else:
					N1_recent = numpy.random.choice( empirical_distributions["N1_recent"] )
				
				if not prior_args["N2_recent"] == "empirical":
					N2_recent = numpy.random.uniform( prior_args["N2_recent"][0], prior_args["N2_recent"][1], 1 )
				else:
					N2_recent = numpy.random.choice( empirical_distributions["N2_recent"] )
				
				Nghost = 0
				
				if not prior_args["Nanc12"] == "empirical":
					Nanc12 = numpy.random.uniform( prior_args["Nanc12"][0], prior_args["Nanc12"][1], 1 )
				else:
					Nanc12 = numpy.random.choice( empirical_distributions["Nanc12"] )
					
				Nanc12ghost = Nanc12
				
			elif prior_args["N_priors_shape"] == "lognormal":

				if not prior_args["N1_recent"] == "empirical":
					N1_recent = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_recent"][0], prior_args["N1_recent"][1] )
				else:
					N1_recent = numpy.random.choice( empirical_distributions["N1_recent"] )
#				N1_recent = numpy.random.uniform( prior_args["N1_recent"][0], prior_args["N1_recent"][1], 1 )

				if not prior_args["N2_recent"] == "empirical":
					N2_recent = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N2_recent"][0], prior_args["N2_recent"][1] )
				else:
					N2_recent = numpy.random.choice( empirical_distributions["N2_recent"] )
				
				Nghost = 0
				
				if not prior_args["Nanc12"] == "empirical":
					Nanc12 = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["Nanc12"][0], prior_args["Nanc12"][1] )
				else:
					Nanc12 = numpy.random.choice( empirical_distributions["Nanc12"] )
				
#				Nanc12 = numpy.random.uniform( prior_args["Nanc12"][0], prior_args["Nanc12"][1], 1 )
				Nanc12ghost = Nanc12
				
			
			if not prior_args["m12_scale"] == "empirical":
				m12_scale = numpy.random.uniform( prior_args["m12_scale"][0], prior_args["m12_scale"][1], 1 )
			else:
				m12_scale = numpy.random.choice( empirical_distributions["m12_scale"] )
			if not prior_args["m21_scale"] == "empirical":
				m21_scale = numpy.random.uniform( prior_args["m21_scale"][0], prior_args["m21_scale"][1], 1 )
			else:
				m21_scale = numpy.random.choice( empirical_distributions["m21_scale"] )
				
			m1ghost_scale = 0
			m2ghost_scale = 0
			
				
			# draw interdependent priors (Time-priors)
			
			if not prior_args["T_merge_2_into_1"] == "empirical":
				T_merge_2_into_1 = numpy.random.uniform( prior_args["T_merge_2_into_1"][0],  prior_args["T_merge_2_into_1"][1], 1)
			else:
				T_merge_2_into_1 = numpy.random.choice( empirical_distributions["T_merge_2_into_1"] )	
			T_merge_1_into_ghost = T_merge_2_into_1
			
			T_m_ghost_stop = 0
			T_m_ghost_start = 0
			
			"CAREFUL: T_m_stop NOT YET READY FOR EMPIRICAL DISTRIBUTIONS!!!"			
			try: # the special flag is caught here:
				T_m_stop = draw_from_maxcutoff_flat_range( prior_args["T_m_stop"][0], prior_args["T_m_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m_stop = T_merge_2_into_1
					
			if not prior_args["T_m_start"] == "empirical":
				T_m_start = draw_from_maxcutoff_flat_range( prior_args["T_m_start"][0], prior_args["T_m_start"][1], T_m_stop) 
			else:
				T_m_start = numpy.random.choice( empirical_distributions["T_m_start"] )
			
			if not prior_args["T_N1_ancient"] == "empirical":
				T_N1_ancient = draw_from_maxcutoff_flat_range( prior_args["T_N1_ancient"][0], prior_args["T_N1_ancient"][1], T_merge_2_into_1)
			else:
				T_N1_ancient = numpy.random.choice( empirical_distributions["T_N1_ancient"] )
			
			if not prior_args["T_N2_ancient"] == "empirical":
				T_N2_ancient = draw_from_maxcutoff_flat_range( prior_args["T_N2_ancient"][0], prior_args["T_N2_ancient"][1], T_merge_2_into_1)
			else:
				T_N2_ancient = numpy.random.choice( empirical_distributions["T_N2_ancient"] )
			
			## get the N_ancient sizes and catch special flag which will also affect growth rate!! (has to be defined as zero); constant size means growth rate zero!
			# the exponential growth rate parameter, which depends on N1_recent, N1_ancient and T_N1_ancient. This is NOT A FREE PARAMETER but necessary to supply it to ms to run!
			# end_size = start_size * exp^-growth_rate*end_time ### end_size and end_time are both BACKWARD in time!! 
			# thus growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent) ## natural logarithm
			try: # the special flag is caught here:
				if prior_args["N_priors_shape"] == "flat":
					N1_ancient = numpy.random.uniform( prior_args["N1_ancient"][0], prior_args["N1_ancient"][1], 1 )
				elif prior_args["N_priors_shape"] == "lognormal":
					N1_ancient = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N1_ancient"][0], prior_args["N1_ancient"][1] )
				growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent)
				growth_rate_1 = str( growth_rate_1[0] )
			except ValueError: # this error means that N1_ancient has no prior range but is set relative to N1 (identical, smaller or larger) OR even should be drawn from an empirical distribution
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
				elif prior_args["N1_ancient"] == "empirical":
					N1_ancient = numpy.random.choice( empirical_distributions["N1_ancient"] )
					growth_rate_1 = -(1.0/T_N1_ancient)*numpy.log(N1_ancient/N1_recent)
					growth_rate_1 = str( growth_rate_1 )
					
			try: # the special flag is caught here:
				if prior_args["N_priors_shape"] == "flat":
					N2_ancient = numpy.random.uniform( prior_args["N2_ancient"][0], prior_args["N2_ancient"][1], 1 )
				elif prior_args["N_priors_shape"] == "lognormal":
					N2_ancient = draw_from_truncated_lognormal( prior_args["N_mode_of_lognormal"][0], prior_args["N_sigma"][0], prior_args["N2_ancient"][0], prior_args["N2_ancient"][1] )
				growth_rate_2 = -(1.0/T_N2_ancient)*numpy.log(N2_ancient/N2_recent)
				growth_rate_2 = str( growth_rate_2[0] )	
			except ValueError:  # this error means that N2_ancient has no prior range but is set relative to N2 (identical, smaller or larger) OR even should be drawn from an empirical distribution
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
				elif prior_args["N2_ancient"] == "empirical":
					N1_ancient = numpy.random.choice( empirical_distributions["N2_ancient"] )
					growth_rate_2 = -(1.0/T_N2_ancient)*numpy.log(N2_ancient/N2_recent)
					growth_rate_2 = str( growth_rate_2 )

			# now the locus-specific stuff:
			# get the migration rates first, this may be clever for speed if all the draws are done before the stdout loop, but have to test that
			nloci = len( loci_properties )			
			
			if migration_shape == "homogenous":
				m12_list = numpy.array( [1]*nloci ) * m12_scale
				m21_list = numpy.array( [1]*nloci ) * m21_scale
				m1ghost_list = numpy.array( [1]*nloci ) * m1ghost_scale
				m2ghost_list = numpy.array( [1]*nloci ) * m2ghost_scale
				
				m12beta1 = m12beta2 = m21beta1 = m21beta2 = m1ghostbeta1 = m1ghostbeta2 = m2ghostbeta1 = m2ghostbeta2 = "NA"
				m12_prop_mig = m21_prop_mig = m1ghost_prop_mig = m2ghost_prop_mig = "NA"
				
				
			elif migration_shape == "heterogenous_beta":
				
				# hyperpriors:
				m12beta1 = numpy.random.uniform( prior_args["m12beta1"][0], prior_args["m12beta1"][1], 1 )
				m12beta2 = numpy.random.uniform( prior_args["m12beta2"][0], prior_args["m12beta2"][1], 1 )
				m21beta1 = numpy.random.uniform( prior_args["m21beta1"][0], prior_args["m21beta1"][1], 1 )
				m21beta2 = numpy.random.uniform( prior_args["m21beta2"][0], prior_args["m21beta2"][1], 1 )
				m1ghostbeta1 = 1
				m1ghostbeta2 = 1e-12
				m2ghostbeta1 = 1
				m2ghostbeta2 = 1e-12
			
				m12_list = numpy.random.beta( m12beta1 , m12beta2, nloci) * m12_scale		
				m21_list = numpy.random.beta( m21beta1 , m21beta2, nloci) * m21_scale
				m1ghost_list = numpy.random.beta( m1ghostbeta1 , m1ghostbeta2, nloci) * m1ghost_scale		
				m2ghost_list = numpy.random.beta( m2ghostbeta1 , m2ghostbeta2, nloci) * m2ghost_scale
				
				m12_prop_mig = m21_prop_mig = m1ghost_prop_mig = m2ghost_prop_mig = "NA"
				
			elif migration_shape == "heterogenous_binomial":
				
				m12_prop_mig = numpy.random.uniform( prior_args["m12_prop_mig"][0], prior_args["m12_prop_mig"][1], 1 )
				m21_prop_mig = numpy.random.uniform( prior_args["m21_prop_mig"][0], prior_args["m21_prop_mig"][1], 1 )
				m1ghost_prop_mig = 0 #numpy.random.uniform( prior_args["m1ghost_prop_mig"][0], prior_args["m1ghost_prop_mig"][1], 1 )
				m2ghost_prop_mig = 0 #numpy.random.uniform( prior_args["m2ghost_prop_mig"][0], prior_args["m2ghost_prop_mig"][1], 1 )
				
				m12_list = numpy.random.binomial( 1 , m12_prop_mig, nloci) * m12_scale		
				m21_list = numpy.random.binomial( 1 , m21_prop_mig, nloci) * m21_scale
				m1ghost_list = numpy.random.binomial( 1 , m1ghost_prop_mig, nloci) * m1ghost_scale		
				m2ghost_list = numpy.random.binomial( 1 , m2ghost_prop_mig, nloci) * m2ghost_scale
				
				m12beta1 = m12beta2 = m21beta1 = m21beta2 = m1ghostbeta1 = m1ghostbeta2 = m2ghostbeta1 = m2ghostbeta2 = "NA"

			parameters = [str(round_except_TypeError(x, 5)) for x in [N1_recent, N2_recent, N1_ancient, N2_ancient, Nghost, Nanc12, Nanc12ghost, T_merge_2_into_1, T_merge_1_into_ghost, T_N1_ancient, T_N2_ancient, T_m_start, T_m_stop, T_m_ghost_start, T_m_ghost_stop, m12_scale, m21_scale, m1ghost_scale, m2ghost_scale, m12beta1, m12beta2, m21beta1, m21beta2, m1ghostbeta1, m1ghostbeta2, m2ghostbeta1, m2ghostbeta2, m12_prop_mig, m21_prop_mig, m1ghost_prop_mig, m2ghost_prop_mig] ]

			[N1_recent, N2_recent, N1_ancient, N2_ancient, Nghost, Nanc12, Nanc12ghost, T_merge_2_into_1, T_merge_1_into_ghost, T_N1_ancient, T_N2_ancient, T_m_start, T_m_stop, T_m_ghost_start, T_m_ghost_stop, m12_scale, m21_scale, m1ghost_scale, m2ghost_scale, m12beta1, m12beta2, m21beta1, m21beta2, m1ghostbeta1, m1ghostbeta2, m2ghostbeta1, m2ghostbeta2, m12_prop_mig, m21_prop_mig, m1ghost_prop_mig, m2ghost_prop_mig] = parameters 			

			# write prior realisations to parameters.txt
			PARAMETERS_OUT.write( "\t".join( parameters  ) + "\n" )
						
			# now the locus-specific stuff:
			migrates = []			
			for idx in range( nloci ):

				migrates.append( [str(round_except_TypeError(x, 5)) for x in [ m12_list[idx], m21_list[idx], m1ghost_list[idx], m2ghost_list[idx] ]] )
		
			for idx, locus in enumerate(loci_properties):

				[locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, total_nsam_per_locus], [m12, m21, m1ghost, m2ghost] = loci_properties[idx], migrates[idx] 	

				# writes directly to stdout -> output can be piped to msnsam 
#				sys.stdout.write(
#				" ".join( [total_nsam_per_locus, locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, N1_recent, N2_recent, Nghost, T_m_stop, T_m_stop, T_m_start, m12, T_m_start, m21, T_m_ghost_stop, T_m_ghost_stop, T_m_ghost_start, m1ghost, T_m_ghost_start, m2ghost, T_merge_2_into_1, T_merge_2_into_1, Nanc12, T_merge_2_into_1, T_merge_1_into_ghost, T_merge_1_into_ghost, Nanc12ghost] ) + "\n" )
								
				# new:
				sys.stdout.write(
				" ".join( [total_nsam_per_locus, locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2,
				N1_recent, 
				growth_rate_1,
				N2_recent,
				growth_rate_2,
				Nghost,
				T_N1_ancient, N1_ancient,
				T_N2_ancient, N2_ancient,
				T_m_start, m12,
				T_m_start, m21,
				T_m_stop,
				T_m_stop,
				T_m_ghost_start, m1ghost,
				T_m_ghost_start, m2ghost,
				T_m_ghost_stop,
				T_m_ghost_stop,
				T_merge_2_into_1,
				T_merge_2_into_1, Nanc12,
				T_merge_2_into_1,
				T_merge_1_into_ghost,
				T_merge_1_into_ghost, Nanc12ghost ]) + "\n")
				
## tranlated:
##msnsam_model_cmd=
"""-t tbs -r tbs tbs -I 3 tbs tbs 0 0 
-en 0.0 1 tbs  
-eg 0.0 1 tbs 
-en 0.0 2 tbs 
-eg 0.0 2 tbs
-en 0.0 3 tbs 
-en tbs 1 tbs
-en tbs 2 tbs
-em tbs 1 2 tbs
-em tbs 2 1 tbs
-em tbs 1 2 0
-em tbs 2 1 0
-em tbs 1 3 tbs
-em tbs 2 3 tbs
-em tbs 1 3 0
-em tbs 2 3 0
-ej tbs 2 1
-en tbs 1 tbs
-eM tbs 0
-ej tbs 1 3
-en tbs 3 tbs
"""



###
def round_except_TypeError(invalue, n_digits):
	
	try:
		outvalue = round(invalue, n_digits)
	except TypeError:
		outvalue = invalue
	return outvalue
	

def read_bpfile(bpfile):
	
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
	
	nreps = ghostswitch = None
	
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
				elif parameter == "ghost":
					if entry == "True":
						ghostswitch = True
					else:
						ghostswitch = False
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
						if parameter == "T_m_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						elif parameter == "T_m_ghost_stop" and entry == "T_merge_2_into_1":
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

#	[ N1, N2, Nghost, T_m_stop, T_m_start, T_m_ghost_stop, T_m_ghost_start, m12_scale, m21_scale, m1ghost_scale, m2ghost_scale, T_merge_2_into_1, Nanc12, T_merge_1_into_ghost, Nanc12ghost, m12beta1, m12beta2, m21beta1, m21beta2, m1ghostbeta1, m1ghostbeta2, m2ghostbeta1, m2ghostbeta2 ]

	if nreps == None:
		print "argfile is missing specification of nreps"
		exit()
		
	if ghostswitch == None:
		print "argfile is missing specification of ghostswitch"
		exit()
	
	for i in [ "N1_recent", "N2_recent", "N1_ancient", "N2_ancient", "T_N1_ancient", "T_N2_ancient" ,"Nghost", "T_m_stop", "T_m_start", "T_m_ghost_stop", "T_m_ghost_start", "m12_scale", "m21_scale", "m1ghost_scale", "m2ghost_scale", "T_merge_2_into_1", "Nanc12", "T_merge_1_into_ghost", "Nanc12ghost", "m12beta1", "m12beta2", "m21beta1", "m21beta2", "m1ghostbeta1", "m1ghostbeta2", "m2ghostbeta1", "m2ghostbeta2", "migration_shape", "m12_prop_mig", "m21_prop_mig", "m1ghost_prop_mig", "m2ghost_prop_mig",'N_priors_shape' ]:
		if i not in prior_args.keys():
			print "argfile is missing specification of ", i
			exit()
			
	return nreps, prior_args, ghostswitch	

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


def read_empir_distributions (emp_distr_file):
	
	distributions = []
	with open( emp_distr_file , "r") as INF:
		header = INF.readline().strip("\n").split("\t")[1:]
		for d in header:
			distributions.append([])	
		for line in INF:
			fields = line.strip("\n").split("\t")[1:]
			for idx,f in enumerate(fields):
				distributions[idx].append( round(float(f),5) )
	
	distr_dict = {}
	for idx,h in enumerate(header):
		distr_dict[h] = distributions[idx]
				
#	print "read ", len(distributions), " empirical distributions:"
#	print header
	return distr_dict	

####################			
### MAIN

args = 	get_commandline_arguments ()

loci_properties = read_bpfile(args.bpfile)
#print loci_properties

n_iterations, prior_args, ghostswitch = read_argfile(args.argfile)
#print n_iterations, prior_args

empirical_distributions = read_empir_distributions (args.emp_distr)	

#modify_spinput_txt(n_iterations)

if ghostswitch == True:
	draw_and_pipe_priors_incl_ghost (n_iterations, loci_properties, prior_args)
else:
	draw_and_pipe_priors_and_emp_distr_ghost_off (n_iterations, loci_properties, prior_args, empirical_distributions)







