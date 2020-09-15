#!/usr/local/bin/python
# Python 2.7.6
# 11 Feb 2016
# 

# Mathias Scharmann

# a common issue:
# IOError: [Errno 32] Broken pipe
# is caused by mis-specification of $msnsam_arg2!

"""
generalised msnsam tbs string (paste to a spreadsheet for reading):
#
total number of samples	theta	rho and nsites	how many samples?	pop size	pop size	pop size	stop migration	stop migration	start migration	start migration	stop migration	stop migration	start migration	start migration	merge 2 and 1	set size of ancestor12	set all migration rates to zero	merge into Nanc123	set size of ancestor123
msnsam tbs $msnsam_arg2	-t tbs	-r tbs tbs 	-I 3 tbs tbs tbs 0 	-n 1 tbs 	-n 2 tbs 	-n 3 tbs 	-em tbs 1 2 0 	-em tbs 2 1 0 	-em tbs 1 2 tbs 	-em tbs 2 1 tbs 	-em tbs 2 3 0 	-em tbs 3 2 0	-em tbs 2 3 tbs 	-em tbs 3 2 tbs	-ej tbs 2 1	-en tbs 1 tbs	-eM tbs 0	-ej tbs 1 3	-en tbs 3 tbs
																			
total_nsam_per_locus	locus_theta	 locus_rho,  locus_nsites	 nsam_1,  nsam_2,  nsam_3	 N1	 N2	 N3	 T_m12_stop	 T_m12_stop	 T_m12_start,  m12	 T_m12_start,  m21	 T_m23_stop	 T_m23_stop	 T_m23_start,  m23	 T_m23_start,  m32	 T_merge_2_into_1	T_merge_2_into_1, Nanc12	 T_merge_2_into_1	 T_merge_1_into_3	 T_merge_1_into_3, Nanc123

#############

## this is: bash_wrapper.sim3pop.sh
## for EULER

ori_dir=..
bpfile=final_filtered.allsites.pres0.9.minDP6.vcf.genepop.txt.fixedtheta.bpfile.txt
spinput=final_filtered.allsites.pres0.9.minDP6.vcf.genepop.txt.spinput.txt


toolbase=/cluster/home/schamath/coalescent_simulations/ABC_tools


# common msnsam string to all models (drop of parameters is achieved by setting them zero or equal to a redundant one)
msnsam_model_cmd="-t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 2 3 0 -em tbs 3 2 0 -em tbs 2 3 tbs -em tbs 3 2 tbs -ej tbs 2 1 -en tbs 1 tbs -eM tbs 0 -ej tbs 1 3 -en tbs 3 tbs"



#for scenario in allo_heterogenous allo_homogenous island_heterogenous_and_ghost island_heterogenous island_homogenous_and_ghost island_homogenous iso sym_heterogenous sym_homogenous ;
#for scenario in iso island_hom island_hetbin sym_hom sym_hetbin allo_hom allo_hetbin ;
scenario=3pop

argfile=argfile.${scenario}.txt
	
module load python/2.7
		
mkdir temp_${scenario}_${LSB_JOBINDEX}
cd temp_${scenario}_${LSB_JOBINDEX}

cp ${ori_dir}/${bpfile} ./
cp ${ori_dir}/${argfile} ./
cp ${ori_dir}/${spinput} ./

mv $spinput spinput.txt

# build msnsam_arg2 from the arg and spinput files (nloci * nreps); cannot be a tbs argument to msnsam!:
nloci=$(head -2 spinput.txt | tail -1 )
nreps=$( cat $argfile | grep "nreps" | awk 'FS=" = " {print $3}' )
msnsam_arg2=$(( $nloci * $nreps ))

# make sure spinput also contains correct nreps from argfile: this is actually only relevant for correct output of the progress report:
head -n -3 spinput.txt > firstpart
tail -2 spinput.txt > lastpart
cat firstpart <( echo $nreps ) lastpart > spinput.txt
rm firstpart lastpart

python $toolbase/draw_ms_priors.3pop.py -bpfile $bpfile -argfile $argfile  | $toolbase/msnsam tbs $msnsam_arg2 $msnsam_model_cmd | python $toolbase/ms2stats.3pop.counter.py | python $toolbase/ms2stats.3pop.stats.2DSFS.py


"""

# example:
# pseudocode:

#### HEAD


import sys
import os
import argparse
import numpy

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

	
def draw_from_maxcutoff_flat_range( lowerbound, upperbound, maxcutoff ):
	
	x = maxcutoff + 1
	while x > maxcutoff:
		x = numpy.random.uniform( lowerbound, upperbound, 1 )
	
	return x


def draw_and_pipe_priors_3pops (n_iterations, loci_properties, prior_args):
	
	migration_shape = prior_args["migration_shape"]
	
	loci_properties_unmod = loci_properties
	
	with open("parameters.txt", "w") as PARAMETERS_OUT:
		PARAMETERS_OUT.write( "\t".join( [ "N1","N2","N3","Nanc12","Nanc123","T_m12_start","T_m12_stop","T_m23_start","T_m23_stop","T_merge_2_into_1","T_merge_1_into_3","m12_scale","m21_scale","m23_scale","m32_scale","m12_prop_mig","m21_prop_mig","m23_prop_mig","m32_prop_mig","m12beta1","m12beta2","m21beta1","m21beta2","m23beta1","m23beta2","m32beta1","m32beta2" ] ) + "\n" )		
		for i in range(n_iterations) : # argfile specifies how many sets!
			
			# draw independent priors:
			N1 = numpy.random.uniform( prior_args["N1"][0], prior_args["N1"][1], 1 )
			N2 = numpy.random.uniform( prior_args["N2"][0], prior_args["N2"][1], 1 )
			N3 = numpy.random.uniform( prior_args["N3"][0], prior_args["N3"][1], 1 )
			Nanc12 = numpy.random.uniform( prior_args["Nanc12"][0], prior_args["Nanc12"][1], 1 )
			Nanc123 = numpy.random.uniform( prior_args["Nanc123"][0], prior_args["Nanc123"][1], 1 )	
			
			m12_scale = numpy.random.uniform( prior_args["m12_scale"][0], prior_args["m12_scale"][1], 1 )
			m21_scale = numpy.random.uniform( prior_args["m21_scale"][0], prior_args["m21_scale"][1], 1 )
			m23_scale = numpy.random.uniform( prior_args["m23_scale"][0], prior_args["m23_scale"][1], 1 )
			m32_scale = numpy.random.uniform( prior_args["m32_scale"][0], prior_args["m32_scale"][1], 1 )
	
			# draw interdependent priors (Time-priors)
			T_merge_1_into_3 = numpy.random.uniform( prior_args["T_merge_1_into_3"][0], prior_args["T_merge_1_into_3"][1], 1 )
			T_merge_2_into_1 = draw_from_maxcutoff_flat_range( prior_args["T_merge_2_into_1"][0],  prior_args["T_merge_2_into_1"][1], T_merge_1_into_3 )
						
			try: # for the general case where T_m12_stop is free and only has upper constraint = T_merge_2_into_1
				T_m12_stop = draw_from_maxcutoff_flat_range( prior_args["T_m12_stop"][0], prior_args["T_m12_stop"][1], T_merge_2_into_1 )
			except ValueError: # for the special case where ghost migration always starts at speciation of 1and2 :
				T_m12_stop = T_merge_2_into_1
			T_m12_start = draw_from_maxcutoff_flat_range( prior_args["T_m12_start"][0], prior_args["T_m12_start"][1], T_m12_stop )
			
			try: # for the general case where T_m23_stop is free and only has upper constraint = T_merge_2_into_1
				T_m23_stop = draw_from_maxcutoff_flat_range( prior_args["T_m23_stop"][0], prior_args["T_m23_stop"][1], T_merge_2_into_1 )
			except ValueError: # this error means that T_m23_stop has nor prior range but is set to be identical to T_merge_2_into_1!
				T_m23_stop = T_merge_2_into_1
			T_m23_start = draw_from_maxcutoff_flat_range( prior_args["T_m23_start"][0], prior_args["T_m23_start"][1], T_m23_stop) 
			
			# now the locus-specific stuff:
			# get the migration rates first, this may be clever for speed if all the draws are done before the stdout loop, but have to test that
			nloci = len( loci_properties )			
						
			if migration_shape == "homogenous":
				m12_list = numpy.array( [1]*nloci ) * m12_scale
				m21_list = numpy.array( [1]*nloci ) * m21_scale
				m23_list = numpy.array( [1]*nloci ) * m23_scale
				m32_list = numpy.array( [1]*nloci ) * m32_scale
				
				m12beta1 = m12beta2 = m21beta1 = m21beta2 = m23beta1 = m23beta2 = m32beta1 = m32beta2 = "NA"
				m12_prop_mig = m21_prop_mig = m23_prop_mig = m32_prop_mig = "NA"
				
			elif migration_shape == "heterogenous_beta":
				
				# hyperpriors:
				m12beta1 = numpy.random.uniform( prior_args["m12beta1"][0], prior_args["m12beta1"][1], 1 )
				m12beta2 = numpy.random.uniform( prior_args["m12beta2"][0], prior_args["m12beta2"][1], 1 )
				m21beta1 = numpy.random.uniform( prior_args["m21beta1"][0], prior_args["m21beta1"][1], 1 )
				m21beta2 = numpy.random.uniform( prior_args["m21beta2"][0], prior_args["m21beta2"][1], 1 )
				m23beta1 = numpy.random.uniform( prior_args["m23beta1"][0], prior_args["m23beta1"][1], 1 )
				m23beta2 = numpy.random.uniform( prior_args["m23beta2"][0], prior_args["m23beta2"][1], 1 )
				m32beta1 = numpy.random.uniform( prior_args["m32beta1"][0], prior_args["m32beta1"][1], 1 )
				m32beta2 = numpy.random.uniform( prior_args["m32beta2"][0], prior_args["m32beta2"][1], 1 )
			
				m12_list = numpy.random.beta( m12beta1 , m12beta2, nloci) * m12_scale		
				m21_list = numpy.random.beta( m21beta1 , m21beta2, nloci) * m21_scale
				m23_list = numpy.random.beta( m23beta1 , m23beta2, nloci) * m23_scale		
				m32_list = numpy.random.beta( m32beta1 , m32beta2, nloci) * m32_scale
				
				m12_prop_mig = m21_prop_mig = m23_prop_mig = m32_prop_mig = "NA"
				
			elif migration_shape == "heterogenous_binomial":
				
				m12_prop_mig = numpy.random.uniform( prior_args["m12_prop_mig"][0], prior_args["m12_prop_mig"][1], 1 )
				m21_prop_mig = numpy.random.uniform( prior_args["m21_prop_mig"][0], prior_args["m21_prop_mig"][1], 1 )
				m23_prop_mig = numpy.random.uniform( prior_args["m23_prop_mig"][0], prior_args["m23_prop_mig"][1], 1 )
				m32_prop_mig = numpy.random.uniform( prior_args["m32_prop_mig"][0], prior_args["m32_prop_mig"][1], 1 )
				
				m12_list = numpy.random.binomial( 1 , m12_prop_mig, nloci) * m12_scale		
				m21_list = numpy.random.binomial( 1 , m21_prop_mig, nloci) * m21_scale
				m23_list = numpy.random.binomial( 1 , m23_prop_mig, nloci) * m23_scale		
				m32_list = numpy.random.binomial( 1 , m32_prop_mig, nloci) * m32_scale
				
				m12beta1 = m12beta2 = m21beta1 = m21beta2 = m23beta1 = m23beta2 = m32beta1 = m32beta2 = "NA"

			parameters = [str(round_except_TypeError(x, 5)) for x in [N1,N2,N3,Nanc12,Nanc123,T_m12_start,T_m12_stop,T_m23_start,T_m23_stop,T_merge_2_into_1,T_merge_1_into_3,m12_scale,m21_scale,m23_scale,m32_scale,m12_prop_mig,m21_prop_mig,m23_prop_mig,m32_prop_mig,m12beta1,m12beta2,m21beta1,m21beta2,m23beta1,m23beta2,m32beta1,m32beta2] ]

			[N1,N2,N3,Nanc12,Nanc123,T_m12_start,T_m12_stop,T_m23_start,T_m23_stop,T_merge_2_into_1,T_merge_1_into_3,m12_scale,m21_scale,m23_scale,m32_scale,m12_prop_mig,m21_prop_mig,m23_prop_mig,m32_prop_mig,m12beta1,m12beta2,m21beta1,m21beta2,m23beta1,m23beta2,m32beta1,m32beta2] = parameters 				

			# write prior realisations to parameters.txt
			PARAMETERS_OUT.write( "\t".join( parameters  ) + "\n" )
						
			# now the locus-specific stuff:
			migrates = []			
			for idx in range( nloci ):

				migrates.append( [str(round_except_TypeError(x, 5)) for x in [ m12_list[idx], m21_list[idx], m23_list[idx], m32_list[idx] ]] )
		
			for idx, locus in enumerate(loci_properties):

				[locus_theta, locus_rho, locus_nsites, nsam_1, nsam_2, nsam_3, total_nsam_per_locus], [m12, m21, m23, m32] = loci_properties[idx], migrates[idx] 	

				# writes directly to stdout -> output can be piped to msnsam 
				sys.stdout.write(
				" ".join( [total_nsam_per_locus, locus_theta, locus_rho, locus_nsites,	nsam_1, nsam_2, nsam_3,	N1,	N2, N3, T_m12_stop, T_m12_stop, T_m12_start, m12, T_m12_start, m21, T_m23_stop, T_m23_stop, T_m23_start, m23, T_m23_start, m32, T_merge_2_into_1, T_merge_2_into_1, Nanc12, T_merge_2_into_1, T_merge_1_into_3, T_merge_1_into_3, Nanc123] ) + "\n" )

##



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

				else:	
					try:		
						prior_args[parameter] = [ float(x) for x in entry.split(",") ]
					except ValueError:
						if parameter == "T_m12_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						elif parameter == "T_m23_stop" and entry == "T_merge_2_into_1":
							prior_args[parameter] = entry
						else:							
							print "argfile format is offended by ", parameter
							exit()

	if nreps == None:
		print "argfile is missing specification of nreps"
		exit()
	
	for i in ["m12beta1","m12beta2","m12_prop_mig","m12_scale","m21beta1","m21beta2","m21_prop_mig","m21_scale","m23beta1","m23beta2","m23_prop_mig","m23_scale","m32beta1","m32beta2","m32_prop_mig","m32_scale","migration_shape","N1","N2","N3","Nanc123","Nanc12","T_m12_start","T_m12_stop","T_m23_start","T_m23_stop","T_merge_1_into_3","T_merge_2_into_1"]:
		if i not in prior_args.keys():
			print "argfile is missing specification of ", i
			exit()
			
	return nreps, prior_args
	

####################			
### MAIN

args = 	get_commandline_arguments ()

loci_properties = read_bpfile(args.bpfile)
#print loci_properties

n_iterations, prior_args = read_argfile(args.argfile)
#print n_iterations, prior_args

draw_and_pipe_priors_3pops (n_iterations, loci_properties, prior_args)








