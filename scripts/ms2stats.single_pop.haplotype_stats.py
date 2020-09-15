#!/usr/local/bin/python
# Python 3
# 2020
#
# ms2stats.arbitrary_npop.counter.py
#
# Mathias Scharmann

# INPUTS: 	- stdin in ms format
#			- "spinput.txt"-file: specifies number of loci, sample sizes and site-lengths for each locus

# OUTPUT: 	- "ABCstat.txt"-file with selected popgen summary statistics

# example usage:
"""
just counts ms output and sends to stdout

python $HOME/tools/draw_ms_priors.1.py -bpfile $bpfile -argfile $argfile  | $toolbase/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd | python ms2stats.counter.py | python ms2stats.stats.py


"""

##############
# HEAD
##############

import sys
import numpy as np

##############
# FUNCTIONS
##############

def determine_npop_from_spinput ():
	
	line_cnt = 0
	
	with open("spinput.txt", "r") as INFILE:
		INFILE.readline() # discard first line (empty)
		n_loci_per_replicate = int( INFILE.readline().strip("\n") )
		for line in INFILE:
			if len(line) > 1:
				# count non-empty lines
				line_cnt += 1
	npop = ( float(line_cnt-2)/n_loci_per_replicate ) -1
	if not npop.is_integer():
		print ("failed to determine npop from spinput.txt")
		exit()
	return int(npop)

def read_spinput(npop):
	
	nsam_locilengths = []
	
	with open("spinput.txt", "r") as INFILE:
		INFILE.readline() # discard first line (empty)
		n_loci_per_replicate = int( INFILE.readline().strip("\n") )
		
		loci_cnt = 0
		while loci_cnt < n_loci_per_replicate:
			loc_sampling_scheme = []
			for p in range(npop):				
				nsam_p = int( INFILE.readline().strip("\n") )
				loc_sampling_scheme.append( nsam_p )
			locus_length = int( INFILE.readline().strip("\n") )
			loc_sampling_scheme.append( locus_length )
			nsam_locilengths.append( loc_sampling_scheme )
			
			loci_cnt += 1
		
		n_replicates = int( INFILE.readline().strip("\n") )
	
	return n_replicates, nsam_locilengths
	
def ms2counts(ms_out, nsam_and_locuslength):


	# catch loci without segsites but return fixed counts!
	try:

		sites_cnts = []
		nsites = len( ms_out[0] )

		for s in range( nsites ):
			totalsam = 0
			s_counts = []
			for pop in range(len(nsam_and_locuslength)-1):
				nsam = nsam_and_locuslength[pop]
	#			print "nsam", nsam 

				pop_samples = ms_out[totalsam:totalsam+nsam]
		#		print pop_samples
				totalsam += nsam
			
				resh = [ [row[idx] for row in pop_samples]  for idx in range(len(pop_samples[0])) ] 
	#			print len(resh)
	
				cnt1 = resh[s].count("1")
				cnt0 = nsam - cnt1
				s_counts.append( cnt1 )
				s_counts.append( cnt0 )
			sites_cnts.append( s_counts )
	
	except IndexError:

		fixed = [ ]
		for pop in range(len(nsam_and_locuslength)-1):
			nsam = nsam_and_locuslength[pop]
			fixed.append(nsam)
			fixed.append(0)
		sites_cnts = [ fixed ]
		
#	print "newmethod", sites_cnts
	"""
	final structure of output is:
	site_1_pop_1_count_1 site_1_pop_1_count_0 site_1_pop_2_count_1 site_1_pop_2_count_0 site_1_pop_3_count_1 site_1_pop_3_count_0
	site_2_pop_1_count_1 site_2_pop_1_count_0 site_2_pop_2_count_1 site_2_pop_2_count_0 site_2_pop_3_count_1 site_2_pop_3_count_0
	site_3_pop_1_count_1 site_3_pop_1_count_0 site_3_pop_2_count_1 site_3_pop_2_count_0 site_3_pop_3_count_1 site_3_pop_3_count_0
	"""
	return sites_cnts

def count_pairwise_differences (a,b):
	
	diffs = 0
	for idx, x in enumerate(a):
		if x != b[idx]:
			diffs += 1
	
	return(float(diffs))

def make_haplotype_stats (samples):
	
	if len(samples) > 0:
		resdict = {h : samples.count(h) for h in set(samples)}
	else:
		samples = ["0"]
		resdict = {"0":1}
	n_samples = float(len(samples))
	n_haplotypes = float(len(resdict.keys()))

	freqs = sorted([ x/n_samples for x in resdict.values() ], key = float)
	outstats = []
	if n_haplotypes > 2:
		outstats += [n_haplotypes, freqs[-1], freqs[-2], freqs[0]]
	elif n_haplotypes == 2:
		outstats += [2.0, freqs[-1], freqs[-2], 0.0]
	else: # n_haplotypes <= 1:
		outstats += [1.0, 1.0, 0.0, 0.0]
					
	h_diversity = 1.0 - sum( [f**2 for f in freqs ] ) ## == expected individual heterozygosity if were diploid
	outstats += [h_diversity]
					
	# number of differences between haplotypes
	pwdict = {}
	if n_haplotypes > 1:
		for a in resdict.keys():
			for b in resdict.keys():
				if a != b:
					pair = ".".join(sorted([a,b]))
					if not pair in pwdict.keys():
						pwdict[pair] = count_pairwise_differences (a,b)
	else:
		pwdict[0] = 0
	
	outstats += [ max(pwdict.values()), sum( pwdict.values() )/len(pwdict.values()) ]
	
	# which haplotypes are the two most common ones?
	if n_haplotypes > 1:
		first = [k for k,v in resdict.items() if v == max(resdict.values())]
#		print first
#		print sorted(resdict.values(), key = int)[-2]
		second = [k for k,v in resdict.items() if v == sorted(resdict.values(), key = int)[-2] ]
#		print second
		outstats += [ pwdict[".".join(sorted([first[0],second[-1]]))] ]
	else:
		outstats += [0.0]
	
	# calculate LD => mean r2
	# instead of loop over all pairs of SNPs, we take a random sample of SNP pairs, up to 50:
	haplotype_length = len(samples[0])
	possible_pairs = set()
	for a in range(haplotype_length):
		for b in range(haplotype_length):
			if a != b:	
				pair = ".".join(sorted([str(a),str(b)]))
				possible_pairs.add(pair)
	try:
		selected_pairs = np.random.choice(list(possible_pairs), 50, replace = False)
	except ValueError:
		selected_pairs = list(possible_pairs)
#	print selected_pairs

	LD_list = []
	for p in selected_pairs:
		[a,b] = [int(x) for x in p.split(".")]
		#print pair
		a_allele = "0"
		b_allele = "0"
		pa = [x[a] for x in samples].count(a_allele)/n_samples
		pb = [x[b] for x in samples].count(b_allele)/n_samples
		if not pa <= 0.5: # ensure that these are minor allele frequencies
			pa = 1.0 - pa
			a_allele = "1"
		if not pb <= 0.5:
			pb = 1.0 - pb
			b_allele = "1"
		pab = len([x for x in samples if x[a] == a_allele and x[b] == b_allele])/n_samples
				
		r2 = ((pab - (pa*pb))**2)/float( pa*(1-pa)*pb*(1-pb) )
		# r2_max, conditional on where in allele frequency space we are:
		if pa >= pb: # section S6
			r2_max = ((1-pa)*pb)/(pa*(1-pb))
		else: # section S7
			r2_max = (pa*(1-pb))/((1-pa)*pb)
						
#		print pa, pb, pab, r2
		LD_list.append( r2/r2_max )
#		print r2/r2_max, pa, pb, pab, r2, r2_max
	
	outstats.append(np.mean(LD_list))

	return( outstats )

	
	

##############
# MAIN
##############


if __name__ == "__main__":

	header = ["n_haplotypes_mean", "freq_hist_bin1_mean", "freq_hist_bin2_mean", "freq_hist_last_bin_mean", "haplotype_diverstiy_mean", "mean_haplotype_distance_mean", "max_haplotype_distance_mean", "most_common_haplotypes_distance_mean", "LD_mean"]

	with open("ABChaplostat.txt", "w") as F:
 		F.write("\t".join(header) + "\n")

	
	npop = determine_npop_from_spinput ()

	n_replicates, nsam_locilengths = read_spinput(npop) 

	n_loci = len( nsam_locilengths )

	locilengths = [ x[-1] for x in nsam_locilengths ]

	locus_cnt = 0
	replicate_cnt = 0
	
	slinecnt = 0
	
	samples = []
	loci_counts_of_replicate = []
	
	preceding_line = "placeholderstring"
	
	outlines = []

			
	# discard first 3 lines:
	sys.stdin.readline()
	sys.stdin.readline()
	sys.stdin.readline()
# 	first_locus_cmd = sys.stdin.readline().split("	")
# 	# nsam_locuslength = [ int(first_locus_cmd[idx]) for idx in [5,6,4] ]
# 	nsam_locuslength = nsam_locilengths[ 0 ]
# 	locilengths.append( nsam_locuslength[2] )
	statscollector = []	
	for line in sys.stdin:
#		print line
		
		try:
			first_char = int( line[0] ) # this means it is a sample line (most lines are)
			samples.append( line.strip("\n") )
			slinecnt += 1
#			print "read chrom", slinecnt
    		
		except ValueError: # first char is not an int, so a character line
			
			if line.startswith("//"):
				# this means that a new locus starts
				locus_cnt += 1
				samples = []
				# nsam_locuslength = [ int( line.split("	")[idx] ) for idx in [5,6,4] ]
				# locilengths.append( nsam_locuslength[2] )	
				
									
			elif len(line) == 1:	
				# empty line : means the previous locus is finished OR that there were no samples in this locus because of segsites: 0
				if len(preceding_line) == 1:
					# skip, because locus was already finished with preceding empty line
					continue				
					
				else:
#					print "locus_finished, counting", len(samples)
			
					nsam_locuslength = nsam_locilengths[locus_cnt-1] # zero-based index!
#					site_cnts = ms2counts(samples, nsam_locuslength)
					
					haplotype_stats = make_haplotype_stats (samples)
					statscollector.append( haplotype_stats )
#					print statscollector
										
					# and send to stdout
					#sys.stdout.write( " ".join( [ " ".join( [str(cnt) for cnt in SNP] ) for SNP in site_cnts] ) + "\n")

					if (locus_cnt ) == n_loci :
						# this means that one replicate == set of x loci is full and we have to harvest it
# 						print "replicate_finished, harvesting stats"
# 						print locus_cnt, len(loci_counts_of_replicate)
						
						#print statscollector
						sumstats = [ np.nanmean([x[idx] for x in statscollector]) for idx in range(len(statscollector[0])) ]
						with open("ABChaplostat.txt", "a") as F:
 							F.write("\t".join([str(x) for x in sumstats]) + "\n")
#						print "statscollector reset"
						statscollector = []
						locus_cnt = 0

						
			else:
				continue
				
		preceding_line = line			
	
	# end of stdin: if stdin does not end with an empty line (which is the regular ms output), replicate harvest is not triggered! So do it here.
	if len(line) > 1:

		nsam_locuslength = nsam_locilengths[locus_cnt-1] # zero-based index!
		#print statscollector
		haplotype_stats = make_haplotype_stats (samples)
		statscollector.append( haplotype_stats )
		sumstats = [ np.nanmean([x[idx] for x in statscollector]) for idx in range(len(statscollector[0])) ]
		with open("ABChaplostat.txt", "a") as F:
			F.write("\t".join([str(x) for x in sumstats]) + "\n")
