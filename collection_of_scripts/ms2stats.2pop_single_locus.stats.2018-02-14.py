#!/usr/local/bin/python
# Python 2.7
# 14 Feb 2018
#
# ms2stats.2pop_single_locus.stats.2018-02-14.py
# 
# Mathias Scharmann
#
# FOR SINGLE LOCI! some stats are designed for samples/sets of loci, not single loci, e.g. standard deviations, but they
# are simply returned as zero
# switched off and returned as "na" here are:
# Fst_beta1, Fst_beta2, lower_limit, upperlimit = scipy.stats.beta.fit( Fst_list, floc=0, fscale=1) 
# piA_piB_Rsq = scipy.stats.pearsonr(piA_list, piB_list)[0]
#
# for single loci datasets, its possible that Tajima's D stats are empty when there there are no segsits, and return as numpy.float64("--").
# I let them return as an R-friendly "na" instead.


# INPUTS: 	- stdin in ms format
#			- "spinput.txt"-file: specifies number of loci, sample sizes and site-lengths for each locus

# OUTPUT: 	- "ABCstat.txt"-file with selected popgen summary statistics

# modular version: is c. 15% faster than the single script that counts AND claculates stats
# 
# example usage:
# or simply cat a file into it:
"""
cp Contig_removed_indels_names.fasta.genepop.txt.spinput.txt spinput.txt
cat Contig_removed_indels_names.fasta.genepop.txt.ms.txt | python $toolbase/ms2stats.counter.py | python ../tools_used/ms2stats.2pop_single_locus.stats.2018-02-14.py
mv ABCstat.txt Contig_removed_indels_names.fasta.ABCstat.txt
"""

##############
# HEAD
##############
from datetime import datetime

import sys
import numpy as np
import scipy.stats


##############
# FUNCTIONS
##############

def read_spinput():
	
	nsam_locilengths = []
	
	with open("spinput.txt", "r") as INFILE:
		INFILE.readline() # discard first line (empty)
		n_loci_per_replicate = int( INFILE.readline().strip("\n") )
		
		loci_cnt = 0
		while loci_cnt < n_loci_per_replicate:
			nsam1 = int( INFILE.readline().strip("\n") )
			nsam2 = int( INFILE.readline().strip("\n") )
			locus_length = int( INFILE.readline().strip("\n") )
			nsam_locilengths.append( [ nsam1, nsam2, locus_length ] )
			
			loci_cnt += 1
		
		n_replicates = int( INFILE.readline().strip("\n") )
	
	return n_replicates, nsam_locilengths

def get_freqs_from_count (all_loci_cnts):
	
	# has the structure:
	# [[pop1_freq1, pop2_freq1, all_pops_freq1],[pop1_freq1, pop2_freq1, all_pops_freq1]]
	
	# returns the global minor allele frequency (global MAF)
		
	all_loci_freqs = []
	for locus_idx in range(len(all_loci_cnts)):
		freqlist = []
		
		for s in all_loci_cnts[locus_idx]:
			try:
				freq1_1 = s[0] / ( s[0] + s[1])
			except ZeroDivisionError:
				freq1_1 = 0.0
			try:
				freq2_1 = s[2] / ( s[2] + s[3])
			except ZeroDivisionError:
				freq2_1 = 0.0
			try:
				freqall_1 = (s[0] + s[2]) / sum(s)
				if freqall_1 > 0.5:
					# here the minor allele is determined globally and the pop-wise AFs are correspondingly polarised
					# that means although AF is globally minor, a pop may have AF > 0.5!
					freqall_1 = 1.0 - freqall_1
					freq1_1 = 1.0 - freq1_1
					freq2_1 = 1.0 - freq2_1
				
			except ZeroDivisionError:
				freqall_1 = 0.0

			freqlist.append( [ freq1_1, freq2_1, freqall_1] )
		all_loci_freqs.append( freqlist )
	
	return all_loci_freqs

def count_segsites_stats (all_loci_freqs):
	
	# input has the structure:
	# [[pop1_freq1, pop2_freq1, all_pops_freq1],[pop1_freq1, pop2_freq1, all_pops_freq1]]
	
	ssglob_list = []		# segregating sites global
	ssA_list = []			# segregating sites A
	ssB_list = []			# segregating sites B
	
	paA_list = []			# private alleles A
	paB_list = []			# private alleles B
	
	srecf_list = []				# sites reciprocally fixed between A and B
		
	for locus in all_loci_freqs:
		ssglob = 0
		ssA = 0
		ssB = 0
		
		paA = 0
		paB = 0
		
		srecf = 0 # sites reciprocally fixed
				
		for site in locus: # check the global MAFs
			if site[0] > 0:
				ssA += 1
				if site[1] == 0:
					paA += 1
			if site[1] > 0:
				ssB += 1
				if site[0] == 0:
					paB += 1	
			if site[0] == 1 and site[1] == 0 or site[0] == 0 and site[1] == 1:
				srecf += 1
			if site[2] > 0: 
				ssglob += 1
			
		ssA_list.append(ssA)
		ssB_list.append(ssB)
		ssglob_list.append(ssglob)
		
		paA_list.append(paA)
		paB_list.append(paB)
		
		srecf_list.append(srecf)
		
	[ssA_avg, ssB_avg, ssglob_avg, paA_avg, paB_avg, srecf_avg] = [np.mean(x) for x in [ssA_list, ssB_list, ssglob_list, paA_list, paB_list, srecf_list]]
	[ssA_std, ssB_std, ssglob_std, paA_std, paB_std, srecf_std] = [np.std(x) for x in [ssA_list, ssB_list, ssglob_list, paA_list, paB_list, srecf_list]]
	
	segsite_stats = [ssA_avg, ssA_std, ssB_avg, ssB_std, ssglob_avg, ssglob_std, paA_avg, paA_std, paB_avg, paB_std, srecf_avg, srecf_std]
	
	return segsite_stats, ssA_list, ssB_list

def make_watterson_theta_and_tajimas_D (pi_per_site, nsites, segsites, samplesize):
	

	#	implemented as on wikipedia page, reference:
	#	Tajima (1989) Genetics 123 (3): 585-95
	# 	TajD is undefined (numpy.nan) if there are no segsites
	
	def n_minus_1_harmonic_number (x):
		harmsum = 0.0
		for k in range(1,int(x)): # range(1,x) expands to 1,2,...,x-1
			harmsum += 1.0/k
		return harmsum
	
	def n_minus_1_harmonic_number_squared (x):
		hsq = 0.0
		for k in range(1,int(x)): # range(1,x) expands to 1,2,...,x-1
			hsq += 1.0/(k*k)
		return hsq
	
	def calc_TajD_constant (segsites, samplesize):		
		a1 = n_minus_1_harmonic_number(samplesize)
		a2 = n_minus_1_harmonic_number_squared(samplesize)
		b1 = (samplesize + 1.0) / (3*(samplesize-1.0))
		b2 = (2*(samplesize*samplesize+samplesize+3.0))/ (9*samplesize*(samplesize-1.0))	
		c1 = b1 - 1.0/a1
		c2 = b2 - ((samplesize+2)/(a1*samplesize)) + (a2/(a1*a1))
		e1 = c1/a1
		e2 = c2/((a1*a1)+a2)	
		C = np.sqrt((e1*segsites)+(e2*segsites*(segsites-1)))	
		return C
	
	theta_watterson = segsites / n_minus_1_harmonic_number (samplesize)
	TajD_constant = calc_TajD_constant (segsites, samplesize)	
	tajd = ((pi_per_site*nsites) - theta_watterson) / TajD_constant
	
	return theta_watterson, tajd


def get_stats( replicates_counts, locilengths ):
	
	# input has the structure (example: a replicate with 2 loci of which the first had 3 seqgsites and the 2nd had 2 segsites):s
	# [ [ [cnt1_1, cnt1_0, cnt2_1, cnt2_0 ],[cnt1_1, cnt1_0, cnt2_1, cnt2_0 ],[cnt1_1, cnt1_0, cnt2_1, cnt2_0 ] ], [ [cnt1_1, cnt1_0, cnt2_1, cnt2_0 ],[cnt1_1, cnt1_0, cnt2_1, cnt2_0 ] ] ]
	
	all_loci_cnts = [ np.asarray(x, dtype = "float" ) for x in replicates_counts ]	
	
	all_loci_freqs = get_freqs_from_count (all_loci_cnts)
	
	segsite_stats, ssA_list, ssB_list = count_segsites_stats (all_loci_freqs)
	
	n_list = [ [ sum(y) for y in x ] for x in replicates_counts ]
	pi_list = get_pi ( all_loci_freqs, locilengths, 2, n_list )
	
	n_list = [ [ sum([y[0], y[1]]) for y in x ] for x in replicates_counts ]
	piA_list = get_pi ( all_loci_freqs, locilengths, 0, n_list )
	
	n_list = [ [ sum([y[2], y[3]]) for y in x ] for x in replicates_counts ]
	piB_list = get_pi ( all_loci_freqs, locilengths, 1, n_list )
		
	dxy_list = make_dxy (all_loci_freqs, locilengths)
	
	netDiv_list = dxy_list - ((piA_list + piB_list) / 2.0 ) 

	"""
	netDiv = dxy - mean([piA, piB])
	Nei & Li (1979): eq. 25 
	there: delta, also called Da, Dm etc.
	"""
	
	Fst_list = WC_Fst_per_locus (all_loci_freqs, all_loci_cnts)
		
	Fst_list[Fst_list <= 0.0] = 0.0001
	Fst_list[Fst_list >= 1.0] = 0.9999
	
	pi_mean = np.mean( pi_list )
	pi_std = np.std( pi_list )
	piA_mean = np.mean( piA_list )
	piA_std = np.std( piA_list )
	piB_mean = np.mean( piB_list )
	piB_std = np.std( piB_list )
	dxy_mean = np.mean( dxy_list )
	dxy_std = np.std( dxy_list )
	netDiv_mean = np.mean( netDiv_list )
	netDiv_std = np.std( netDiv_list )
	Fst_mean = np.mean( Fst_list )
	Fst_std = np.std( Fst_list )
	
	###
	## now watterson_theta, tajimas_D per population:
	# needs data from (by index): piA_list[locus_idx], locilengths[locus_idx], ssA_list[locus_idx], samplesize_A
	# samplesize_A: sum( all_loci_cnts[locus_idx][0][:2] )
	# samplesize_B: sum( all_loci_cnts[locus_idx][0][2:4] )
	
	tw_tjD_list = [ make_watterson_theta_and_tajimas_D ( piA_list[locus_idx], locilengths[locus_idx], ssA_list[locus_idx], sum( all_loci_cnts[locus_idx][0][:2] )  ) for locus_idx in range(len(locilengths)) ]
	thetaW_A_list = [x[0] for x in tw_tjD_list]
	TajD_A_list = [x[1] for x in tw_tjD_list]

	tw_tjD_list = [ make_watterson_theta_and_tajimas_D ( piB_list[locus_idx], locilengths[locus_idx], ssB_list[locus_idx], sum( all_loci_cnts[locus_idx][0][2:4] )  ) for locus_idx in range(len(locilengths)) ]
	thetaW_B_list = [x[0] for x in tw_tjD_list]
	TajD_B_list = [x[1] for x in tw_tjD_list]
	
	thetaW_A_avg = np.mean( thetaW_A_list )
	thetaW_A_std = np.std( thetaW_A_list )
	thetaW_B_avg = np.mean( thetaW_B_list )
	thetaW_B_std = np.std( thetaW_B_list )
	
	# TajD is undefined (numpy.nan) if there are no segsites, so we just drop numpy.nan :
	[TajD_A_avg, TajD_B_avg] = [ np.mean( np.ma.masked_array( x , np.isnan( x )) ) for x in [ TajD_A_list, TajD_B_list] ]
	[TajD_A_std, TajD_B_std] = [ np.std( np.ma.masked_array( x , np.isnan( x )) ) for x in [ TajD_A_list, TajD_B_list] ]
	
	# for single loci datasets, its possible that these stats are empty when there there are no segsits, and return as "--".
	# let them return as an R-friendly "na" instead:
	if "--" in str(TajD_A_avg): TajD_A_avg = "na"
	if "--" in str(TajD_B_avg): TajD_B_avg = "na"
	if "--" in str(TajD_A_std): TajD_A_std = "na"
	if "--" in str(TajD_B_std): TajD_B_std = "na"	
	####	
	
	
#	["ssA_avg", "ssA_std", "ssB_avg", "ssB_std", "ssglob_avg", "ssglob_std", "paA_avg", "paA_std", "paB_avg", "paB_std", "srecf_avg", "srecf_std", "thetaW_A_avg", "thetaW_A_std", "thetaW_B_avg", "thetaW_B_std", "TajD_A_avg", "TajD_B_avg", "TajD_A_std", "TajD_B_std"]
	
#		[ssA_avg, ssA_std, ssB_avg, ssB_std, ssglob_avg, ssglob_std, paA_avg, paA_std, paB_avg, paB_std, srecf_avg, srecf_std, thetaW_A_avg, thetaW_A_std, thetaW_B_avg, thetaW_B_std, TajD_A_avg, TajD_B_avg, TajD_A_std, TajD_B_std]
	
	
	
	# fit a beta distrib and force it to acknowledge that the input data is already in the interval 0,1
#	Fst_beta1, Fst_beta2, lower_limit, upperlimit = scipy.stats.beta.fit( Fst_list, floc=0, fscale=1) 	
	Fst_beta1, Fst_beta2, lower_limit, upperlimit = ["na","na","na","na"]	

	# percentiles:
	dxy_perc05, dxy_perc25, dxy_perc50, dxy_perc75, dxy_perc95 = np.percentile(dxy_list, [5, 25, 50, 75, 95])
	netDiv_perc05, netDiv_perc25, netDiv_perc50, netDiv_perc75, netDiv_perc95 = np.percentile(netDiv_list, [5, 25, 50, 75, 95])
	Fst_perc05, Fst_perc25, Fst_perc50, Fst_perc75, Fst_perc95 = np.percentile(Fst_list, [5, 25, 50, 75, 95])
	
	# correlaction R2 between the pi of A vs B:
#	piA_piB_Rsq = scipy.stats.pearsonr(piA_list, piB_list)[0]
	piA_piB_Rsq = "na"
	
	# local maxima in distrib: was not useful in classification according to random forest 
# 	dxy_peaks = find_peaks_in_hist (dxy_list)
#	Fst_peaks = find_peaks_in_hist (Fst_list)

	# bimodality?
	dxy_hist_lowest_bins, dxy_hist_subhighest_bins, dxy_hist_highest_bins = simple_bimodality_check (dxy_list, [np.min(dxy_list), np.max(dxy_list)] )
	netDiv_hist_lowest_bins, netDiv_hist_subhighest_bins, netDiv_hist_highest_bins = simple_bimodality_check (netDiv_list, [np.min(netDiv_list), np.max(netDiv_list)] )
	Fst_hist_lowest_bins, Fst_hist_subhighest_bins, Fst_hist_highest_bins = simple_bimodality_check (Fst_list, [0.0,1.0])
	
	## get 2D SFS / joint site frequency spectrum JSFS 2DSFS
	JSFS = make_JSFS (all_loci_freqs) # returns a list of the bin densities in a folded 2D SFS, length of 55 bins for a 2D hist of 10x10 fields
	
	# finish
#	stats = [ pi_mean, pi_std, piA_mean, piA_std, piB_mean, piB_std, dxy_mean, dxy_std, netDiv_mean, netDiv_std, Fst_mean, Fst_std, Fst_beta1, Fst_beta2, dxy_perc05, dxy_perc25, dxy_perc50, dxy_perc75, dxy_perc95, netDiv_perc05, netDiv_perc25, netDiv_perc50, netDiv_perc75, netDiv_perc95, Fst_perc05, Fst_perc25, Fst_perc50, Fst_perc75, Fst_perc95, piA_piB_Rsq, dxy_peaks, Fst_peaks, dxy_hist_lowest_bins, dxy_hist_subhighest_bins, dxy_hist_highest_bins, netDiv_hist_lowest_bins, netDiv_hist_subhighest_bins, netDiv_hist_highest_bins, Fst_hist_lowest_bins, Fst_hist_subhighest_bins, Fst_hist_highest_bins ] + JSFS
	
	stats = segsite_stats + [thetaW_A_avg, thetaW_A_std, thetaW_B_avg, thetaW_B_std, TajD_A_avg, TajD_A_std, TajD_B_avg, TajD_B_std] + [ pi_mean, pi_std, piA_mean, piA_std, piB_mean, piB_std, dxy_mean, dxy_std, netDiv_mean, netDiv_std, Fst_mean, Fst_std, Fst_beta1, Fst_beta2, dxy_perc05, dxy_perc25, dxy_perc50, dxy_perc75, dxy_perc95, netDiv_perc05, netDiv_perc25, netDiv_perc50, netDiv_perc75, netDiv_perc95, Fst_perc05, Fst_perc25, Fst_perc50, Fst_perc75, Fst_perc95, piA_piB_Rsq, dxy_hist_lowest_bins, dxy_hist_subhighest_bins, dxy_hist_highest_bins, netDiv_hist_lowest_bins, netDiv_hist_subhighest_bins, netDiv_hist_highest_bins, Fst_hist_lowest_bins, Fst_hist_subhighest_bins, Fst_hist_highest_bins ] + JSFS

# 	array_to_file ("pi_tot.txt", pi_list)
# 	array_to_file ("dxy.txt", dxy_list)
# 	array_to_file ("netDiv.txt", netDiv_list)
# 	array_to_file ("FstWC.txt", Fst_list)
	
	return stats

def array_to_file (outfname, outarray):
	
	outlines = "\n".join( [ str(x) for x in outarray ] )
	with open(outfname, "w") as OUTF:
		OUTF.write(outlines)
		
	
def simple_bimodality_check (mydat, myrange):
	
	"""
	returns mean bin frequencies for 3 groups of bins in a normalised histogram
	"""
	hist_freqs = np.histogram(mydat, bins = 30, range = myrange)[0]/float(len(mydat)) # normalised histogram

	a = np.mean( hist_freqs[:3] ) 		# the 3 lowest bins			: _hist_lowest_bins
	b = np.mean( hist_freqs[24:27] ) 	# the 3 bins left of highest: _hist_subhighest_bins
	c = np.mean( hist_freqs[27:] )		# the 3 highest bins		: _hist_highest_bins
	
# 	print hist_freqs[:3]
# 	print hist_freqs[24:27]
#  	print hist_freqs[27:]
#  	print hist_freqs

	return	a, b, c
	

def find_peaks_in_hist (mydat):
	
	"""
	- bins data into histogram with 30 bins
	- finds the number of peaks (local maxima) in the histogram
	- each peak needs a prominence of at least 5% of the maximum bin count in the histogram (to drop noise / insignificant peaks)
	"""
	
	my_hist = np.histogram( mydat, bins = 30)

	data = list( my_hist[0] )
	data.insert(0,0)
	data.insert(len(data),data[len(data)-1])
	x = my_hist[1]

	c = (np.diff( np.sign( np.diff(data) )) < 0).nonzero()[0] + 1 # indices of local maxima (always includes first bin)

	y_values_of_maxima = [data[idx] for idx in c]

	peaks_in_hist = len( [ x for x in y_values_of_maxima if x > 0.05 * max(y_values_of_maxima) ] )
	
	return peaks_in_hist
	

def make_dxy (all_loci_freqs, locilengths):
	
	"""
	Nei & Li (1979): unnumbered eq., between eqs. 24 and 25 
	Nei 1987 eq 10.20 
	for a single biallelic SNP (1,2) in two pops it simplifies to: dxy = pop1_freq1 * ( 1.0 - pop2_freq1 ) + pop2_freq1 * ( 1.0 - pop1_freq1 )
	"""
	
	loci_dxy_list = []
	for locus_idx in range(len(all_loci_freqs)):
		dxysum = 0.0
		for s in all_loci_freqs[locus_idx]:
			dxy = s[0] * ( 1.0 - s[1] ) + s[1] * ( 1.0 - s[0] )
			dxysum += dxy
		dxy = dxysum / locilengths[locus_idx]
		loci_dxy_list.append( dxy )
		
	return np.array( loci_dxy_list )
	
def get_pi(all_loci_freqs, locilengths, idx, n_list):
	
	"""
	for biallelic SNPs:
	nucleotide diversity pi = sum_over_all_SNPs (freq_p * (1-freq_p) * 2 * (n/(n-1.0)) ) / locus_length
	Nei & Li (1979)
	Gillespie (2004) Population Genetics: A concise guide. 2nd Ed., p. 45 top
	=> not exhaustively building all pairs of sampled alleles and counting segsites , but using the allele frequencies.
	(faster but identical result with correction factor n/(n-1) )
	"""

	loci_pi_list = []
	for locus_idx in range(len(all_loci_freqs )):
#		print locus_idx
		pisum = 0.0
		for site_idx in range( len( all_loci_freqs[locus_idx])):
			pi = all_loci_freqs[locus_idx][site_idx][idx] * ( 1 - all_loci_freqs[locus_idx][site_idx][idx] ) * 2 * (n_list[locus_idx][site_idx]/( n_list[locus_idx][site_idx]-1.0 ))
			pisum += pi
		pi = pisum / locilengths[locus_idx]
		loci_pi_list.append( pi )
		
	return np.array( loci_pi_list )

def WC_Fst_per_locus (all_loci_freqs, all_loci_cnts):
	
	"""
	- returns Weir & Cockerham's 1984 Fst estimator theta_hat averaged over the number of SNPs per locus
	- negative estimates are returned as 0.0, because that is what they mean.
	- if no SNPs in a locus, nan is returned
	- Heterozygosities are NOT taken from actual genotypes, but calculated under HWE!
	- that's because ms Hudson simulator returns chromosomes, not diploid genotypes
	"""
	
	loci_Fst_list = []
	for locus_idx in range(len(all_loci_freqs)):
		Fst_sum = 0.0
		for site_idx in range( len( all_loci_freqs[locus_idx])):
			p_A = all_loci_freqs[locus_idx][site_idx][0]
			p_B = all_loci_freqs[locus_idx][site_idx][1] 
			n_A = sum( all_loci_cnts[locus_idx][site_idx][:2] )
			n_B = sum( all_loci_cnts[locus_idx][site_idx][2:] ) 
			h_A = 1.0 - p_A**2 - (1.0-p_A)**2 # HWE for proportion of heterozygotes
			h_B = 1.0 - p_B**2 - (1.0-p_B)**2 # HWE for proportion of heterozygotes

			Fst = compute_WC84_Fst_per_site (p_A, p_B, n_A, n_B, h_A, h_B)
			Fst_sum += Fst
		Fst = Fst_sum / len( all_loci_freqs[locus_idx] )

		loci_Fst_list.append( Fst )
	
	# drop nans:
	loci_Fst_list = [ x for x in loci_Fst_list[:] if np.isfinite(x) ]

	# set negative Fsts to zero:
	loci_Fst_list = [ x if x >= 0.0 else 0.0 for x in loci_Fst_list[:] ]
		
	return np.array( loci_Fst_list )

def compute_WC84_Fst_per_site (p_i, p_j, n_i, n_j, h_i, h_j):	
	
	"""theta hat == Fst
	caution: first version was erroneous since not all divisions were floating-point! make sure everything is a float here!
	debugged using pythonised R-code from: Eva Chan 2008 calc_wcFstats.R
	"""
	
	# inputs to calculation:
	# allele freq of allele A in pop i
	p_i = float(p_i)
	p_j = float(p_j)
	
	# sample size in pop i
	n_i = float(n_i)
	n_j = float(n_j)
	
	r = 2.0 # number of populations; here of course 2
	
	# observed proportion of individuals heterozygous for allele A in pop i
	h_i = float(h_i)
	h_j = float(h_j)
	
	
	
	# n = average sample size
	n = (n_i + n_j) / r  # average sample size
	
	 # scv = squared coefficient of variation of sample sizes = square of (standard deviation / mean)
#	scv = ( ( math.sqrt( (1/r)*( (n_i-n)**2 + (n_j-n)**2 ) ) )/n )**2
#	scv = ( ( (n_i**2) + (n_j**2) ) - (n*n*r) ) / ( (n*n) * (r-1) ) # Eva Chan 2008 calc_wcFstats.R 	
#	n_c = n*(1-(scv/r)) 
	
	n_c = ((r*n) - ( ((n_i*n_i)/(r*n)) + ((n_j*n_j)/(r*n)) ) ) / (r - 1.0) # Eva Chan 2008 calc_wcFstats.R
	
	# p = average sample frequency of allele A
#	p = ( (n_i*p_i) + (n_j*p_j) ) / (r*n) # identical to Eva Chans result
	p = (n_i*p_i)/(r*n) + (n_j*p_j)/(r*n) # Eva Chan 2008 calc_wcFstats.R

	# s_2 = the sample variance of allele A frequnecies over populations
#	s_2 = ( (n_i*(p_i-p)**2) + (n_j*(p_j-p)**2) ) / ((r-1)*n)  # identical to Eva Chans result
	s_2 = (n_i*((p_i-p)**2)) / ((r-1)*n) + (n_j*((p_j-p)**2)) / ((r-1.0)*n) # Eva Chan 2008 calc_wcFstats.R
		
	# h = the average heterozygote frequency for allele A
#	h = ( (n_i*h_i) + (n_j*h_j) ) / (r*n) # identical to Eva Chans result
	h = (n_i*h_i)/(r*n) + (n_j*h_j)/(r*n) # Eva Chan 2008 calc_wcFstats.R
	
	# now the "three quantities"
#	a = (n / n_c )*(s_2 - (1 / (n-1) ) * ( p*(1-p) - ( (r -1)/r )*s_2 ) - ( (1/4)*h) ) # identical to Eva Chans result
	a = (n/n_c) * ( s_2 - ((1.0/(n-1.0))*((p*(1.0-p)) - (((r-1.0)/r)*s_2) - ((1.0/4.0)*h))) ) # Eva Chan 2008 calc_wcFstats.R
	
#	b = ( n/(n-1) ) * ( p*(1-p) - ((r-1)/r)*s_2 - ( (2*n - 1)/(4*n) )*h ) # identical to Eva Chans result
	b = (n/(n-1.0)) * ((p*(1.0-p)) - (((r-1.0)/r)*s_2) - ((((2.0*n)-1.0)/(4.0*n))*h)) # Eva Chan 2008 calc_wcFstats.R
	
	c = h/2.0 
	
	# finally, the Fst estimator theta; it is undefined for (a + b + c = 0) caused by fixation of the same allele in all samples
	if (a + b + c) == 0.0:
		theta_hat = float("-nan")
	else:
		theta_hat = a / (a + b + c)
	
	return theta_hat

def make_JSFS (all_loci_freqs):
	
	# all_loci_freqs has the structure:
	# [[[pop1_freq1, pop2_freq1, all_pops_freq1],[pop1_freq1, pop2_freq1, all_pops_freq1]],[[pop1_freq1, pop2_freq1, all_pops_freq1]]]
		
#	print "doing JSFS"
	# remove locus-structure from the SNP freq list:
	SNP_freqs = [item for sublist in all_loci_freqs for item in sublist]
	
	# remove invariant sites:
	SNP_freqs = [ x for x in SNP_freqs if x[2] > 0 ]
	
	MAFs_1 = [ x[0] for x in SNP_freqs ]
	MAFs_2 = [ x[1] for x in SNP_freqs ]
	
	JSFS = np.histogram2d(  MAFs_1 , MAFs_2, bins=[10,10], range = [[0.0, 1.0], [0.0, 1.0]])[0]
	totalbinsum = sum(sum(JSFS))
	JSFS_normed = JSFS / totalbinsum

#    	for row in JSFS_normed:
# 		print "\t".join( [ str(x) for x in row] )
	
	folded_JSFS = fold_2DSFS (JSFS_normed)
	
	# now flatten it:
	flat_folded_JSFS = [ str(round(item, 6)) for sublist in folded_JSFS for item in sublist]
#	print len(flat_folded_JSFS) 
		
	return flat_folded_JSFS

def fold_2DSFS (unfolded_2DSFS):

	"""
	from dadi code:
	Folded frequency spectrum
    
	The folded fs assumes that information on which allele is ancestral or
	derived is unavailable. Thus the fs is in terms of minor allele 
	frequency.  Note that this makes the fs into a "triangular" array.
    
	"""
	
	# reverse the unfolded spectrum along both axes and add it to the un-reversed spectrum:
	reversed_unfolded_2DSFS = np.flipud(np.fliplr(unfolded_2DSFS))
	not_yet_masked_2dSFS = unfolded_2DSFS + reversed_unfolded_2DSFS
#	print sum([sum(x) for x in not_yet_masked_2dSFS ]) # sum must be 2.0
	
# 	for row in not_yet_masked_2dSFS:
# 		print "\t".join( [ str(x) for x in row] )
	
#	print not_yet_masked_2dSFS

	bins_of_unfolded = len( unfolded_2DSFS[0] )
	bins_of_folded = bins_of_unfolded -1
    
    # now we construct a 2D list of only those entries (coordinates) that are allowed in the folded spectrum:
	folded_2DSFS = []
	# loop over the two dimensions:
	for i in range( bins_of_unfolded ):
		folded_2DSFS.append([])
		for j in range( bins_of_unfolded ) :	
			if (i + j ) < bins_of_folded : # this gices the diagonal of the array
#				print i,j
				folded_2DSFS[i].append( not_yet_masked_2dSFS[i][j] )
			elif (i + j ) == bins_of_folded :
				folded_2DSFS[i].append( unfolded_2DSFS[i][j] )
	
# 	for row in folded_2DSFS:
# 		print "\t".join( [ str(x) for x in row] )
	
#	print sum([sum(x) for x in folded_2DSFS]) # this sum must be 1.0 if the folded is correctly constructed
#	print sum([len(x) for x in folded_2DSFS]) # number of fields / entries in the folded

	return folded_2DSFS
	

##############
# MAIN
##############

if __name__ == "__main__":

	n_replicates, nsam_locilengths = read_spinput() 

	n_loci = len( nsam_locilengths )

	locilengths = [ x[2] for x in nsam_locilengths ]
	
	startTime = datetime.now()
	
	f = open("progress.log.txt", "w")
	f.write( "\t".join( [str(x) for x in [len( locilengths ), n_loci * n_replicates, "started:", startTime ]]) +"\n" )
	f.close()
	
	locus_cnt = 0
	replicate_cnt = 0
	
	slinecnt = 0
	
	samples = []
	loci_counts_of_replicate = []
	
	outlines = []
#	ABCstats_header = [ "pi_mean", "pi_std", "piA_mean", "piA_std", "piB_mean", "piB_std", "dxy_mean", "dxy_std", "netDiv_mean", "netDiv_std", "Fst_mean", "Fst_std", "Fst_beta1", "Fst_beta2", "dxy_perc05", "dxy_perc25", "dxy_perc50", "dxy_perc75", "dxy_perc95", "netDiv_perc05", "netDiv_perc25", "netDiv_perc50", "netDiv_perc75", "netDiv_perc95", "Fst_perc05", "Fst_perc25", "Fst_perc50", "Fst_perc75", "Fst_perc95",  "piA_piB_Rsq", "dxy_peaks", "Fst_peaks", "dxy_hist_lowest_bins", "dxy_hist_subhighest_bins", "dxy_hist_highest_bins", "netDiv_hist_lowest_bins", "netDiv_hist_subhighest_bins", "netDiv_hist_highest_bins", "Fst_hist_lowest_bins", "Fst_hist_subhighest_bins", "Fst_hist_highest_bins" ] + [ "2DSFS_bin_" + str(i) for i in range(55) ]
	
	ABCstats_header = ["ssA_avg", "ssA_std", "ssB_avg", "ssB_std", "ssglob_avg", "ssglob_std", "paA_avg", "paA_std", "paB_avg", "paB_std", "srecf_avg", "srecf_std", "thetaW_A_avg", "thetaW_A_std", "thetaW_B_avg", "thetaW_B_std", "TajD_A_avg", "TajD_A_std", "TajD_B_avg", "TajD_B_std"] + [ "pi_mean", "pi_std", "piA_mean", "piA_std", "piB_mean", "piB_std", "dxy_mean", "dxy_std", "netDiv_mean", "netDiv_std", "Fst_mean", "Fst_std", "Fst_beta1", "Fst_beta2", "dxy_perc05", "dxy_perc25", "dxy_perc50", "dxy_perc75", "dxy_perc95", "netDiv_perc05", "netDiv_perc25", "netDiv_perc50", "netDiv_perc75", "netDiv_perc95", "Fst_perc05", "Fst_perc25", "Fst_perc50", "Fst_perc75", "Fst_perc95",  "piA_piB_Rsq", "dxy_hist_lowest_bins", "dxy_hist_subhighest_bins", "dxy_hist_highest_bins", "netDiv_hist_lowest_bins", "netDiv_hist_subhighest_bins", "netDiv_hist_highest_bins", "Fst_hist_lowest_bins", "Fst_hist_subhighest_bins", "Fst_hist_highest_bins" ] + [ "2DSFS_bin_" + str(i) for i in range(55) ]
	
	with open("ABCstat.txt", "w") as OUTFILE:
		OUTFILE.write("\t".join(ABCstats_header) + "\n" )
				
		for line in sys.stdin:
#			print line
			
			try:
				first_char = int( line[0] ) # this means it is a sample line (most lines are)
				indat = [ int(x) for x in line.strip("\n").split(" ") ]
#				indat = [ x.split(",") for x in line.strip("\n").split("-")]
#				loci_counts_of_replicate.append( [[int(a) for a in b] for b in indat ] )
				loci_counts_of_replicate.append( [indat[x:x+4] for x in range(0, len(indat), 4) ] )
			
#				slinecnt += 1
#				print "read counts", slinecnt
    			
			except ValueError: # it must be an empty line , indicates end of replicate
				
				if replicate_cnt < n_replicates:
				
	#				print loci_counts_of_replicate
				
					# this means that one replicate == set of x loci is full and we have to harvest it
# 					print "replicate_finished, harvesting stats"
# 					print len(loci_counts_of_replicate)
						
					replicate_cnt += 1
					# sometimes the beta-fit fails, so we except it here and append a line full of NA to the output instead (a lost replicate)
					try:
						stats = get_stats (loci_counts_of_replicate, locilengths )	
					except scipy.stats.distributions.FitSolverError:
						stats = ["na"]*len(ABCstats_header)
					#print stats
					outlines.append( "\t".join( [ str(x) for x in stats ] ) )
					loci_counts_of_replicate = []
					if ( replicate_cnt / 100.0 ).is_integer():
						OUTFILE.write("\n".join( outlines )+"\n")
						outlines = []
 						f = open("progress.log.txt", "a")
 						f.write( "\t".join( [str(x) for x in ["replicate", replicate_cnt, "at", (datetime.now() - startTime) ]]) +"\n" )
 						f.close()
	
	# write one last time if not already written the results:
	if len(outlines) > 0:
		with open("ABCstat.txt", "a") as OUTFILE:
			OUTFILE.write("\n".join( outlines )+"\n")
		with open("progress.log.txt", "a") as f:
			f.write( "\t".join( [str(x) for x in ["replicate", replicate_cnt, "at", (datetime.now() - startTime) ]]) +"\n" )
		
	f = open("progress.log.txt", "a")
	f.write( "\t".join( [str(x) for x in ["Successful and done!", datetime.now() - startTime ]]) +"\n" )
	f.close()
			
