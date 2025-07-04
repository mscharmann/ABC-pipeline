#!/usr/local/bin/python
# Python 2.7
# genepop_to_mscalc.arbitrary_npop.py
# 2020
# Mathias Scharmann

"""

THIS version does an important sanity-check: each contig must have at least 2 samples; otherwise contig is discarded!

# any contig without VCF variants but present in the 'contig_lengths' inopu gets added here with FULL/maximum possible sample size!

- converts a genepopfile for X populations into formats readible by mscalc; i.e. the spinput.txt and a genotypes file in the format of Hudson's ms coalescent simulator.
- names of SNPs in the genepopfile need to strictly follow the convention: locusID_SNPposition
- this script DOES catch uneven missingness within RADtags, which is a problem for ms Hudson and mscalc! samples with unusual missingness are dropped from that RADtag

- Allows calculation of summary statistics from observed datasets with exactly the same method as the simulation summary statistics.


outputs:
	x.ms.txt		genotypes in MS Hudson format
	x.spinput.txt	ready to use for calculation of observed data statistics;
					for msnsam simulation stats calculated by mscalc using this file modify:
					 - pre-ultimate line to nreps, e.g. 10000
					 - ultimate line to name of file containing msnsam output, e.g. myfifo
	x.bpfile.txt	ready to use in msnsam simulation


"""
import sys
import os
import argparse
import random

# diese Funktion gibt command line arguments mit flag an das script weiter und checkt ob die Angaben plausibel / vollstandig sind
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("-contig_lengths", help="lengths of each locus, separate for rho (==total sites) and for mu", type=extant_file, required=True, metavar="FILE")
	parser.add_argument("-gp", help="genepop infile", type=extant_file, required=True, metavar="FILE")
	parser.add_argument("-mu", help="per site per generation mutation rate mu", default="7.5e-9")
	parser.add_argument("-rho", help="per site per generation recombination rate rho (probability of rec at a site)", default="7.5e-9")
	# default mu value is from Krasovec et al 2018 "The Mutation Rate and the Age of the Sex Chromosomes in Silene latifolia" Current Biology.
	parser.add_argument("-haplopops", help="comma-separated list of pops from which to sample only 1 chromosome per locus")
	
	args = parser.parse_args()

	# finish
	return args

#######

# diese Funktion checkt ob ein file existiert und stoppt das script wenn ein file nicht exisitiert
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print ("Error: {0} does not exist".format(x))
		exit()
	x = str(x)
	return x

###########

class Vividict(dict):
	def __missing__(self, key):
		value = self[key] = type(self)()
		return value
    
def de_vividict_loop (outer_dict1):
	
	out_dict = {}
	for key1, inner_dict1 in outer_dict1.items():
		if isinstance  (inner_dict1, Vividict):
			out_dict[key1] = dict(inner_dict1)
			outer_dict2 = inner_dict1
			for key2, inner_dict2 in outer_dict2.items():
				if isinstance  (inner_dict2, Vividict):
					out_dict[key1][key2] = dict(inner_dict2)
					outer_dict3 = inner_dict2
					for key3, inner_dict3 in outer_dict3.items():
						if isinstance  (inner_dict3, Vividict):
							out_dict[key1][key2][key3] = dict(inner_dict3)
	return out_dict	

#################

def get_nested_keys(dct):
	"""
	Safely get keys from nested dictionaries at two levels deep.
	Returns empty list if structure is missing or empty.
	"""
	first_level = next(iter(dct.values()), None)
	if first_level is None:
		return []
	second_level = next(iter(first_level.values()), None)
	if second_level is None:
		return []
	return list(second_level.keys())


def read_genepop (genepopfile):
	
	print ("reading genepop to dictionary")

	with open(genepopfile, "r") as INFILE:
		header = INFILE.readline() # remove first line
		hps = header.split(", ")
		pops_named_ordered = []
		for i in hps:
			if "pop_" in i:
				p = i.split("=")[1].strip()
				pops_named_ordered.append(p)
		print (pops_named_ordered)
		
		loci_ordered = INFILE.readline().strip("\n").split(",")
		uberdict = {}
		popcount = 0
		for line in INFILE:
			items = line.strip("\n").split("\t")
			if "pop" in items:
				population = pops_named_ordered[ popcount ]
				popcount += 1
				continue
			sample = "".join(items[:1]).strip(",")
			try:
				uberdict[population][sample] = items[1:]
			except KeyError:
				uberdict[population] = {sample:items[1:]}
# 				except KeyError:
# 					uberdict[population] = {sample:items[1:]}
		
	genotypes = {}
	for pop in uberdict.keys():
		genotypes[pop] = {}
		for sample in uberdict[pop].keys():
			print (sample)
			genotypes[pop][sample] = {}
			for idx in range(len(loci_ordered)): # previously used the zip() function here, but it is MUCH slower!
				locus = loci_ordered[idx]
				genotype = uberdict[pop][sample][idx]
#				print genotype
				genotypes[pop][sample][locus] = genotype
	
		print ("read genepop to dict, now filtering")

	
	# if within an individual, the sites of a RAD-tag have unequal presence or absence, then mscalc will crash:
	# "error in reading dataset 0 : cannot read datasetfile prout"
	# this can happen when e.g. the RAD-tag is longer than the read-length and was not mapped for that individual in some positions!
	# so samples with this type of problem need be excluded altogether from that RAD -> replace with missing data the offending genotypes!
	# the error will however only be thrown if this concerns SNPs! fixed sites can have unequal missingness
	
	replace_count = 0
	for pop in genotypes.keys():
		for sample in genotypes[pop].keys():
		
			locus_genotypes = [] # initiate
			locus_sites = [] # initiate
			previous_locus = loci_ordered[0].split("-")[0] # initiate
			for i in loci_ordered:
				this_locus = i.split("-")[0]
				
				if this_locus == previous_locus:
					locus_genotypes.append( genotypes[pop][sample][i] )
					locus_sites.append( i.split("-")[1]  ) # record sites of the locus
#					print locus_sites
				else:
# 					print sample
#   					print previous_locus
#   					print locus_genotypes

 					# evaluate collected genotypes from the PREVIOUS locus:
					pres_cnt = sum([ 1 for x in locus_genotypes if x != "0000"]) 
					abs_cnt = locus_genotypes.count("0000")
#  					print pres_cnt
#  					print abs_cnt
	
					if pres_cnt < len(locus_sites):
						if abs_cnt < len(locus_sites):
							# this means that not all SNPs in a RAD are together present or absent					
							for site in locus_sites:
								genotypes[pop][sample]["{0}-{1}".format(previous_locus, site)] = "0000" 	
								replace_count += 1			
 					# reset collectors: but not to empty but to the values from the current loop! otherwise, loci/sites will be skipped 
 					# => the first site in each locus will be skipped, which creates problems for the replacing with missingness and hence counting later on, 
 					# since RADtag presence counter checks only the first site per RADtag...
					locus_genotypes = [ genotypes[pop][sample][i] ] # reset collector
					locus_sites = [ i.split("-")[1] ]				
 				
				# update in each round
				previous_locus = this_locus
 
	print ("warning: {0} genotypes for contigs where sites within them had inconsistent number of missing data where replaced with missing data, because downstream cannot deal with this!".format(replace_count ) )
		
	return genotypes, loci_ordered, pops_named_ordered

def genotype_stats(genotypes, loci_ordered, pseudo_haploid_pops):
	
	print ("counting pres/abs of data")
	
	diplopops = [x for x in genotypes.keys() if not x in pseudo_haploid_pops]
	
	print("regular diploid pops: ", diplopops)
	print("pseudo-haploid pops: ", pseudo_haploid_pops)
	
	# this is for the spinput.txt
	loci_unique = {} # separate loci from their SNPs
	for locus in loci_ordered:
		uloc = locus.split("-")[0]
		snp = locus.split("-")[1]
		try:
			loci_unique[uloc].append(snp)
		except KeyError:
			loci_unique[uloc] = [snp]
	
	print ("first stage complete")
	
	data_pres_count = {} # this works because all SNPs/sites in one locus have same number of available genotypes (or at least I thought so and will make it happen -- see filtering function above.)
	
#	sample_256-rafflesiana-Brunei-18
#	print genotypes["pop_2"]["sample_256-rafflesiana-Brunei-18"]["22304_L114-20"]
	
	for pop in diplopops:
		data_pres_count[pop] = {}
		
		for locus, sites in loci_unique.items():
			data_pres_count[pop][locus] = 0
			
			for sample in genotypes[pop].keys():
				try:
					if genotypes[pop][sample]["{0}-{1}".format(locus, sites[0])] != "0000": # 0000 means missing data in this genepop format
#					print genotypes[pop][sample][locus]
						data_pres_count[pop][locus] += 2 # count 2 for each genotype because of 2 observed chromosomes
				except KeyError:
#					print "genotype missing"
					genotypes[pop][sample]["{0}-{1}".format(locus, sites[0])] = "0000"


	for pop in pseudo_haploid_pops:
		data_pres_count[pop] = {}
		
		for locus, sites in loci_unique.items():
			data_pres_count[pop][locus] = 0
			
			for sample in genotypes[pop].keys():
				try:
					if genotypes[pop][sample]["{0}-{1}".format(locus, sites[0])] != "0000": # 0000 means missing data in this genepop format
#					print genotypes[pop][sample][locus]
						data_pres_count[pop][locus] += 1 # count 1 for each genotype because of 1 observed chromosomes
				except KeyError:
#					print "genotype missing"
					genotypes[pop][sample]["{0}-{1}".format(locus, sites[0])] = "0000"

	
	print ("counted pres/abs")
					
	# bring data into a structure accessible for outputting ms-format:
	# it needs 5 dimensions...
	
	ms_datadict = {}
	for locus in loci_unique.keys():	# first level: locus
		ms_datadict[locus] = {}
		for pop in pseudo_haploid_pops:	# 2 level: populations
			ms_datadict[locus][pop] = {}
			for sample in genotypes[pop].keys():	# 3 level: samples
				ms_datadict[locus][pop][sample] = {}
				for snp in loci_unique[locus]:	# 4 level: the snps
					ms_datadict[locus][pop][sample][snp] = []	
					chroms = [genotypes[pop][sample]["{0}-{1}".format(locus, snp)][:2], genotypes[pop][sample]["{0}-{1}".format(locus, snp)][2:]]	# 5 level: the chromosomes, but omit missing genotypes "0000"
					if chroms != ["00","00"]:
						ms_datadict[locus][pop][sample][snp].append( random.choice(chroms) )
		for pop in diplopops:	# 2 level: populations
			ms_datadict[locus][pop] = {}
			for sample in genotypes[pop].keys():	# 3 level: samples
				ms_datadict[locus][pop][sample] = {}
				for snp in loci_unique[locus]:	# 4 level: the snps
					ms_datadict[locus][pop][sample][snp] = []	
					for chrom in [genotypes[pop][sample]["{0}-{1}".format(locus, snp)][:2], genotypes[pop][sample]["{0}-{1}".format(locus, snp)][2:]]:	# 5 level: the chromosomes, but omit missing genotypes "0000"
						if chrom != "00":
							ms_datadict[locus][pop][sample][snp].append(chrom)
	
	
	
	print ("made ms_datadict with 5 dimensions")
							
	out_msdict = Vividict()	
	# now re-assign allele identifiers into 0 or 1:
	# first, make set() of alleles per snp per locus, take into account all samples and populations:
	for locus in loci_unique.keys():
		first_level = next(iter(ms_datadict[locus].values()))
		second_level = next(iter(first_level.values()))
		for snp in second_level.keys():
			allels = set()
			for pop in genotypes.keys():
				for sample in genotypes[pop].keys():
					for chrom in ms_datadict[locus][pop][sample][snp]:
					 allels.add(chrom)
			
			# now, replace these with 1 or 0, or omit SNP if invariant in all samples and pops
			if len(allels) == 2:
				allels_dict = {}
				allels_dict[list(allels)[0]] = "0"
				allels_dict[list(allels)[1]] = "1"
			
				for pop in diplopops:
					for sample in genotypes[pop].keys():
						for chrom in [0,1]:
							# populate !
							if len(ms_datadict[locus][pop][sample][snp]) > 0: # last time getting rid of missing data/empty chrom lists
								out_msdict[locus][pop][sample][snp][chrom] = allels_dict[ ms_datadict[locus][pop][sample][snp][chrom] ] # looks up the allele and replaces allelename with ms-code
				for pop in pseudo_haploid_pops:
					for sample in genotypes[pop].keys():
						for chrom in [0]:
							# populate !
							if len(ms_datadict[locus][pop][sample][snp]) > 0: # last time getting rid of missing data/empty chrom lists
								out_msdict[locus][pop][sample][snp][chrom] = allels_dict[ ms_datadict[locus][pop][sample][snp][chrom] ] # looks up the allele and replaces allelename with ms-code


	return data_pres_count, loci_unique, out_msdict
	
def make_spinput(spinput_outfile, ms_outfile, gstats, loci_ordered, pop_order, contig_len_dict):
	
	print ("making spinput.txt")
	
	lines = []
	for locus in loci_ordered: # keep loci in order
		for pop in pop_order:
			try: ## loci without VCF variants get added here with FULL/maximum possible sample size!
				lines.append(gstats[pop][locus]) # the counts of number of samples per population, multiplied by 2 for number of chromosomes!!
			except KeyError:
				lines.append( max(gstats[pop].values()) )	
		locus_length = str( contig_len_dict[locus][1] ) # ### number of mutateable sites per contig
		lines.append( locus_length )
		
	printlines = "\n".join([str(x) for x in lines])
	
	with open(spinput_outfile, "w") as OUTFILE: 
		OUTFILE.write("\n")
		OUTFILE.write(str(len( loci_ordered )) + "\n") # the number of loci
		OUTFILE.write(printlines + "\n")
		OUTFILE.write("1\n") # the nreps; 1 since this is only 1 observation
		OUTFILE.write(ms_outfile + "\n" + "\n") # name of the file with data in ms format
		OUTFILE.close()
	
def make_ms(ms_outfile, out_msdict, loci_ordered, pop_order, pseudo_haploid_pops):

	print ("making ms output")
	
	lines = []
	for locus in loci_ordered:  # keep loci in order
		lines.append("// observed locus {0}".format(locus))
		
		# some loci may not contain any segregating sites either because they were not polymorphic or because the sites were lost in filtering.
		keys = get_nested_keys(out_msdict[locus])
		lines.append("segsites: {0}".format(len(keys)))
		lines.append("positions: {0}".format(" ".join(str(x) for x in keys) if keys else "-"))
		
		# the core: 0 or 1 for alleles
		for pop in pop_order:
			if pop in out_msdict[locus].keys():
				for sample in out_msdict[locus][pop].keys():
					if pop not in pseudo_haploid_pops:
						for chrom in [0,1]:
							chromlinelist = []
		#					print out_msdict[locus][pop][sample].keys()
							for snp in out_msdict[locus][pop][sample].keys():
								chromlinelist.append(out_msdict[locus][pop][sample][snp][chrom])
							chromline = "".join([str(x) for x in chromlinelist])				
							lines.append(chromline)
					else: # is a pseudo-haploid pop
						for chrom in [0]:
							chromlinelist = []
		#					print out_msdict[locus][pop][sample].keys()
							for snp in out_msdict[locus][pop][sample].keys():
								chromlinelist.append(out_msdict[locus][pop][sample][snp][chrom])
							chromline = "".join([str(x) for x in chromlinelist])				
							lines.append(chromline)
						
				
		lines.append("")
		
	printlines = "\n".join([str(x) for x in lines])
	
	with open(ms_outfile, "w") as OUTFILE: 
		OUTFILE.write("this is observed data\n")
		OUTFILE.write("no random seeds available\n")
		OUTFILE.write("\n")
		OUTFILE.write(printlines + "\n")

	
def make_bpfile (pop_order, gstats, bpoutfile, nreps, loci_ordered, mu, rho, contig_len_dict ):
	
	print ("making bpfile")
	nloci = len(loci_ordered)
	
	# give advice on msnam argument 2 
	msnsamarg = int(nloci) * int(nreps)
	print ("msnsam argument 2 must be:\t{0}".format(msnsamarg))
	
	# make the output
	
	outlines = []
	
#	line1 = "#Sp1={0} Sp2={1}\n".format(sp1, sp2)
	line1 = "# total pops: " + str(len(pop_order)) + " , ordered exactly as in genepop input\n" 	
	line2 = "\t".join( [str(contig_len_dict[l][1]) for l in loci_ordered] ) + "\n"  ### number of mutateable sites per contig
#	print (loci_ordered)
#	line2 = "\t".join( [ x.split("L")[1] for x in loci_ordered] ) + "\n"
	

	outlines.append(line1)
	outlines.append(line2)
	
	for pop in pop_order:	
		samples_pop = []
		for loc in loci_ordered:  # keep loci in order
			try: ## loci without VCF variants get added here with FULL/maximum possible sample size!
				samples_pop.append(str(gstats[pop][loc]))
			except KeyError:
				samples_pop.append(str( max(gstats[pop].values()) ) )
		outline = "\t".join(samples_pop) + "\n"
		outlines.append( outline )
	
	# per-population per-locus mutation rate theta is scaled to number of mutatebale sites per contig (may or may not be the same as total length)
	theta_values = [] # 
	for loc_len in [ contig_len_dict[l][1] for l in loci_ordered]:  # keep loci in order  # keep loci in order
		theta_values.append( str( loc_len * mu * 4 * 100000 ) )	# here mutation rate mu is multiplied by 4 Ne_ref and the number of sites per contig	
	outlines.append( "\t".join(theta_values) + "\n" )
	
	# intra-locus recombination rate rho is scaled to FULL length of the contig
	rho_values = []
	for loc_len in [ contig_len_dict[l][0] for l in loci_ordered]:
		rho_values.append( str( loc_len * rho * 4 * 100000 ) )
	outlines.append( "\t".join(rho_values) + "\n" )
	
	with open(bpoutfile, "w") as OUTFILE: 
		OUTFILE.write("".join( outlines ))

def read_contig_lengths (INFILE):
	
	contig_lengths_dict = {}
	with open(INFILE, "r") as INF:
		INF.readline() # discard header
		for line in INF:
			if len(line) > 2:
				fields = line.strip("\n").split()
				contig_lengths_dict[fields[0]] = [ int( fields[1] ), int( fields[2] ) ]
	print ("read ", INFILE)
	return contig_lengths_dict

############ MAIN

def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]
		
	args = 	get_commandline_arguments ()
	
	# for simulations, define here the spontaneous mutation rate mu
	# per site per generation per chromosome (NOT the population mutation rate but 'inidividual')
	# e.g. for Brassicaceae: 6.51548E-09 from De La Torre et al. (2017)
	mu = float(args.mu)
	
	# for simulations, define here the spontaneous intra-locus recombination rate rho
	# per site per generation per chromosome (NOT the population recombination rate but 'inidividual')
	rho = float(args.rho)
	
	# The mutation rate is scaled to a different number of sites than the rho or total contig lengths; this is useful when only 4-fold sites are desired. (mu_per_contig = mu * length_4-FDGS).
	# BUT the recombination rate should reflect that these sites were actually part of a contig that is much longer than the number of 4-FDGS
	# and hence more likely to recombine than the culled set of sites would suggest
	contig_len_dict = read_contig_lengths ( args.contig_lengths )
	
#	print(contig_len_dict)
	
	genepopfile = args.gp
	
	try:
		pseudo_haploid_pops = args.haplopops.split(",")
	except AttributeError:
		pseudo_haploid_pops = []
	if pseudo_haploid_pops == ["None"]:
		pseudo_haploid_pops = []
	
	
	nreps = "1" ## dummy value to be written to bpfile and spinput; suitable to calculate the obs data but can be changed later for simulations. 
		
	spinput_outfile = "{0}.spinput.txt".format(genepopfile)
	ms_outfile = "{0}.ms.txt".format(genepopfile)
	bpoutfile = "{0}.bpfile.txt".format(genepopfile)
	
	(genotypes, loci_ordered, pop_order) = read_genepop(genepopfile)
	# pop_x is the order of appearance in the genepopfile == consistent with genepop header!
 	
 	#print (loci_ordered)
 	
	(gstats, loci_unique, out_msdict) = genotype_stats(genotypes, loci_ordered, pseudo_haploid_pops)
	# pop_x is still the order of appearance in the genepopfile
	
	print(gstats)
	
	#loci_ordered = sorted( set( loci_unique.keys() + contig_len_dict.keys()) )
	
	# identify loci that have zero samples left in one or more populations:
	bad_loci = []
	for pop,stats in gstats.items():
		bs = [k for k,v in stats.items() if v < 2] 
		bad_loci += bs

	print("WARNING: " + str( len(set(bad_loci))) + " contigs had no samples left after cleaning for incomplete missingness along the contig; will be dropped")
#	print("HERE WE ARE")
	
	# we get rid of the bad contigs and sort the remainder
	loci_ordered = sorted( list( set(contig_len_dict.keys()) - set(bad_loci) ) )
	
	print ( "\n".join(bad_loci) )
	
	print ("effective number of contigs from contig_lengths file but then dropped here:	" + str( len( [x for x in bad_loci if x in set(contig_len_dict.keys())] ) ) )
	
	print (pop_order)
	
	# all these functions produce outputs where the order of populations is the same as in the genepop input
	make_bpfile (pop_order, gstats, bpoutfile, nreps, loci_ordered, mu, rho, contig_len_dict)
 	
	make_spinput(spinput_outfile, ms_outfile, gstats, loci_ordered, pop_order, contig_len_dict)
 	
	make_ms(ms_outfile, out_msdict, loci_ordered, pop_order, pseudo_haploid_pops) 
		
	print ("Done!")
	
	
if __name__ == "__main__":
	main(sys.argv[1:])
