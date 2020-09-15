#!/usr/local/bin/python
# Python 2.7.6
# vcf_to_genepop.py
# 12 Feb 2018
# Mathias Scharmann


# usage example
# python vcf_to_genepop.py --vcf bla.vcf --popfile popfile.txt

"""
DOES NOT DEAL WITH HETEROZYGOTE AND AMBIGUOUS BASES!!
"""


import argparse
import os

########################## HEAD

# this function checks if file exists and exits script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format( x )
		exit()
	x = str(x)
	return x

# breaks script if non-UNIX linebreaks in input file popmap
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
		exit()

######

# parses arguments from command line
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--fasta", required=True, type=extant_file,
		help="name/path of fasta alignment input file", metavar="FILE")
	parser.add_argument("--popfile", required=True, type=extant_file,
		help="name/path of population file", metavar="FILE")


# 	parser.add_argument("--format", required=True, help="the format of the vcf: stacks or freebayes")
# 	
 	args = parser.parse_args()
		
	# finish
	return args

################################## CORE

def read_pop_info (INFILE):
	
	popdict = {}
	pop_order_in_outfile = []
	with open(INFILE, "r") as infile:
		for line in infile:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				if fields[1] in popdict.keys():
					popdict[fields[1]].append(fields[2])
				else: 
					popdict[fields[1]] = [fields[2]]
					pop_order_in_outfile.append( fields[1] )

	return popdict, pop_order_in_outfile


def parse_fasta(fasta_file, idv_seqs_dict ):
	
	outlines = []

	with open(fasta_file, "r") as INFILE:
	
		first_id = INFILE.readline().strip("\n")
		outlines.append(first_id.strip(">"))
		seq = ""
	
		for line in INFILE:
			line_clean = line.strip("\n")
			if line_clean.startswith(">"):
				outlines.append(seq)
				outlines.append(line_clean.strip(">"))
				seq = ""
			else:
				if len(line_clean) > 0:
					seq += line_clean
	
	# append last seq
	outlines.append(seq)
	
	fasta_dict = {}
	for i in xrange(0,len(outlines),2):
		fasta_dict[outlines[i]] = outlines[i+1]
	
	samples = idv_seqs_dict.keys()
	
	# initialise output collectors:
	genotypes_dict = {}
	for s in samples:
		genotypes_dict[s] = []
	
	loci_list = []
	snpcnt = 0
	for i in range(len(fasta_dict.values()[0])):
		totalnucs = list( set( [ x[i] for x in fasta_dict.values() ] ) )
		if len(totalnucs) > 2:
			continue # exclude variants that are not biallelic SNPs; we can't simulate anything else in ms anyway!
		snpcnt += 1
		locusname = "fasta_locus" + "-" + str(snpcnt)  # this name format is essential for conversion script to ms Hudson
		loci_list.append( locusname )

		for sample in samples:
			genocode = ""
			for chrom in [0,1]:
				base = fasta_dict[ idv_seqs_dict[sample][chrom] ][i]
				if base in ["A", "T", "G", "C"]:
					idx = totalnucs.index(base) + 1
					genocode += "0"+str(idx)
				else:
					genocode += "00"
#			print genocode
			genotypes_dict[sample].append( genocode )
				
	return genotypes_dict, loci_list


def make_genepop (genotypes_dict, loci_list, popdict, outfilename, pop_order_in_outfile):
	
	"""	
	headerline = meaningless
	lociline = all loci_IDs separated by ","
	pop
	indiv_ID,<tab>two_digits_per_allele<tab>two_digits_per_allele...
	indiv_ID,<tab>two_digits_per_allele<tab>two_digits_per_allele...
	.
	.
	.
	pop
	indiv_ID,<tab>two_digits_per_allele<tab>two_digits_per_allele...
	indiv_ID,<tab>two_digits_per_allele<tab>two_digits_per_allele...
	.
	.
	.
	"""
	headerline = "## genepop converted from fasta alignment for conversion to ms Hudson formats" + "	n_loci:	" + str( len(loci_list) ) + "	pops:	" + str(len( popdict.keys() )) + "	samples:	" + str( sum( [ len(x) for x in popdict.values() ] ) ) + "	pop_1 = " + pop_order_in_outfile[0] + "	pop_2 = " + pop_order_in_outfile[1]
	lociline = ",".join(loci_list)
	
	outlines = [headerline, lociline]
	for pop in pop_order_in_outfile:
		outlines.append("pop")
		for sample in popdict[pop]:
			try:
				outlines.append( sample + ",\t" + "\t".join( genotypes_dict[sample] ) )
			except KeyError:
				print "Warning: {0} on popmap but not found in vcf data".format(sample)
				continue
	
	with open(outfilename, "w") as OUTFILE:
		OUTFILE.write("\n".join( outlines ))
			
################################## MAIN


args = 	get_commandline_arguments ()
	
#print args
popdict, pop_order_in_outfile = read_pop_info(args.popfile)


idv_seqs_dict_prelim = {}
with open(args.popfile, "r") as INF:
	for line in INF:
		if len(line) > 1:
			fields = line.strip("\n").split("\t")
			try:
				idv_seqs_dict_prelim[fields[2]].append( fields[0] )
			except KeyError:
				idv_seqs_dict_prelim[fields[2]] = [fields[0]]

# diploidize
idv_seqs_dict = {}
for k,v in idv_seqs_dict_prelim.items():
	if len(v) < 2:
		idv_seqs_dict[k] = [v[0],v[0]]
	else:
		idv_seqs_dict[k] = v
			
		
#print popdict
genotypes_dict, loci_list = parse_fasta(args.fasta, idv_seqs_dict)
print "sites:	", len(loci_list)

outfilename = args.fasta + ".genepop.txt"

make_genepop (genotypes_dict, loci_list, popdict, outfilename, pop_order_in_outfile)
	
print "Done!"
	

