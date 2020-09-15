## this script automatically decides which simulate.sh to use
# argfile => manually edited by user!

## FUNCTIONAL for 1, 2 and 3 pop models.

configfile: 'config.yaml'

import shutil, os

# calculate number of chunks and set up a list for the number of chunks:
n_chunks = int( int(config["total_sim_per_scenario"]) / int(config["n_sim_per_chunk"]) )

chunks = range(1,n_chunks+1)

print(config["scenarios"])

def choose_sim_script (popmapfile):
	
	# simulation script must be chosen according to the number and topology of pops desired.
	
	# look at popmap and count npop
	popv = []
	with open(popmapfile, "r") as I:
		for line in I:
			fields = line.strip("\n").split()
			popv.append(fields[1])
	
	npop = len(set(popv))
	
	# sanity-check (superficially) that the argfiles/models are compatible with the popmap (npop 1, 2, or 3).
	# fill a (global) dictionary to connect scenarions with their appropriate simulation.sh script
	scenario_to_simscript_dict = {}
	
	if npop == 1:
		for k,v in config["scenarios"].items():
			if not "argfiles/1pop" in v:
				print("Error! ", v, " can not be used for 1-population models. Popmap incompatible with chosen argfile (model).")
				exit()
			scenario_to_simscript_dict[k] = "simulate.one_pop.2020-09-14.sh"
	
	elif npop == 2:
		for k,v in config["scenarios"].items():
			if not "argfiles/2pop" in v:
				print("Error! ", v, " can not be used for 2-population models. Popmap incompatible with chosen argfile (model).")
				exit()
			scenario_to_simscript_dict[k] = "simulate.two_pops.2020-09-14.sh"
	
	elif npop == 3:
		for k,v in config["scenarios"].items():
			if not "argfiles/3pop" in v:
				print("Error! ", v, " can not be used for 3-population models. Popmap incompatible with chosen argfile (model).")
				exit()
			topo = v.split(".")[1]
			if topo == "1-3--2":
				scenario_to_simscript_dict[k] = "simulate.1-3--2.2020-09-15.sh"
			elif topo == "2-1--3":
				scenario_to_simscript_dict[k] = "simulate.2-1--3.2020-09-15.sh"
			elif topo == "2-3--1":
				scenario_to_simscript_dict[k] = "simulate.2-3--1.2020-09-15.sh"
			else:
				print("Error! ", v, " not a valid topology, re-name and change the argfile")
				exit()	
	return scenario_to_simscript_dict


scenario_to_simscript_dict = choose_sim_script (config["popmap"])
print(scenario_to_simscript_dict)	


rule all: ## the all rule MUST appear at the top !!
	input:
		expand("results/ABCstat.{scenario}.txt", scenario=config["scenarios"]),
		expand("results/parameters.{scenario}.txt", scenario=config["scenarios"]),
		expand("{vcf}.ABCstat.txt", vcf=config["vcf"])
#	run: # cleanup
#		shutil.rmtree('tempdir/')



rule collect_chunks:
	input:
		expand("tempdir/temp_{scenario}_{chunk}/ABCstat.txt", scenario=config["scenarios"], chunk=chunks),
		expand("tempdir/temp_{scenario}_{chunk}/parameters.txt", scenario=config["scenarios"], chunk=chunks)
	output:
		expand("results/ABCstat.{scenario}.txt", scenario=config["scenarios"]),
		expand("results/parameters.{scenario}.txt", scenario=config["scenarios"])
	run:
		for my_scenario in config["scenarios"]:
			print(my_scenario)
			with open("results/ABCstat." + my_scenario + ".txt", "w") as O:
				with open("tempdir/temp_" + my_scenario + "_1/ABCstat.txt", "r") as I:
					O.write(I.read())
				for c in ["tempdir/temp_" + my_scenario + "_" + str(x) + "/ABCstat.txt" for x in range(2,n_chunks + 1)]:
					print(c)
					with open(c, "r") as I:
						I.readline()
						O.write(I.read())
			with open("results/parameters." + my_scenario + ".txt", "w") as O:
				with open("tempdir/temp_" + my_scenario + "_1/parameters.txt", "r") as I:
					O.write(I.read())
				for c in ["tempdir/temp_" + my_scenario + "_" + str(x) + "/parameters.txt" for x in range(2,n_chunks + 1)]:
					print(c)
					with open(c, "r") as I:
						I.readline()
						O.write(I.read())
				

rule convert_inputs:
	output:
		expand("{vcf}.ABCstat.txt", vcf=config["vcf"])
	params:
		vcf = config["vcf"],
		popmap = config["popmap"],
		contig_lengths_file = config["contig_lengths_file"],
		mu = config["mu"],
		rho = config["rho"]
	conda:
		"abc_sims_numpy_scipy.yml"
	threads: 2
	shell:
		"""
		bash ./scripts/convert_inputs.sh {params.vcf} {params.popmap} {params.contig_lengths_file} {params.mu} {params.rho}
		"""


rule make_ABCstat_chunks:
	input:
		expand("{vcf}.ABCstat.txt", vcf=config["vcf"]),
		lambda wildcards: config["scenarios"][wildcards.scenario] # here we get which argfile to use per scenario from the config.yml
	output:
		o1 = temp( "tempdir/temp_{scenario}_{chunk}/ABCstat.txt"),
		o2 = temp( "tempdir/temp_{scenario}_{chunk}/parameters.txt")
	params:
		vcf = config["vcf"],
		n_sim_per_chunk = config["n_sim_per_chunk"],
		the_script = lambda wildcards: scenario_to_simscript_dict[wildcards.scenario],
		relative_Ne_distribution = config["relative_Ne_distribution"]
	threads: 2
	shell:
		"""
		echo {input[1]} {params.the_script}
		
		mkdir -p tempdir/temp_{wildcards.scenario}_{wildcards.chunk}
		
		cp scripts/{params.the_script} tempdir/temp_{wildcards.scenario}_{wildcards.chunk}/
		
		cat {input[1]} | sed 's/NREPS/{params.n_sim_per_chunk}/g' > tempdir/temp_{wildcards.scenario}_{wildcards.chunk}/argfile.txt
		
		cp {params.vcf}.genepop.txt.bpfile.txt tempdir/temp_{wildcards.scenario}_{wildcards.chunk}/bpfile.txt
		cp {params.vcf}.genepop.txt.spinput.txt tempdir/temp_{wildcards.scenario}_{wildcards.chunk}/spinput.txt
		
		cp {params.relative_Ne_distribution} tempdir/temp_{wildcards.scenario}_{wildcards.chunk}/relative_Ne_distribution.txt
		
		cd tempdir/temp_{wildcards.scenario}_{wildcards.chunk}
		
		bash {params.the_script} bpfile.txt spinput.txt argfile.txt relative_Ne_distribution.txt
		
		# cleanup
		rm {params.the_script} bpfile.txt spinput.txt argfile.txt relative_Ne_distribution.txt progress.log.txt seedms
				
		"""
