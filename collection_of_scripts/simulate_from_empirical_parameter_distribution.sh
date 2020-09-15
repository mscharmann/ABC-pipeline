# getting the inputs right:
# bpfile
# spinput file

# the way to get it done:
# open all scenario scripts, argfiles, the spinput and this master script at once!

# set the dates in all filenames and references thereto? (scenario scripts)?
# set msnsam_arg2 ? 5124*10000 = 51240000
# set n_datasets in the argfiles? 10000
# set n_datasets in the spinput? 10000

## convert bpfiles to fixed theta (µ = 1e-08 per site per generation, Ne = 100000, locus_theta = 4*Ne*µ*nsites)
# with open("final_filtered.allsites.pres0.9.minDP6.vcf.genepop.txt.bpfile.txt", "r") as INFILE:
# 	head = INFILE.readline().strip("\n")
# 	locilengths = INFILE.readline().strip("\n").split("\t")
# 	sp1_nsam = INFILE.readline().strip("\n")
# 	sp2_nsam = INFILE.readline().strip("\n")
# 	thetas = INFILE.readline().strip("\n")
# 	rholine = INFILE.readline().strip("\n") 
# 
# thetas_new = [ str( 1e-08 * 4 * 100000 * int(x) ) for x in locilengths ]
# 
# 
# outlines = [head, "\t".join(locilengths) ,sp1_nsam, sp2_nsam,  "\t".join(thetas_new), rholine]
# 
# with open("fixed_theta_bpfile.txt", "w") as OUTFILE:
# 	OUTFILE.write("\n".join(outlines))

LSB_JOBINDEX=1
ori_dir=..
bpfile=Contig_removed_indels_names.fasta.genepop.txt.bpfile.txt
spinput=Contig_removed_indels_names.fasta.genepop.txt.spinput.txt
empdistrfile=param_estim_results.Dianthus_exponential_growth.posterior_distributions.txt

toolbase=/cluster/home/schamath/coalescent_simulations/ABC_tools

#msnsam_arg2=51240000
# 1000*5000 sims on each process on srv1
#msnsam_arg2=10000000

# common msnsam string to all models (drop of parameters is achieved by setting them zero or equal to a redundant one)
msnsam_model_cmd="-t tbs -r tbs tbs -I 3 tbs tbs 0 0
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
-en tbs 3 tbs"

#for scenario in allo_heterogenous allo_homogenous island_heterogenous_and_ghost island_heterogenous island_homogenous_and_ghost island_homogenous iso sym_heterogenous sym_homogenous ;
#for scenario in iso island_hom island_hetbin sym_hom sym_hetbin allo_hom allo_hetbin ;
scenario=Dianthus_exponential_growth
argfile=argfile.${scenario}.txt

module load python/2.7

mkdir temp_${scenario}_${LSB_JOBINDEX}
cd temp_${scenario}_${LSB_JOBINDEX}

cp ${ori_dir}/${bpfile} ./
cp ${ori_dir}/${argfile} ./
cp ${ori_dir}/${spinput} ./
cp ${ori_dir}/${empdistrfile} ./

mv $spinput spinput.txt

# build msnsam_arg2 from the arg and spinput files (nloci * nreps); cannot be a tbs argument to msnsam!:
nloci=$(head -2 spinput.txt | tail -1 )
nreps=$( cat $argfile | grep "nreps" | awk 'FS=" = " {print $3}' )
msnsam_arg2=$(( $nloci * $nreps ))

# make sure spinput also contains correct nreps from argfile: this is actually only relevant for correct output of the progress report:
head -n -3 spinput.txt > firstpart
tail -2 spinput.txt > lastpart
echo $nreps > middlepart
cat firstpart middlepart lastpart > spinput.txt
rm firstpart middlepart lastpart

python /cluster/project/gdc/people/schamath/Dianthus_CEN_introgression_test/tools_used/draw_ms_priors_from_empirical_distribution.2018-02-15.py -bpfile $bpfile -argfile $argfile -emp_distr ${empdistrfile} | $toolbase/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd | python $toolbase/ms2stats.counter.py | python $toolbase/ms2stats.2pop_single_locus.stats.2018-02-14.py



