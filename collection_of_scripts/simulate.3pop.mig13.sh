## this is: bash_wrapper.sim3pop.sh
## for EULER

ori_dir=../
bpfile=loci_5000.l_1000bases.50chromperpop.mu_1.rho_1.3pop.bpfile.txt
spinput=loci_5000.l_1000bases.50chromperpop.3pop.spinput.txt


toolbase=/cluster/home/schamath/coalescent_simulations/ABC_tools


# common msnsam string to all models (drop of parameters is achieved by setting them zero or equal to a redundant one)
msnsam_model_cmd_mig23="-t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 2 3 0 -em tbs 3 2 0 -em tbs 2 3 tbs -em tbs 3 2 tbs -ej tbs 2 1 -en tbs 1 tbs -eM tbs 0 -ej tbs 1 3 -en tbs 3 tbs"

# to have the migration between 1 - 3 instead of 2 - 3, we change only the four migration (-em time reciever donor migrant_rate) commandos to ms while keeping all the prior sampling and other scripts identical:
msnsam_model_cmd_mig13="-t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 1 3 tbs -em tbs 3 1 tbs -ej tbs 2 1 -en tbs 1 tbs -eM tbs 0 -ej tbs 1 3 -en tbs 3 tbs"



#for scenario in allo_heterogenous allo_homogenous island_heterogenous_and_ghost island_heterogenous island_homogenous_and_ghost island_homogenous iso sym_heterogenous sym_homogenous ;
#for scenario in iso island_hom island_hetbin sym_hom sym_hetbin allo_hom allo_hetbin ;
scenario=3pop.mig13

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
echo $nreps > middlepart
cat firstpart middlepart lastpart > spinput.txt
rm firstpart middlepart lastpart

python $toolbase/draw_ms_priors.3pop.py -bpfile $bpfile -argfile $argfile  | $toolbase/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd_mig13 | python $toolbase/ms2stats.3pop.counter.py | python $toolbase/ms2stats.3pop.stats.2DSFS.v2.py


##################

