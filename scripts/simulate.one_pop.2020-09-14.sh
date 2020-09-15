# simulate.sh

# This script executes the following pipeline:
# draw_priors | simulate_coalescent_samples | record_summary_statistics

toolbase=../../scripts

# collect trailing arguments and validate
bpfile=$1
spinput=$2
argfile=$3
Ndist=$4

if [ -f ./${bpfile} ]
then
	echo "bpfile found, proceeding."
else
	echo "bpfile (1st trailing argument) not found, dying."
	exit
fi

if [ -f ./${spinput} ]
then
	echo "spinput file found, proceeding."
else
	echo "spinput file (2nd trailing argument) not found, dying."
	exit
fi

if [ -f ./${argfile} ]
then
	echo "argfile found, proceeding."
else
	echo "argfile not found, dying. Make an argfile and give its path as 3rd trailing argument."
	exit
fi

if [ -f ./${Ndist} ]
then
	echo "Ndist found, proceeding."
else
	echo "Ndist not found, dying."
	exit
fi


# common msnsam string to all models (drop of parameters is achieved by setting them zero or equal to a redundant one)
# verbose of msnsam_model_cmd
v="-t tbs[theta] -r tbs[rec_rate] tbs[nsites] -I 1[npop] tbs[nsam1] 0[symm_migration_matrix]
-en 0.0[present time] 1[popID] tbs[N1_recent] 
-eg 0.0[present time] 1[popID] tbs[p1_growth_rate_during_recent_phase !unfree: constant OR exponential rate given by N1_recent and N1_ancient; controlled in argfile!] 
-en tbs[T_N1_interm] 1[popID] tbs[N1_interm]
-en tbs[T_N1_ancient] 1[popID] tbs[N1_ancient]
-L"



# the actual (non-verbose) msnsam_model_cmd: copy verbose version and then in TextWrangler grep-replace \[.*?\] with nothing: 
msnsam_model_cmd="-t tbs -r tbs tbs -I 1 tbs 0
-en 0.0 1 tbs 
-eg 0.0 1 tbs 
-en tbs 1 tbs
-en tbs 1 tbs
-L"
	
# build msnsam_arg2 from the arg and spinput files (nloci * nreps); cannot be a tbs argument to msnsam!:
nloci=$(head -2 spinput.txt | tail -1 )
nreps=$( cat $argfile | grep "nreps" | awk 'FS=" = " {print $3}' )
msnsam_arg2=$(( $nloci * $nreps ))

# make sure spinput also contains correct nreps from argfile: this is only relevant for correct output of the progress report (the simulations would run fine without)
head -n -3 spinput.txt > firstpart
tail -2 spinput.txt > lastpart
echo $nreps > middlepart
cat firstpart middlepart lastpart > spinput.txt
rm firstpart middlepart lastpart


# the pipeline itself
# python $toolbase/draw_ms_priors.1_pop.2020-09-14.py -bpfile $bpfile -argfile $argfile -relative_N $Ndist | $toolbase/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd | tee >(python -W ignore $toolbase/ms2stats.single_pop.haplotype_stats.py) | python $toolbase/ms2stats.arbitrary_npop.counter.record_coaltime.py | python -W ignore $toolbase/ms2stats.arbitrary_npop.stats.py

python $toolbase/draw_ms_priors.1_pop.2020-09-14.py -bpfile $bpfile -argfile $argfile -relative_N $Ndist | $toolbase/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd | python $toolbase/ms2stats.arbitrary_npop.counter.record_coaltime.py | python -W ignore $toolbase/ms2stats.arbitrary_npop.stats.py




