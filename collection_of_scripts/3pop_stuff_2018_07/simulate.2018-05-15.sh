# simulate.sh

# This script executes the following pipeline:
# draw_priors | simulate_coalescent_samples | record_summary_statistics

# it is designed for an LSB batch submission system on a computer cluster, requiring the variable $LSB_JOBINDEX

# example use: bash simulate.2018-05-15.sh Cochlearia.inputFormat.noExcessiveContigs.2018-04-10.vcf.gz.genepop.txt.bpfile.txt Cochlearia.inputFormat.noExcessiveContigs.2018-04-10.vcf.gz.genepop.txt.spinput.txt Cochlearia.test

# toolbase is the ABSOLUTE path to where you keep the scripts
toolbase=/cluster/project/gdc/people/schamath/Cochlearia/tools


# collect trailing arguments and validate
bpfile=$1
spinput=$2
scenario=$3

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

argfile=argfile.${scenario}.txt

if [ -f ./argfile.${scenario}.txt ]
then
	echo "argfile.${scenario}.txt found, proceeding."
else
	echo "argfile.${scenario}.txt not found, dying. Make an argfile.SCENARIO.txt and give SCENARIO as 3rd trailing argument."
	exit
fi


# common msnsam string to all models (drop of parameters is achieved by setting them zero or equal to a redundant one)
# verbose of msnsam_model_cmd
v="-t tbs[theta] -r tbs[rec_rate] tbs[nsites] -I 3[npop] tbs[nsam1] tbs[nsam2] tbs[nsam3] 0[symm_migration_matrix]
-en 0.0[present time] 1[popID] tbs[N1_recent] 
-eg 0.0[present time] 1[popID] tbs[growth_rate_1 !unfree: constant OR exponential rate given by N1_recent and N1_ancient; controlled in argfile!] 
-en 0.0[present time] 2[popID] tbs[N2_recent] 
-eg 0.0[present time] 2[popID] tbs[growth_rate_2 !unfree: constant OR exponential rate given by N2_recent and N2_ancient; controlled in argfile!]
-en 0.0[present time] 3[popID] tbs[N3_recent]
-eg 0.0[present time] 3[popID] tbs[growth_rate_3 !unfree: constant OR exponential rate given by N3_recent and N3_ancient; controlled in argfile!] 
-en tbs[T_N1_ancient] 1[popID] tbs[N1_ancient]
-en tbs[T_N2_ancient] 2[popID] tbs[N2_ancient]
-en tbs[T_N3_ancient] 3[popID] tbs[N3_ancient]
-em tbs[T_m12_start] 1[recipient_ID] 2[donor_ID] tbs[m12 !given by hyper-parameters m12_scale and m12_prop_mig!]
-em tbs[T_m12_stop] 1[recipient_ID] 2[donor_ID] 0[mig_rate]
-em tbs[T_m13_start] 1[recipient_ID] 3[donor_ID] tbs[m13 !given by hyper-parameters m13_scale and m13_prop_mig!]
-em tbs[T_m13_stop] 1[recipient_ID] 3[donor_ID] 0[mig_rate]
-em tbs[T_m21_start] 2[recipient_ID] 1[donor_ID] tbs[m21 !given by hyper-parameters m21_scale and m21_prop_mig!]
-em tbs[T_m21_stop] 2[recipient_ID] 1[donor_ID] 0[mig_rate]
-em tbs[T_m23_start] 2[recipient_ID] 3[donor_ID] tbs[m23 !given by hyper-parameters m23_scale and m23_prop_mig!]
-em tbs[T_m23_stop] 2[recipient_ID] 3[donor_ID] 0[mig_rate]
-em tbs[T_m31_start] 3[recipient_ID] 1[donor_ID] tbs[m31 !given by hyper-parameters m31_scale and m31_prop_mig!]
-em tbs[T_m31_stop] 3[recipient_ID] 1[donor_ID] 0[mig_rate]
-em tbs[T_m32_start] 3[recipient_ID] 2[donor_ID] tbs[m32 !given by hyper-parameters m32_scale and m32_prop_mig!]
-em tbs[T_m32_stop] 3[recipient_ID] 2[donor_ID] 0[mig_rate]
-ej tbs[T_merge_2_into_1] 2[popID] 1[popID]
-en tbs[T_merge_2_into_1] 1[popID] tbs[Nanc12 !constant!]
-eM tbs[T_merge_2_into_1] 0[symm_migration_matrix, i.e. all migration stops at tau_merge_2_into_1]
-ej tbs[T_merge_anc12_into_3] 1[popID] 3[popID]
-en tbs[T_merge_anc12_into_3] 3[popID] tbs[Nanc123 !constant!]"



# the actual (non-verbose) msnsam_model_cmd: copy verbose version and then in TextWrangler grep-replace \[.*?\] with nothing: 
msnsam_model_cmd="-t tbs -r tbs tbs -I 3 tbs tbs tbs 0
-en 0.0 1 tbs 
-eg 0.0 1 tbs 
-en 0.0 2 tbs 
-eg 0.0 2 tbs
-en 0.0 3 tbs
-eg 0.0 3 tbs 
-en tbs 1 tbs
-en tbs 2 tbs
-en tbs 3 tbs
-em tbs 1 2 tbs
-em tbs 1 2 0
-em tbs 1 3 tbs
-em tbs 1 3 0
-em tbs 2 1 tbs
-em tbs 2 1 0
-em tbs 2 3 tbs
-em tbs 2 3 0
-em tbs 3 1 tbs
-em tbs 3 1 0
-em tbs 3 2 tbs
-em tbs 3 2 0
-ej tbs 2 1
-en tbs 1 tbs
-eM tbs 0
-ej tbs 1 3
-en tbs 3 tbs"

	
module load python/2.7

ori_dir=..

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

# make sure spinput also contains correct nreps from argfile: this is only relevant for correct output of the progress report (the simulations would run fine without)
head -n -3 spinput.txt > firstpart
tail -2 spinput.txt > lastpart
echo $nreps > middlepart
cat firstpart middlepart lastpart > spinput.txt
rm firstpart middlepart lastpart


# the pipeline itself
python $toolbase/draw_ms_priors.Cochlearia.2018-05-15.py -bpfile $bpfile -argfile $argfile | $toolbase/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd | python $toolbase/ms2stats.arbitrary_npop.counter.py | python $toolbase/ms2stats.arbitrary_npop.stats.py





