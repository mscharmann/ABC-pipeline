
vcf=$1
popmap=$2
contig_lengths_file=$3
mu=$4
rho=$5
pshp=$6

toolbase=./scripts

##
python $toolbase/vcf_to_genepop.aribtrary_npop.py --vcf ${vcf} --popfile ${popmap} 

## 
python $toolbase/genepop_to_mscalc.arbitrary_npop.2025-06-03.py -gp ${vcf}.genepop.txt -contig_lengths ${contig_lengths_file} -mu ${mu} -rho ${rho} -haplopops ${pshp}

cp ${vcf}.genepop.txt.spinput.txt spinput.txt
cat ${vcf}.genepop.txt.ms.txt | python $toolbase/ms2stats.arbitrary_npop.counter.py | python -W ignore $toolbase/ms2stats.arbitrary_npop.stats_and_pidistr.py

mv pi_distr_mean_1.txt ${vcf}.pi_distr_mean_1.txt
mv ABCstat.txt ${vcf}.ABCstat.txt
rm progress.log.txt spinput.txt

