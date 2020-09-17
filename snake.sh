snakemake --cluster 'bsub -W 4:0 -R "rusage[mem=1000]" -n 2' --jobs 100
