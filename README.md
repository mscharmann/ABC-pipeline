# ABC-pipeline

This is a workflow and collection of scripts to estimate the demographic history of populations from sequencing data. It is computationally feasible for datasets from 2 to c. 15,000 loci, if the loci are relatively short (such as <5kb). This includes RNA-seq, ddRAD and other reduced-representation datasets, or WGS data if subsampled to a number of relevent regions.

Almost all steps and scripts can be modified and adopted to test many different questions. This pipeline is derived from the pioneering work of Camille Roux (LINK). My focus was to allow non-programmers like myself to more easily modify the pipeline for particular projects. Hence important parts are written in python (2.7 and/or 3.x). Although code readability and flexibility come at the cost of computational efficiency, it may overall be more efficient to spend more CPU hours but less hands-on working/coding hours.

## Main steps of the workflow
- convert observed data to pipeline inputs: genotypes and fixed model parameters of the sampling scheme (number of populations and samples, number of loci, missingness of loci, length of loci in bp)
- calculate summary statistics for the observed data
- define additional fixed model parameters: spontaneous mutation rate per generation per site, spontaneous recombination rate per generation per site
- define models, free model parameters, and free parameter prior distributions
- simulate the coalescent under the models and parameters chosen, record summary statistics
- compare fit of observed data to the models (model choice)
- fit parameters of one or more models to the observed data (parameter estimation)
- sanity-check the joint posterior distribution of model parameters (posterior-predictive checks)



# Example Workflow & Tutorial

## 0. installation / satisfy dependencies

install conda / bioconda

## 1. convert observed data to pipeline inputs
clone this repository
```
git clone https://github.com/mscharmann/ABC-pipeline

```

make a conda environment that contains snakemake, numpy and scipy
```
cd ABC-pipeline
conda env create --file ABC_sims.yml
conda activate ABC_sims
```

a quick example of the pipeline can then be run locally like so:

```
snakemake --cores 10
```

else on SLURM cluster can do:
```
snakemake --cluster "sbatch" --jobs 100
```

else on LSF cluster can do (give 2 cores per job for 2- and 3-pop models, give only 1 core per job for 1-pop models):
```
snakemake --cluster 'bsub -W 4:0 -R "rusage[mem=1000]" -n 2' --jobs 100
```

or else on LSF better submit a master job for snakemake, named snake.sh and containing the command above, and then submit:
```
bsub -J "master" -W 12:0 -n 1 -R "rusage[mem=2000]" -M 48000000 < snake.sh
```


afterwards, can remove the conda environment
```
conda deactivate
conda remove --name ABC_sims --all
```


### Inputs
- Filtered population genotypes (diploid) in VCF format. (Allo-polyploids may be dealt with as described by the excellent Roux & Pannell paper)

Only bi-allelic SNPs are possible. A highly recommended VCF quality filtering workflow is described here: (Puritz link). An especially important filter is the removal of sites with excessive heterozygosity (Hardy-Weinberg test), which are frequent artefacts in RAD-seq and genome skimming data when paralogous (multiplied) sequences are erroneously collapsed onto one locus in the reference.
	
To characterise the distribution of coalescent times along a genome, many short and unlinked loci are ideal. If you have denovo assembled contigs from RAD-data or otherwise contigs that are <10kb long, you are fine. If you have a proper genome assembly, I recommend to analyse windows of a few hundred bp each. Space the non-overlapping windows such that there is no significant LD between them.

In the VCF, the contigs (column XY) need to be named following a special format: Append to the end of each contig name the length in bp, separated from the rest of the name by an underscore and a capital L. Example: "contigXY_L100". The length of the contigs is an input parameter for simulations, and for the calculation of any per-site summary statistics such as nucleotide diversity pi.

- a population map

This is a tab-separated file with two columns: First column contains names of samples (exactly as in the VCF), the second column the population to which samples are assigned. No header.
sample_1	population_A
sample_2	population_A
sample_3	population B

## execution

## 2. specify models and priors

## 3. ABC in R: model choice

## 4. ABC in R: parameter estimation







