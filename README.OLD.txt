# ABC-pipeline

This is a workflow and collection of scripts to estimate the demographic history of populations from sequencing data. Almost all steps and scripts can be modified and adopted to test many different questions. This pipeline is derived from the pioneering work of Camille Roux (LINK). My focus was to allow non-programmers like myself to more easily modify the pipeline for particular projects. Hence important parts are written in python 2.7. Although code readability and flexibility come at the cost of computational efficiency, it may overall be more efficient to spend more CPU hours but less hands-on working/coding hours.

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
- msnsam (downloadlink)
- python 2.7 with the following modules: numpy
- a linux cluster with many CPUs.
The speed limiting part of the pipeline is the calculation of summary statistics by python, which is confined to a single CPU. However, an unlimited number of python instances can be run in parallel. Requirements of I/O, disk space and RAM are very minor.


## 1. convert observed data to pipeline inputs

### Inputs
- Filtered population genotypes (diploid) in VCF format. Polyploids can be dealt with as described by the excellent Roux & Pannell paper!

Only bi-allelic SNPs are possible. A highly recommended VCF quality filtering workflow is described here: (Puritz link). An especially important filter is the removal of sites with excessive heterozygosity (Hardy-Weinberg test), which are frequent artefacts in RAD-seq and genome skimming data when paralogous (multiplied) sequences are erroneously collapsed onto one locus in the reference.
	
To characterise the distribution of coalescent times along a genome, many short and unlinked loci are ideal. If you have denovo assembled contigs from RAD-data or otherwise contigs that are <10kb long, you are fine. If you have a proper genome assembly, I recommend to analyse windows of a few hundred bp each. Space the non-overlapping windows such that there is no significant LD between them.

In the VCF, the contigs (column XY) need to be named following a special format: Append to the end of each contig name the length in bp, separated from the rest of the name by an underscore and a capital L. Example: "contigXY_L100". The length of the contigs is an input parameter for simulations, and for the calculation of any per-site summary statistics such as nucleotide diversity pi.

- a population map

This is a tab-separated file with two columns: First column contains names of samples (exactly as in the VCF), the second column the population to which samples are assigned. No header.
sample_1	population_A
sample_2	population_A
sample_3	population B

### execution

vcf_to_genepop
generates:
- genepop format. Historical reaons - not in itself useful.

genepop_to_ms
generates:
- bpfile
- spinput
- ms-formatted observed population samples (not genotypes anymore at this point, since the diploid gentoypes

bash script to with pipeline "cat" obs-data into scripts
generates: ABCstats.txt for the obs data!!


## 2. specify models and priors
This is the most complicated part of setting up your study, because several scripts need to be carefully adjusted and sanity-checked. For this, you must read python code, and you must understand the ms command line. All necessary information is contained in Hudson's ms manual ( https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13 ). It takes time to get used to but it will pay off.
essential files and scripts:
- draw_priors.py => draws from priors and writes the drawn parameters (a) to STDOUT as arguments to be piped into the coalescent simulator ms, and (b) to a file. This script must be altered so that it correctly reads and samples the prior distributions defined in argfile, and (b) the lines written to STDOUT exactly correspond to the ms commandline.
- simulate.sh =>  contains the actual pipe commands. Must be adjusted to execute the correct draw_priors.py version, argfiles, and most importantly, the ms commandline.
- argfiles: specifies prior distributions

Start by building the ms commandline, contained in "simulate.sh". To adjust and test the command, you have to execute only msnsam but not the rest of "simualte.sh", like so:

```
toolbase=.
$toolbase/msnsam/msnsam
```
In this first test, ms nsam should complain and ask for more arguments, displaying a list of options. Next, give it some nsam parameter value (the total number of chromosomes / samples to be generated), a value for the parameter msnsam_arg2 (nreps = number of repetitions; in this pipeline it corresponds to the number of independet loci * the number of iterations), and the simplest possible population command, i.e. only -t theta. This will generate chromosomes from a single neutral, panmictic, haploid (Wright-Fisher) population of constant size:

```
msnsam_arg2=1
msnsam_model_cmd="-t 4.0"

$toolbase/msnsam/msnsam 10 $msnsam_arg2 $msnsam_model_cmd
```
Yeah, coalescent simulations! But in order to efficiently run many simulations with many different values for these parameters (from prior distributions), they must be provided as variables. Fortunately, msnsam (like ms) accepts so-called 'tbs' arguments, i.e. values that are piped into it from STDIN. They will be read strictly in the supplied order:
```
msnsam_arg2=1
msnsam_model_cmd="-t tbs"

echo 10 4.0 | $toolbase/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd
```
Further parameters must be added in most cases (except you were interested in this simple scenario). Please use the ms manual to craft an ms command for your own models / scenarios.

It is easy to get lost when you have many parameters / arguments. I found it helped a lot to create a "verbose" version of the ms command, e.g. "-t tbs[theta]". Each flag and argument are followed by a description of the parameter in square brackets. Even better, you can put a descriptive name into these brackets corresponding exactly to the name of the model parameter that you later want to vary and estimate using ABC. To make a runnning and testable ms command, copy the verbose version and then, in TextWrangler, grep-replace \[.*?\] with nothing. 

When the ms commandline is doing what you want, you need to adjust the prior generator python script to supply the tbs arguments in correct order. Copy the verbose version of your ms command to the python script to keep an eye on the structure of the STDOUT that is aimed for.
Then, delete everything else except the names of parameters that the script is going to produce: This is the structure of the STDOUT.
At the same time, you have to adjust the parameter names in the python script, and how it reads the prior distributions of the parameters from the argfile.txt

To test the prior generator script that you are modifying, get some bpfile and argfiles and
```
python draw_ms_priors.X.py -argfile argfile.X.txt -bpfile XXX.bpfile.txt

python draw_ms_priors.Cochlearia.2018-05-15.py -argfile argfile.Cochlearia.test.txt -bpfile Cochlearia.inputFormat.noExcessiveContigs.2018-04-10.vcf.gz.genepop.txt.bpfile.txt

```

Edit scripts and inout files in this order:
- add all parameters for which priors exist to the argfile
- add all parameters for which priors exist to the list params_to_be_read in the function read_argfile(), and check that draw_ms_priors.py runs it correctly
- edit the function draw_and_pipe_priors()
	- param_log_header : the header for your parameters log file . It is the same as the list params_to_be_read, except for the categorical switches.
	- then edit the core of the function, i.e. how each of the paramaters is drawn from the specifified prior ranges

## 3. simulate
Now all your scripts are in working condition and ready to run. But you still need to parallelise the simulations - for this there are scripts to start and control many parallel instances of the pipeline

## 4. ABC in R: model choice

## 5. ABC in R: parameter estimation

##






