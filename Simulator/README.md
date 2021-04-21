# SimEigenstrat
This is the Simulator we used for all our experiments on simulation data. It generates individual level genotypes in plink .tped and .tfam files under a Balding-Nichols model. We named it `SimEigenstrat` because it is the same simulation model used in [this paper](https://www.nature.com/articles/ng1847).

**To get summary statistics to run ReACt, you will need to install [plink](https://www.cog-genomics.org/plink/) and use command --assoc or --logistic to get the summary statistics to be used as input. Summary statistics used in our manuscript are generated through a standard manner using plink.**

## Compilation
Same as all other ReACt modules, simply download the folder and do 
```
bash Compile /directory/of/where/you/what/the/tool/to/be/
```
then an executable titled `SimEigenstrat` should be created in the designated directory.

## Quick demo for SimEigenstrat
Go to the directory where `SimEigenstrat` is created, and run 
```
./SimEigenstrat Homo 12345 5 0.05 100000 1.2 1000 1000 1 DemoSim_Homo
```
This should not take more than 10 minutes to run on any "normal" desktop or laptop. As a reuslt, it should give us two files 
`DemoSim_Homo.tfam` and `DemoSim_Homo.tped`. They are standard plink plain text files. They should contain genotypes of 100,000 SNPs (out of which 1000 are causal with predefined risk r = 1.2), for 5 stratified populations with fixation index Fst = 0.05, each with 1000 case and 1000 controls -- so 1000*2*5 = 10,000 samples,a and 100,000 SNPs in total. You can verified this by doing
```
$ wc -l DemoSim_Homo.tfam
$ wc -l DemoSim_Homo.tped
```
These should give 
```
10000 DemoSim_Homo.tfam
100000 DemoSim_Homo.tped
```
The files look like
```
$ head -5 DemoSim_Homo.tfam
Pop-1	Sample-ca1	0	0	0	2
Pop-1	Sample-ca2	0	0	0	2
Pop-1	Sample-ca3	0	0	0	2
Pop-1	Sample-ca4	0	0	0	2
Pop-1	Sample-ca5	0	0	0	2
```
So you can easily seperate the populations and cases and controls, it you wish, using plink. DemoSim_Homo.tped will be too wide to show, but it should look like standard plink .tped file with 1 and 2 as genotype encoding. With that, seperate the file by populations, and run standard `--assoc` or `--logistic` (if introduced sample overlap) will give use the input summary statistics we can use for running ReACt. We used this command (with variable parameters for Fst and r, see below) to generate our simulation for experiments for Meta-analysis and group PRS.

Run 
```
./SimEigenstrat Heter 12345 2 0.05 100000 1.2 49000 1000 1000 1000 1000 1 DemoSim_Heter
```
This should be faster than the previous line, as a smaller dataset is being simulated. It should give us `DemoSim_Heter.tfam` and `DemoSim_Heter.tped`. The output should contain 100,000 SNPs (out of which 1000 are causal only in population 1, 1000 are causal only in population 2, 49,000 are causal for both population, all with predefined risk r = 1.2) and only 2 populations, each with 1000 cases and 1000 controls -- so in total 100,000 SNPs and 4000 samples. They are in the same plink-friendly format as the `DemoSim_Homo.tfam` and `DemoSim_Homo.tped`. This is how we generate simulation for our ccGWAS experiments. 

## Parameters to run SimEigenstrat
By specifying flag 'Homo' or 'Heter', it can be used to generate homogeneous and hetergeneous populations, under different level of population stratification. More specifically:

To generate homogeneous populations (simulations used for meta-analysis and group PRS), we can run 
```
./SimEigenstrat Homo seed nPop Fst nSNP r Ngroup nCausal i output
```
with parameters as below:
* **seed**: a sequence of seed for pseudo randomness generating
* **nPop**: number of populations we need (we used nPop = 5, while most of our experiments involved only two)
* **Fst**: predefined level of population stratification (we tested Fst = 0.01/0.05/0.1)
* **nSNP**: total number of SNPs to be generated (we used 100,000)
* **r**: predefined odds risk for causal variants (we tested r = 1.15/1.2/1.3)
* **Ngroup**: number of case/control in each population to be generated. In this model, nPop populations of same size will be generated (we used Ngroup = 1000 upon simulation, and by introducing sample overlap the final sample size will be different depending on the experiment. Please see the manuscript for more details.)
* **nCausal**: Prefefined number of causal SNPs (we used nCausal = 1000. In this model all SNPs are simulated independently, so no matter what number of causal SNPs tested, trend of the results should be similar)
* **i**: integer experiment trail label, this will show up in the SNP ids in the output
* **output**: prefix for output file
* In this case `Homo` is fixed to be the first parameter, which is a model flag for the simulator

For this model, in the simulated SNPs, the first _nCausal_ SNPs will be causal with predefined risk, and the rest _nSNP-nCausal_ will be null SNPs with risk = 1.0.

_Note that for this simulator, `Fst` has to be a value between 0 and 1. If we want to simulate data without any stratification, we will need to simulate a larger population and split it to get two sub populations without stratification._

To generate hetergeneous populations (simulations used for cc-GWAS), we can run 
```
./SimEigenstrat Heter seed 2 Fst nSNP r nCausalShr Ngroup1 nCausal1 Ngroup2 nCausal2 i output
```
It takes more parameters than generating the homogeneous populations, with the extra parameters as below:
* **nCausalShr**: predefined number of causal SNPs **shared** between populations (the stress test SNPs in our manuscript, we used nCausalShr = 49,000)
* **Ngroup1/2**: number of case/control in each population to be generated
* **nCausal1/2**: Prefefined number of causal SNPs _exclusive for each population_ (the trait differencial SNPs in our manuscript, we used nCausal = 1000 for each population, so number of differencial SNPs was 2000 in total)
* In this case `Heter` is fixed to be the first parameter, which is a model flag for the simulator
* Also for simulation for cc-GWAS, `2` is fixed to be the parameter for **nPop**
Under this flag, size of populations and number of causal variants in each population can be different, therefore parameters `Ngroup` and `nCausal` are repeated for twice when two populations are to be generated. The number of case/control in each population and number of causal SNPs can be specified for each population individually.

For this model, in the simulated SNPs, the first _nCausalShr_ SNPs will be causal with predefined risk **in both populations**, followed by _nCausal1_ SNPs having effect `r` in only population 1,  and _nCausal2_ SNPs having effect `r` in only population 2; the rest _nSNP-nCausalShr-nCausal1-nCausal2_ will be null SNPs with risk = 1.0 in both populations.



