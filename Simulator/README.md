# SimEigentrat
This is the Simulator we used for all our experiments on simulation data. It generates individual level genotypes in plink .tped and .tfam files under a Balding-Nichols model. 

**To get summary statistics to run ReACt, you will need to install [plink](https://www.cog-genomics.org/plink/) and use command --assoc or --logistic to get the summary statistics to be used as input. Summary statistics used in our manuscript are generated through a standard manner using plink.**

## Compilation
Same as all other ReACt modules, simply download the folder and do 
```
bash Compile /directory/of/where/you/what/the/tool/to/be/
```
then an exicutable titled `SimEigentrat` should be created in the designated directory.

## Quick demo for SimEigentrat
Go to the directory where `SimEigentrat` is created, and run 
```
./SimEigentrat Homo 12345 5 0.05 ${nSNP} ${r} ${Ngroup} ${nCausal} ${i} ${data}/Test0${i}
```



## Running SimEigentrat
By specifying flag 'Homo' or 'Heter', it can be used to generate homogeneous and hetergeneous populations, under different level of population stratification. More specifically:

To generate homogeneous populations (simulations used for meta-analysis and group PRS), we can run 
```
./SimEigentrat Homo seed nPop Fst nSNP r Ngroup nCausal i output
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

_Note that for this simulator, `Fst` has to be a value between 0 and 1. If we want to simulate data without any stratification, we will need to simulate a larger population and split it to get two sub populations without stratification._

To generate hetergeneous populations (simulations used for cc-GWAS), we can run 
```
./SimEigentrat Heter seed 2 Fst nSNP r nCausalShr Ngroup nCausal Ngroup nCausal i output
```
It takes more parameters than generating the homogeneous populations, with the extra parameters as below:
* **nCausalShr**: predefined number of causal SNPs **shared** between populations (the stress test SNPs in our manuscript, we used nCausalShr = 49,000)
* **Ngroup**: number of case/control in each population to be generated
* **nCausal**: Prefefined number of causal SNPs _exclusive for each population_ (the trait differencial SNPs in our manuscript, we used nCausal = 1000 for each population, so number of differencial SNPs was 2000 in total)
* In this case `Heter` is fixed to be the first parameter, which is a model flag for the simulator
* Also for simulation for cc-GWAS, `2` is fixed to be the parameter for **nPop**
Under this flag, size of populations and number of causal variants in each population can be different, therefore parameters `Ngroup` and `nCausal` are repeated for twice when two populations are to be generated. The number of case/control in each population and number of causal SNPs can be specified for each population individually.



