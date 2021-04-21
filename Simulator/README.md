# SimEigentrat
This is the Simulator we used for all our experiments on simulation data. It generates individual level genotypes in plink .tped and .tfam files under a Balding-Nichols model. 

**To get summary statistics to run ReACt, you will need to install [plink](https://www.cog-genomics.org/plink/) and use command --assoc or --logistic to get the summary statistics to be used as input. Summary statistics used in our manuscript are generated through a standard manner using plink.**

## Compilation
Same as all other ReACt modules, simply download the folder and do 
```
bash Compile /directory/of/where/you/what/the/tool/to/be/
```
then an exicutable titled 'SimEigentrat' should be created in the designated directory.

By specifying flag 'Homo' or 'Heter', it can be used to generate homogeneous and hetergeneous populations, under different level of population stratification. More specifically:

To generate homogeneous populations (simulations used for meta-analysis and group PRS), we can run 
```
./SimEigentrat Homo seed nPop Fst nSNP r Ngroup nCausal i output
```
with parameters as below:
* **seed**: a sequence of seed for pseudo randomness generating
* **nPop**: number of populations we need
* **Fst**: predefined level of population stratification (> 0)
* **nSNP**: total number of SNPs to be generated
* **r**: predefined odds risk for causal variants 
* **Ngroup**: number of case/control in each population to be generated. In this model, nPop populations of same size will be generated.
* **nCausal**: Prefefined number of causal SNPs
* **i**: integer experiment trail label, will show up in the SNP ids in the output
* **output**: prefix for output file
* In this case `Homo` is fixed to be the first parameter, which is a model flag for the simulator

To generate hetergeneous populations (simulations used for cc-GWAS), we can run 
```
./SimEigentrat Heter seed 2 Fst nSNP r nCausalShr Ngroup nCausal Ngroup nCausal i output
```
It takes more parameters than generating the homogeneous populations, with the extra parameters as below:
* **nCausalShr**: predefined number of causal SNPs **shared** between populations (the stress test SNPs in our manuscript)
* **Ngroup**: number of case/control in each population to be generated
* **nCausal**: Prefefined number of causal SNPs _exclusive for each population_ (the trait differencial SNPs in our manuscript)
* In this case `Heter` is fixed to be the first parameter, which is a model flag for the simulator
* Also for simulation for cc-GWAS, `2` is fixed to be the parameter for **nPop**
Under this flag, size of populations and number of causal variants in each population can be different, therefore parameters `Ngroup` and `nCausal` are repeated for twice when two populations are to be generated. The number of case/control in each population and number of causal SNPs can be specified for each population individually.



