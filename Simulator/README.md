# SimEigentrat
This is the Simulator we used for all our experiments on simulation data. It generates individual level genotypes in plink .tped and .tfam files under a Balding-Nichols model. 

**To get summary statistics to run ReACt, you will need to install [plink](https://www.cog-genomics.org/plink/) and use command --assoc or --logistic to get the summary statistics to be used as input.**

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


