
# ReACt 

ReACt(**Re**constructing **A**llelic **C**oun**t**) is a tool built upon our genotype reconstruction framework for case-control GWAS summary statistics.
It includes three modules: Meta-analysis, group PRS and case-case GWAS.

Please find more details in our **[manuscript on BioRxiv](https://www.biorxiv.org/content/10.1101/2021.04.02.438281v2)**.

All three modules accept tab or space separated summary statistics of case-control GWAS as input, with `SNP`, `CHR`, `BP`, `A1`, `A2`, `OR`/`Beta` and `SE` as mandatory fields, and can be run by specifying a plain text file with designated parameters. 

Some GWAS summary statistics include fields specifying the sample sizes for each SNP. If column headers `nCase` and `nControl` are found in the input, ReAct will use values from these two columns instead of sample sizes specified by parameter file for more accurate results. 

The implementation has been tested on Linux (CentOS Linux 7 (Core)) and Mac OS (Mac OS X EI Capitan). 
**This is a rather preliminary implementation. Please contact us if you identify any bug when using this version of ReACt and we will keep improving. Please also let us know if you have any questions regarding this readme file. Thank you for understanding.**

For individual module specifics:

**[Meta-Analysis](#Meta-analysis)** -- [fast forward to demo](#A-quick-demo-of-MetaAnalysis)

**[Group PRS](#Group-PRS)** -- [fast forward to demo](#A-quick-demo-of-GrpPRS)

**[Case-case GWAS](#Case-case-GWAS)** -- [fast forward to demo](#A-quick-demo-of-ccGWAS)

We also provided codes for the simulator we used to generate synthetic data (toy input). Please see the readme file in _Simulator_ for more details. 


_Updated Mar 2021: Added preliminary sample overlap correction for GrpPRS module._

_Updated Apr 17 2021: more robust matrix inverse (svd) in MetaAnalysis module; adjusted iteration steps to accelerate convergence._

_Updated Apr 21 2021: Added codes for simulator and toy input files; added more detailed demo for each module._

# Meta-analysis
## Compilation
**Compilation of MetaAnalysis module now requires installation of GNU Scientific Library (GLS). Please find [here](https://www.gnu.org/software/gsl/) for download and installation.**

Once GLS has been installed, download folder `MetaAnalysis_src` for source code of meta-analysis module. Inside the directory with source code, simply run the command:
```
bash Compile /directory/of/where/you/what/the/tool/to/be/
```
Then an excutable titled `MetaAnalysis` shall be created in the specified directory.

## A quick demo of MetaAnalysis
For a fast run of the demo, we suggest:
Download `ReAct-main.zip`, unzip it, and do below for compilation of our Meta-analysis module:
```
cd ReAct-main/MetaAnalysis_src/
bash Compile ../
cd ../
```
So that the `MetaAnalysis` tool will be in the folder `ReAct-main/`, which is in the same layer as folder `ToyInput/`. **Don't forget to install GLS for this module.**

Then we can easily run MetaAnalysis on the toy input example as below:

The commands
```
echo -e "Input\tToyInput/Toy_Meta.In1,ToyInput/Toy_Meta.In2
CaseInCase\t1000,0,0,1000
CaseInControl\t0,0,0,0
ControlInControl\t1000,0,0,1000
Output\tToy_Meta.out" > par.metaanalysis
```
should give us a parameter file `par.metaanalysis` that looks like
```
Input     ToyInput/Toy_Meta.In1,ToyInput/Toy_Meta.In2
CaseInCase      1000,0,0,1000
CaseInControl   0,0,0,0
ControlInControl        1000,0,0,1000
nFiles  2
Output  Toy_Meta.out
```
then simply run 

```
./MetaAnalysis par.metaanalysis
```
Within 10 seconds on most of the PC, we should be able to get two files `Toy_Meta.out` and `Toy_Meta.out.log` (for this toy input the log should be empty), where the results are in `Toy_Meta.out` and it looks like this:
```
$ head Toy_Meta.out
SNP	CHR	BP	A1	A2	nCase	nControl	OR	SE	Pval
rs1.16780	1	8390000	2	1	2000	2000	0.911469	0.047884	5.2883e-02
rs1.30808	1	15440000	2	1	2000	2000	1.013197	0.061327	8.3071e-01
rs1.33756	1	16910000	2	1	2000	2000	1.033609	0.048559	4.9603e-01
rs1.44462	1	22240000	2	1	2000	2000	1.101519	0.048048	4.4183e-02
rs1.48040	1	24000000	2	1	2000	2000	0.898449	0.044769	1.6760e-02
rs1.53484	1	26710000	1	2	2000	2000	0.993419	0.046420	8.8690e-01
rs1.66016	1	32950000	1	2	2000	2000	0.966554	0.046083	4.6039e-01
rs1.72310	1	36110000	1	2	2000	2000	1.013949	0.045959	7.6310e-01
rs1.80739	1	40300000	1	2	2000	2000	0.971235	0.056950	6.0830e-01
```
We can sort it by the `Pval` column, which will give us
```
$ sort -gk10 Toy_Meta.out| head
SNP	CHR	BP	A1	A2	nCase	nControl	OR	SE	Pval
rs1.892	1	455160	2	1	2000	2000	1.360291	0.045132	9.2514e-12
rs1.334	1	161735	2	1	2000	2000	1.379423	0.047857	1.7995e-11
rs1.574	1	285136	1	2	2000	2000	0.745643	0.045391	1.0053e-10
rs1.480	1	239170	1	2	2000	2000	0.740488	0.046510	1.0484e-10
rs1.78	1	37629	2	1	2000	2000	1.328030	0.045554	4.7345e-10
rs1.782	1	395673	2	1	2000	2000	1.325074	0.045449	5.9031e-10
rs1.433	1	212426	1	2	2000	2000	0.714288	0.054959	9.2290e-10
rs1.179	1	85383	1	2	2000	2000	0.755939	0.045727	9.4317e-10
rs1.890	1	453661	2	1	2000	2000	1.312954	0.045357	1.9365e-09
```
In this toy input, SNP rs1.1-rs1.1000 are all predefined causal SNPs with r = 1.2 (Please see readme of our Simulator for more details), so all of the top SNPs here are ture positive. 

_Note that in our manuscript we reported results based on an average out of 10 iterarions, but for any single experiment the trend should be comparable. Please feel free to further verify power and type I error rate under different conditions using the simulator we provided, and compare with [METAL](https://github.com/statgen/METAL) or [ASSET](https://bioconductor.org/packages/release/bioc/html/ASSET.html). This demo corresponds to results in figure 1,2 and table S1, S2, S3 in our manuscript._ 

## To use MetaAnalysis
To run MetaAnalysis, go into the directory where the excutable locates and do:
```
./MetaAnalysis par.file
```
where `par.file` is a plain text file specifying parameters for the analysis. Mandatory parameters are:
* **Input**: Full path of input summary statistics, separated by comma.
* **Output**: Full path of the output file.
* **CaseInCase**: Flatten matrix of case-case sample overlap across studies. This matrix should have in total *N^2* comma separated elements, where *N* is the number of input studies specified in **Input**. It should be symmetric, and diagonal elements of this matrix should correspond to the sample size for cases in the each input study. See below for an example.
* **ControlInControl**: Flatten matrix of control-control sample overlap across studies.
* **CaseInControl**: Flatten matrix of case-control sample overlap across studies (zero matrix if none). This is the only matrix out of the three that can be asymmetric.
Some optional parameters are:
* **Firth**: By default, [Firth's correction](https://watermark.silverchair.com/80-1-27.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAsUwggLBBgkqhkiG9w0BBwagggKyMIICrgIBADCCAqcGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMWG1Vq2K8Pxv0VOKkAgEQgIICeFDutfgzLqEu5bo1sixxFK6o62rIA6zrnynIjrHtSaPPyrCCjS1pAX-JDE_zK2PzPmgLwggQRu47vB1hxAKW3kSyJiQDqkTuBLGmrI2AyiG3je3qOy16-nRQ_qS62XnvyHrUxAix3ir_nw0GQ6DtZYm9jad93gyJwpZinlojF3SVJTFv1MYO6dJjbtFDC_H27eyoM87lLNADdlQC-5ZdiKsS1LkgOsnpr9zToJsgH__xvewcKYL3UTiie52JWJWKdIEdf23NOSmw_LQEs4JdlVtJ6ZUyUvHpM9kCx0kf26obo_8Hi6Sxoht801CcCtD1xXp9KX-eHRGSDb88uqHOkvvdyHYssA4_e3bhMs0pMSHCbVYUrsgjKki1tKUAhDljIfy3DBldcNtI9cNZqc-kkm2YCz5zA5vpmsRaoXS9y2mKUVGbYfbH5tJnafrvrKrNdsvgnS-dYyho5bhWteJ5_jGdJ5hxElns6Cev_T6Upbl084dM1Gu3KaJcK0UtBcch8akjrnN-R4JHSpkaZEKKW1E6ONmXbM1pJEwFUJgkFatiByhZl8ughPcD-6dxLtf7QHNpqqy03UfBERySDnsJY3Q-Z-ELc6KQM_hEEHwM9PEHUiRc2tOjWwt-4UOsxh4wbrgh33v_RqNTQ1c7KluDgwRA-Xv-LoXgI3w5GT06Cd0MZoKwhV7brOtagGXe7k5BiaHdqBqC3xw37RJ3QaEWjbeP8YHw0RGT-z7e9ESUEiK5rM33jVpQzZXUO06DprzkOeuEYnrQJ8UC6oEdU_pikVAqAN3J0TN8YE6KxmE-3NoYb6sLfhATWel8hhcCa3mCRVRcW1nGLd8o) for logistic regression will be trigered when the difference of sample sized across studies, or between cases and controls, are more than 5 folds (when max sample size/min sample size > 5.0). Specify a positive value for this parameter to change this threshold. 
* **Zthres**: Specify a positive value for *Zthres* to **trigger correction of unknow sample overlap**. If this parameter is not included in the `par.file`, `MetaAnalysis` will take the exact number of sample overlap specified by **CaseInCase**, **ControlInControl** and **CaseInControl** *as is*. In the case where **Zthres** is specified, `MetaAnalysis` will get an estimate for sample overlap using SNPs with |*Z-scores*| < *Zthres*. We adopt this correction sheme from [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation#Sample_Overlap_Correction). Emperically, we suggest choosing *Zthres* = 1.0. Note this parameter will overide the overlap sample sizes provided in **CaseInCase**, **ControlInControl** and **CaseInControl**. Therefore, when the exact sample overlap is known, do not include this parameter in `par.file`.

### Specifying sample size in CaseInCase, ControlInControl and CaseInControl: An example
Suppose we are meta-analyzing three studies, with 1000 cases 1000 controls in study 1; 2000 cases 2000 controls in study 2; and 3000 cases 3000 controls in study 3.
Meanwhile, there are 100 cases shared by study 1 and study 2, 500 controls shared by study 2 and study 3, and 20 cases in study 1 appearing as controls in study 3. 
Then in this case, we shall have three sample overlapping matrices as below:

`Case-Case`

![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Bpmatrix%7D%201000%20%26%20100%20%26%200%20%5C%5C%20100%20%26%202000%20%26%200%20%5C%5C%200%20%26%200%20%26%203000%20%5C%5C%20%5Cend%7Bpmatrix%7D)

`Control-Control`

![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Bpmatrix%7D%201000%20%26%200%20%26%200%20%5C%5C%200%20%26%202000%20%26%20500%20%5C%5C%200%20%26%20500%20%26%203000%20%5C%5C%20%5Cend%7Bpmatrix%7D)

and `Case-Control`

![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Bpmatrix%7D%200%20%26%200%20%26%2020%20%5C%5C%200%20%26%200%20%26%200%20%5C%5C%200%20%26%200%20%26%200%20%5C%5C%20%5Cend%7Bpmatrix%7D)

Then in this case, we should have below in `par.file`:
```
CaseInCase  1000,100,0,100,2000,0,0,0,3000
ControlInControl  1000,0,0,0,2000,500,0,500,3000
CaseInControl 0,0,20,0,0,0,0,0,0
```

## Output
Output file of `MetaAnalysis` includes `SNP`, `CHR`, `BP`, `A1`, `A2`, `nCase` for total number of cases, `nControl` for total number of controls, `OR`, `SE`, and `Pval` for meta-analyzed odds ratio, standard error and p-value.

# Group PRS
## Compilation
Download folder `GrpPRS_src` for source code of group PRS module. Inside the directory with source code, simply run the command:
```
bash Compile /directory/of/where/you/what/the/tool/to/be/
```
Then an excutable titled `GrpPRS` shall be created in the specified directory.

## A quick demo of GrpPRS
Again, for a fast run we suggest:
Download `ReAct-main.zip`, unzip it, and do below for compilation of our Group PRS module:
```
cd ReAct-main/GrpPRS_src/
bash Compile ../
cd ../
```
So that the `GrpPRS` tool will be in the folder `ReAct-main/`, which is in the same layer as folder `ToyInput/`.

Then we can run GrpPRS on the toy input example as below:

The commands
```
echo -e "Target\tToyInput/Toy_GrpPRS.target
Base\tToyInput/Toy_GrpPRS.base
Output\tToy_GrpPRS.out
Pthres\t1e-5
nCase\t1000
nControl\t1000
nBase\t1000,1000
OverlapCases\t0
OverlapControls\t0" > par.grpprs
```
should give us a parameter file `par.grpprs` that looks like
```
Target  ToyInput/Toy_GrpPRS.target
Base    ToyInput/Toy_GrpPRS.base
Output  Toy_GrpPRS.out
Pthres  1e-5
nCase   1000
nControl        1000
nBase   1000,1000
OverlapCases    0
OverlapControls 0
```
then we can run 

```
./GrpPRS par.grpprs
```
For this we should get two files `Toy_GrpPRS.out` and `Toy_GrpPRS.out.log`. Main results are in `Toy_GrpPRS.out`. In this case it looks like this:
```
$ head Toy_GrpPRS.out
InFile	Pthres	nSNPs	CasePRS	ControlPRS	CasePRS_SE	ControlPRS_SE	R2	Pval
ToyInput/Toy_GrpPRS.In1	0.000010	36	0.016312	0.004217	0.017272	0.017407	0.108341	9.5673e-52
```
Note that this toy example is based on a simulation with 1000 causal SNPs shared between base and target, each with a predefined risk r = 1.2 (which is very strong), so we are seeing a very visible seperation between case and control scores from the Pvalue. The log file for toy input should look like this:
```
$ head Toy_GrpPRS.out.log 
Analysis Starts.
P value threshold for base SNPs : 1.00e-05.
36 SNPs below P threshold read from base.
Study ToyInput/Toy_GrpPRS.In1 Finished, 36 SNPs taken for PRS computation.
```

_If you would like to compare the results with [PRSice2](https://www.prsice.info/), please use our simulator to generate individual level genotype. This demo corresponds to results in table 2 in our manuscript._ 

## To use GrpPRS
Similar to MetaAnalysis, go into the directory where the `GrpPRS` locates and run:
```
./GrpPRS par.file
```
Mandatory parameters for `GrpPRS` are:
* **Base**: Full path of the base summary statistcs.
* **Target**: Full path of the target summary statistcs, separated by comma (Supporting group PRS computation for multiple target studies simultaneously).
* **Output**: Full path of the output file.
* **Pthres**: Using SNPs with *p* < *Pthres* in the base summary statistics for PRS computation. Default *Pthres* = 1. If not specified, all SNPs from base will be used.
* **nCase**: Array of case sample size for each target study. Lenth of this array should be the same as number of target sudies specified. 
* **nControl**: Array of control sample size for each target study. Lenth of this array should be the same as number of target sudies specified. 
* **nBase**: Number of cases and controls in the base study, separated by comma
Some optional parameters are:
* **OverlapCases**: Sample overlap for cases between each target study and the base study, separated by comma. Default 0 for all if not specified
* **OverlapControls**: Sample overlap for controls between each target study and the base study, separated by comma. Default 0 for all if not specified
* **Zthres**: Same as **Zthres** for `MetaAnalysis`, triggers correction for unknow sample overlap. Please see [this note](#A-special-note-regarding-the-overlap-correction-for-GrpPRS-and-ccGWAS) before using.

## Output
Output file of `GrpPRS` includes `InFile`, `Pthres`, `nSNPs` for the number of SNPs used in the PRS analysis of this base-target pair, `CasePRS`, `ControlPRS`, `CasePRS_SE` and `ControlPRS_SE` for  group mean PRS of cases and controls and their standard errors; `R2` is the *R^2* value converted from t-test statistics, which herefore, corresponds to the regression *R^2* with only the PRS predictor; `Pval` is the t-test p-value comparing case-control PRS distribution. 

## Note for GrpPRS
`GrpPRS` **does not** automatically prune/clump/lasso the base summary statistics. So we suggest user thinning the base file using some external tool before feeding it as an input to `GrpPRS`.


# Case-case GWAS
## Compilation
Download folder `ccGWAS_src` for source code of case-case GWAS module. Inside the directory with source code, simply run the command:
```
bash Compile /directory/of/where/you/what/the/tool/to/be/
```
Then an excutable titled `ccGWAS` shall be created in the specified directory.

## A quick demo of ccGWAS
Similarly, we suggest:
Download `ReAct-main.zip`, unzip it, and do below for compilation of our ccGWAS module:
```
cd ReAct-main/ccGWAS_src/
bash Compile ../
cd ../
```
So that the `ccGWAS` tool will be in the folder `ReAct-main/`, which is in the same layer as folder `ToyInput/`.

With that we can run ccGWAS on the toy input. We first generate parameter file as below:

```
echo -e "Input\tToyInput/Toy_ccGWAS.In1,ToyInput/Toy_ccGWAS.In2
CaseInCase\t2000,0,0,2000
CaseInControl\t0,0,0,0
ControlInControl\t2000,0,0,2000
Output\tToy_ccGWAS.out" > par.ccgwas
```
which is supposed to give us
```
Input   ToyInput/Toy_ccGWAS.In1,ToyInput/Toy_ccGWAS.In2
CaseInCase      2000,0,0,2000
CaseInControl   0,0,0,0
ControlInControl        2000,0,0,2000
Output  Toy_ccGWAS.out
```
then if we do 
```
./ccGWAS par.ccgwas
```
we will get file `Toy_ccGWAS.out` and `Toy_ccGWAS.out.log` (which is again empty on this toy example). The results should look like below:
```
$ wc -l Toy_ccGWAS.out
100001 Toy_ccGWAS.out
$ head Toy_ccGWAS.out
SNP	CHR	BP	A1	A2	OR	SE	Pval	ControlOR	ControlSE
rs1.15416	1	7700000	2	1	0.926037	0.049281	1.1894e-01	1.015164	0.047930
rs1.46761	1	23390000	2	1	0.966564	0.064864	6.0007e-01	0.911951	0.060649
rs1.55154	1	27580000	2	1	0.961273	0.045885	3.8936e-01	1.073760	0.045726
rs1.55373	1	27690000	2	1	1.003955	0.050324	9.3748e-01	1.032559	0.050186
rs1.70279	1	35180000	1	2	0.997313	0.046590	9.5394e-01	0.860469	0.046735
rs1.71259	1	35680000	1	2	1.017759	0.057443	7.5926e-01	1.508643	0.058011
rs1.71741	1	35920000	1	2	1.001902	0.046481	9.6739e-01	0.830002	0.046866
rs1.75279	1	37650000	2	1	0.995411	0.047411	9.2272e-01	1.139804	0.047244
rs1.97843	1	48960000	1	2	0.959732	0.050914	4.1951e-01	0.515156	0.051061
```
As demonstrated in our manuscript, we suggest filter by the estimated standard deviation of controls. We used 0.05 as a cut off, more specifically, we did:
```
awk '$10<0.05{print $0}' Toy_ccGWAS.out > t && mv t Toy_ccGWAS.out
```
If we look at it now, we can see actually many SNPs will be filtered out under this setting (2000 cases, 2000 controls in each study, predefined r = 1.2). This is why we tested for various sample sizes in our experiment, and why we suggest more controls to be included. If we don't want this to happen, either we include more control samples to get a reliable estimate (smaller ControlSE), or we will have to go lenient on the filtering. 
```
$ wc -l Toy_ccGWAS.out
54703 Toy_ccGWAS.out
```
Now if we sort by the `Pval` column, we get:
```
$ sort -gk8 Toy_ccGWAS.out| head
rs1.50886	1	25436897	2	1	1.381536	0.047724	1.2688e-11	0.363717	0.047411
rs1.50334	1	25165542	2	1	1.368019	0.048733	1.2737e-10	0.799708	0.047389
rs1.50255	1	25123933	1	2	0.749349	0.045594	2.4722e-10	1.667225	0.046662
rs1.50233	1	25112155	1	2	0.744103	0.047230	3.8946e-10	0.930692	0.048594
rs1.50725	1	25354203	1	2	0.760088	0.044730	8.6339e-10	1.275916	0.045050
rs1.50702	1	25340608	1	2	0.759832	0.044886	9.4173e-10	0.964350	0.044984
rs1.50626	1	25303122	2	1	1.314941	0.045056	1.2264e-09	0.919387	0.044754
rs1.50036	1	25017130	2	1	1.333289	0.047671	1.5995e-09	0.978397	0.046487
rs1.8689	1	4354456	1	2	1.335644	0.048424	2.2784e-09	0.961769	0.049879
rs1.50481	1	25231117	1	2	0.753104	0.047711	2.7968e-09	1.254875	0.049409
```
For the setting under which the toy input was simulated, SNP rs1.1-rs1.49000 are stress SNP, SNP rs1.49001-rs1.51000 are trait differencial SNPs, and the rest are null-null SNPs. So in this result, we can see most of the top SNPs picked up are real trait differencial SNPs (true positive), except one of them, rs1.8689, which is a stress testing SNP being falsely picked up. 

_Please feel free to generate more input studies under other conditions using our simulator. This demo corresponds to results in table 5 in our manuscript._ 

## To use ccGWAS
Similarly, we do:
```
./ccGWAS par.file
```
to run `ccGWAS`. For this module, parameters are almost exactly the same as `MetaAnalysis`, except that it does not accept the **Firth** option. Note that for `ccGWAS`, **Input** takes exactly two comma separated files as input. 

## Output
Output file of `ccGWAS` includes `SNP`, `CHR`, `BP`, `A1`, `A2`, `OR`, `SE`, and `Pval` for case-case association odds ratio, standard error and p-value, and `ControlOR`, `ControlSE` for contorl-control assocaition odds ratio and standard error. 

`ControlOR` is an estimate of the stratification effect for each SNP between two input studies, and `ControlSE` is a measurement for the confidence of this estimate. The lower `ControlSE` is, the more confident we are with the estimate. Therefore, we suggest filter the results by `ControlSE` values. Consider 0.05 as an emperical cutoff.

