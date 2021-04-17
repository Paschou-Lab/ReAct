
# ReACt 

ReACt(**Re**constructing **A**llelic **C**oun**t**) is a tool built upon our genotype reconstruction framework for case-control GWAS summary statistics.
It includes three modules: Meta-analysis, group GWAS and case-case GWAS.

Please find more details from our **[manuscript on BioRxiv](https://www.biorxiv.org/content/10.1101/2021.04.02.438281v2)**.

All three modules accept tab or space separated summary statistics of case-control GWAS as input, with `SNP`, `CHR`, `BP`, `A1`, `A2`, `OR`/`Beta` and `SE` as mandatory fields, and can be run by specifying a plain text file with designated parameters. 

Some GWAS summary statistics include fields specifying the sample sizes for each SNP. If column headers `nCase` and `nControl` are found in the input, ReAct will use values from these two columns instead of sample sizes specified by parameter file for more accurate results. 

**This is a preliminary implementation _(ReACt beta)_. Please contact us if you identify any bug when using this version of ReACt and we will keep improving. Please also let us know if you have any qestions regarding this readme file. Thank you for understanding.**

For individual module specifics:

[Meta-Analysis](#Meta-analysis)

[Group PRS](#Group-PRS)

[Case-case GWAS](#Case-case-GWAS)

[Examples to run ReAct](#A-few-examples)


_Updated Mar 2021: Added preliminary sample overlap correction for GrpPRS module._

_Updated Apr 17 2021: more robust matrix inverse (svd) in MetaAnalysis module; adjusted iteration steps to accelerate convergence._

# Meta-analysis
## Compilation
**Compilation of MetaAnalysis module requires installation of GNU Scientific Library (GLS). Please find [here](https://www.gnu.org/software/gsl/) for download and installation.**

Once GLS has been installed, download folder `MetaAnalysis_src` for source code of meta-analysis module. Inside the directory with source code, simply run the command:
```
bash Compile /directory/of/where/you/what/the/tool/to/be/
```
Then an excutable titled `MetaAnalysis` shall be created in the specified directory.

## To run MetaAnalysis
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
* **Firth**: By default, [Firth's correction](https://watermark.silverchair.com/80-1-27.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAsUwggLBBgkqhkiG9w0BBwagggKyMIICrgIBADCCAqcGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMWG1Vq2K8Pxv0VOKkAgEQgIICeFDutfgzLqEu5bo1sixxFK6o62rIA6zrnynIjrHtSaPPyrCCjS1pAX-JDE_zK2PzPmgLwggQRu47vB1hxAKW3kSyJiQDqkTuBLGmrI2AyiG3je3qOy16-nRQ_qS62XnvyHrUxAix3ir_nw0GQ6DtZYm9jad93gyJwpZinlojF3SVJTFv1MYO6dJjbtFDC_H27eyoM87lLNADdlQC-5ZdiKsS1LkgOsnpr9zToJsgH__xvewcKYL3UTiie52JWJWKdIEdf23NOSmw_LQEs4JdlVtJ6ZUyUvHpM9kCx0kf26obo_8Hi6Sxoht801CcCtD1xXp9KX-eHRGSDb88uqHOkvvdyHYssA4_e3bhMs0pMSHCbVYUrsgjKki1tKUAhDljIfy3DBldcNtI9cNZqc-kkm2YCz5zA5vpmsRaoXS9y2mKUVGbYfbH5tJnafrvrKrNdsvgnS-dYyho5bhWteJ5_jGdJ5hxElns6Cev_T6Upbl084dM1Gu3KaJcK0UtBcch8akjrnN-R4JHSpkaZEKKW1E6ONmXbM1pJEwFUJgkFatiByhZl8ughPcD-6dxLtf7QHNpqqy03UfBERySDnsJY3Q-Z-ELc6KQM_hEEHwM9PEHUiRc2tOjWwt-4UOsxh4wbrgh33v_RqNTQ1c7KluDgwRA-Xv-LoXgI3w5GT06Cd0MZoKwhV7brOtagGXe7k5BiaHdqBqC3xw37RJ3QaEWjbeP8YHw0RGT-z7e9ESUEiK5rM33jVpQzZXUO06DprzkOeuEYnrQJ8UC6oEdU_pikVAqAN3J0TN8YE6KxmE-3NoYb6sLfhATWel8hhcCa3mCRVRcW1nGLd8o) for logistic regression will be trigered when the difference of sample sized across studies, or between cases and controls, are more than 5 folds (when max sample size/min sample size > 5.0). Specify a positive value for this parameter to change this threshold. **Note: Even though Firth's correction can to certain extent alleviate the convergence problem in rare-event logistic regressions, we don't recommend using ReACt in any case where the difference of sample sized across studies, or between cases and controls, are more than 50 folds. In those cases, we suggest conventional inverse-variance weighted approach, which will give more robust results.**
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

## To run GrpPRS
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

## To run ccGWAS
Similarly, we do:
```
./ccGWAS par.file
```
to run `ccGWAS`. For this module, parameters are almost exactly the same as `MetaAnalysis`, except that it does not accept the **Firth** option. Note that for `ccGWAS`, **Input** takes exactly two comma separated files as input. 

## Output
Output file of `ccGWAS` includes `SNP`, `CHR`, `BP`, `A1`, `A2`, `OR`, `SE`, and `Pval` for case-case association odds ratio, standard error and p-value, and `ControlOR`, `ControlSE` for contorl-control assocaition odds ratio and standard error. 

`ControlOR` is an estimate of the stratification effect for each SNP between two input studies, and `ControlSE` is a measurement for the confidence of this estimate. The lower `ControlSE` is, the more confident we are with the estimate. Therefore, we suggest filter the results by `ControlSE` values. Consider 0.05 as an emperical cutoff.

## A special note regarding the overlap correction for GrpPRS and ccGWAS
We did implement the sample overlap correction from estimation for both modules (same as `MetaAnalysis`, this can be triggered by specifying the **Zthres** parameter). However, we do not recommend using it with `GrpPRS` or `ccGWAS`. Because this scheme attributes estimated overlap of samples into cases and controls proportionally by their sizes, while in reality we would expect the majority of overlap happening in controls rather than cases. This should not have as much impact on meta-analysis, but can bias the results for `GrpPRS` and greatly hurt the power of ccGWAS, since only cases are considered in this analysis. If the exact number of overlap in cases and controls are known, you can specify them through **OverlapControls** and **OverlapCases** for `GrpPRS`, or **CaseInCase**, **ControlInControl** and **CaseInControl** for `ccGWAS`; If not, but you do expect certain amount of cases to be shared, specifying **Zthres** will generally give a little more conservative result for `ccGWAS`; If you expect only controls but not cases to be shared in `ccGWAS`, we sugest not to use **Zthres**. Instead, since overlap in controls can lead to a smaller `ControlSE` in the result, we suggest use a more stringent threshold for result filtering.

# A few examples
The commands
```
echo -e "Input\tInputSumStat1,InputSumStat2
CaseInCase\t2000,0,0,2000
CaseInControl\t0,0,0,0
ControlInControl\t2000,0,0,2000
Output\tOutputFile" > par.metaanalysis
```
will give you a parameter file `par.metaanalysis` that looks like

```
Input   InputSumStat1,InputSumStat2
CaseInCase      2000,0,0,2000
CaseInControl   0,0,0,0
ControlInControl        2000,0,0,2000
Output  OutputFile
```
then you can run 

```
./MetaAnalysis par.metaanalysis
```
Similarly, you can do 
```
echo -e "Target\tTargetSumStat1,TargetSumStat2
Base\tBaseSumStat
Output\tOutputFile
Pthres\t0.01
nCase\t1000,2000
nControl\t1000,2000
nBase\t3000,3000
OverlapCases\t100,200
OverlapControls\t200,400" > par.grpprs

./GrpPRS par.grpprs
```
and 
```
echo -e "Input\tInputSumStat1,InputSumStat2
CaseInCase\t2000,0,0,2000
CaseInControl\t0,0,0,0
ControlInControl\t2000,0,0,2000
Output\tOutputFile" > par.ccgwas

./ccGWAS par.ccgwas
```
