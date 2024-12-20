
# PhyloFrame


[![DOI](https://zenodo.org/badge/841574410.svg)](https://doi.org/10.5281/zenodo.14034515)

PhyloFrame uses population level information to correct for unequal ancestral representation in genomic training data. It utilizes population resources such as Gnomad, functional interaction networks like HumanBase, and transcriptomic patient data.  

## Installation

You can install the development version of PhyloFrame like so:

```{r example}
library(githubinstall)
githubinstall("PhyloFrame")

```

## Example: Replication of paper results

This is a basic example which shows you how to recreate a run from the paper as well as to run PhyloFrame on a single datatset.
##### Please ensure you have a directory `data-raw` within the project's home directory.

### 1. Set up data
PhyloFrame utilizes two large data sources: Functional interaction networks and enhanced allele frequencies created from Gnomadv4.1. 

#### Functional Interaction Network
Functional interaction networks used in PhyloFrame are from Humanbase and can be downloaded here: https://hb.flatironinstitute.org/download. We use the mammary epithelium, thyroid, and uterine endometrium networks in the associated paper. We use the Full network in our analysis, however any network with the format below can be used.  In this file the columns are listed on Hummanbase as follows: [entrez gene id 1][entrez gene id 2][posterior prob., with known edges set to 1][posterior prob.] We use the [posterior prob., with known edges set to 1] connection in our analysis and drop the fourth column.

```{network file format}    
5983	1663	1	0.217193
5983	1645	1	0.0558941
5983	2255	1	0.0264535
5983	9401	1	0.513181
```
Once the network is downloaded please ensure it is unzipped and in the `data-raw/` directory within the project's home directory. Then run the script `run_network_annotation.R` with the name of the network file as an argument. For example for a network file titled uterine:
```
Rscript run_network_annotation.R uterine
```
The annotated network it written to a file with same network name with the appended *_symbol.tsv* (ex: uterine_symbol.tsv) in the `data-raw/` directory where it can be read into PhyloFrame.

#### Enhanced Allele Frequencies
Previously calculated enhanced allele frequencies (EAF) are provided with paper's associated Source Data. Please keep in mind this is a large file ~28Gb. EAFs in this file were caluclated with 8 ancestries from Gnomadv4.1 exome files. Please download the file and place it in the `data-raw/gnomad4.1_annotated/` directory within the project's home directory with the file name `gnomad.exomes.all.tsv`. If you would like to calculate your own EAFs please see the section below: [Enhanced Allele Frequency Creation](https://github.com/leslie-smith1112/PhyloFrame/blob/main/README.md#enhanced-allele-frequency-creation).

### Example Run

* Expression matrices for TCGA diseases with classification task are that are used in this project are kept as R-data. To see all associated R-data please load the package and run data(package = "PhyloFrame") * 

```{r example}
library(PhyloFrame)
## basic example code
devtools::load_all()
## For reproduction of PhyloFrame results on multiple TCGA datasets
## ARGUEMENTS: Define disease (breast, uterine, thyroid), name of output directory (will be in results/), and whether would would like to create new training batches
devtools::load_all()
main(breast, TCGA_Breast, FALSE)

## To run PhyloFrame on a single dataset, define a list of training samples (samples not in training set are automatically used as test samples), expression matrix with prediction task column names as "subtype", the name of the directory you want your output in (will be in results/USER_DEFINED_RESULTS_DIR, and the name of the training set (this is used to name output files)
## Here we use the BRCA TCGA samples with EUR ancestry as samples for training data

train_samples <- samples_ancestry$patient[samples_ancestry$consensus_ancestry == "eur"]
single_expr_driver(expression_breast,train_samples, "breast", "TCGA_BRCA", "EUR")
```


## Expected Output 
As PhyloFrame trains on many small, independent, single ancestry datasets, reproduction of PhyloFrame results creates directories seperating each model's ancestry training data and the model number within that ancestry. For example within the user defined results directory, the EUR model 1 results will be in: `model_runs/phyloFrame/eur/model_1`. Within this directory test results on other ancestry samples, as well as the model and its gene signature. Associated benchmark model results will be in `model_runs/benchmark/eur/model_1`.

The results shown in this paper were run on a single core and utilized 92GB. The model runs in approximatley 2 hours and 40 minutes.  
For single runs (runs trained on a single datasets not seperated into smaller dataset batches) results will be in the user defined results directory in directories `phyloFrame/` and `benchmark/`. Wihtin the directory will be the test set results as well as model and model gene signature.

The single model was run on a single core and utilizes 87GB. The model runs in approximatley 1 hour and 40 minutes. 

## Enhanced Allele Frequency Creation 
In this version of PhyloFrame Enhanced Allele Frequencies (EAFs) are calculated from population specific allele frequencies in Gnomadv4.1. In order to calculate your own EAFs, you will need vcf files from whichever population database you are using with allele frequencies for the populations you would like to include in the calculation. An example of the expected vcf format similar is the following: 
```
chr10   45366   rs1554737603    G       C       .       AC0;AS_VQSR   AC=0;AN=30;AF=0.00000;AC_fin=0;AF_fin=0.00000;AN_fin=2;AC_mid=0;AF_mid=0.00000;AN_mid=2;AC_sas=0;AF_sas=0.00000;AN_sas=2  

```
We provide a Snakemake file in the repository that may help you create your own EAF file. The Snakefile and associated scripts will help to get started but will likely need to be heavily edited to reflect things like:
  1. Where you are downloading your vcf files from - *edit Snakefile*
  2. If the chromosomes are in seperate files (currently what the Snakefile expects) or if they are all in a single file - *edit Snakefile*
  3. Which ancestries you are parsing from the vcf - *edit vcf_parsing/main.cpp* AND *code/gnomadV4_EAF_Calculation.R*
  4. Where you want to write output files - *edit Snakefile*

To use your calculated EAF in a PhyloFrame run please make sure it is in the `data-raw/gnomad4.1_annotated` directory with the file name: `gnomad.exomes.all.tsv`. Or change where the file is being read from in the function `load_EAF()` in the `load_large_data.R` script. 
