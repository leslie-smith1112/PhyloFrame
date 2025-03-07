
# PhyloFrame


[![DOI](https://zenodo.org/badge/841574410.svg)](https://doi.org/10.5281/zenodo.14034515)

PhyloFrame uses population level information to correct for unequal ancestral representation in genomic training data. It utilizes population resources such as Gnomad, functional interaction networks like HumanBase, and transcriptomic patient data.  

## Installation

You can install the development version of PhyloFrame by cloning this repository. Please note that this repository uses git lfs to store package data.

## Replication of paper results

This is a basic example which shows you how to recreate a run from the paper as well as to run PhyloFrame on a single datatset.

### 1. Set up data
PhyloFrame utilizes two large data sources: Functional interaction networks and enhanced allele frequencies for all chromosoems created from Gnomadv4.1. These are not available withing the R package data and must be downloaded and set up seperatley.

Expression data and sample batches used are included in the package data.

#### Functional Interaction Network
Functional interaction networks used in PhyloFrame are from Humanbase and can be downloaded here: https://hb.flatironinstitute.org/download. We use the mammary epithelium, thyroid, and uterine endometrium networks in the associated paper. We use the Full network in our analysis, however any network with the format below can be used.  In this file the columns are listed on Hummanbase as follows: [entrez gene id 1][entrez gene id 2][posterior prob., with known edges set to 1][posterior prob.] We use the [posterior prob., with known edges set to 1] connection in our analysis and drop the fourth column.

```{network file format}    
5983	1663	1	0.217193
5983	1645	1	0.0558941
5983	2255	1	0.0264535
5983	9401	1	0.513181
```
After selecting the networks - copy the links and put it in the hb_links=("") array within the run_network_annotation.sh script. Example: 

```{nerwork links array}
hb_links=("https://s3-us-west-2.amazonaws.com/humanbase/networks/blood.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/adrenal_gland.gz")
```
This script will download the network, placing it in the `data-raw/` folder, unzip it and annotate it. The annotated network is written to a file with same network name with the appended *_symbol.tsv* (ex: uterine_symbol.tsv) in the `data-raw/` directory where it can be read into PhyloFrame.

Please note that memory usage for run_network_annotation.sh is currently set very high for the annotation of all humanbase networks - depending on how many networks you are annotating memory can be lowered.


#### Enhanced Allele Frequencies
Previously calculated enhanced allele frequencies (EAF) are provided with paper's associated Source Data and can be downloaded from 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14180045.svg)](https://doi.org/10.5281/zenodo.14180045). Please keep in mind this is a large file ~28Gb. EAFs in this file were caluclated with 8 ancestries from Gnomadv4.1 exome files. Please download the file and place it in the `data-raw/` directory within the project's home directory with the file name `gnomad.exomes.all.tsv`. If you would like to calculate your own EAFs please see the section below: [Enhanced Allele Frequency Creation](https://github.com/leslie-smith1112/PhyloFrame/blob/main/README.md#enhanced-allele-frequency-creation).

### 2. Example Run

* Expression matrices for TCGA diseases with classification task are that are used in this project are kept as R-data. To see all associated R-data please load the package and run data(package = "PhyloFrame") * 
If you would like to run PhyloFrame on a slurm-based HPC, you can use `run.sh` - edit email and account information as needed. If you would like to make the PhyloFrame calls directly you can do so as follows:

```{r example}
## basic example code

## load PhyloFrame library
devtools::load_all()

## 2.1 For reproduction of PhyloFrame results on multiple TCGA datasets
## ARGUEMENTS: Define disease (breast, uterine, thyroid), name of output directory (will be in `results/`), and whether would would like to create new training batches (to use training batches from paper leave the third argument set to FALSE - sample batch lists are kept in `data-raw/$DISEASE$_samples/` Example: `brca_samples`)
main("breast", "TCGA_Breast", FALSE)

## 2.2 To run PhyloFrame on a single dataset, define a list of training samples (samples not in training set are automatically used as test samples), expression matrix with prediction task column names as "subtype", the name of the directory you want your output in (will be within `results/`, and the name of the training set (this is used to name output files)
## Here we use the BRCA TCGA samples with EUR ancestry samples for training data

train_samples <- ancestry$patient[ancestry$consensus_ancestry == "eur"]
single_expr_driver(expression_breast,train_samples, "breast", "TCGA_BRCA", "EUR")
```


## Expected Output 
As PhyloFrame trains on many small, independent, single ancestry datasets, reproduction of PhyloFrame results creates directories seperating each model's ancestry training data and the model number within that ancestry. For example within the user defined results directory, the EUR model 1 results will be in: `model_runs/phyloFrame/eur/model_1`. Within this directory test results on other ancestry samples, as well as the model and its gene signature. Associated benchmark model results will be in `model_runs/benchmark/eur/model_1`.

The results shown in this paper were run on a single core and utilized 92GB. The model runs in approximatley 2 hours and 40 minutes.  
For single runs (runs trained on a single datasets not seperated into smaller dataset batches) results will be in the user defined results directory in directories `phyloFrame/` and `benchmark/`. Wihtin the directory will be the test set results as well as model and model gene signature.

The single model was run on a single core and utilizes 87GB. The model runs in approximatley 1 hour and 40 minutes. 

#### To get results by batch
To get metrics by ancestry batch (there are multiple batches per ancestry), run the `run_predict_batches.sh` script. This script expects to run on BRCA, UCEC, and THCA results at once - you will need to have already run all three diseases and define the result directories you used for each cancer within the script. 

## Enhanced Allele Frequency Creation 
In this version of PhyloFrame Enhanced Allele Frequencies (EAFs) are calculated from population specific allele frequencies in Gnomadv4.1. In order to calculate your own EAFs, you will need vcf files from whichever population database you are using with allele frequencies for the populations you would like to include in the calculation. An example of the expected vcf format similar is the following: 
```
chr10   45366   rs1554737603    G       C       .       AC0;AS_VQSR   AC=0;AN=30;AF=0.00000;AC_fin=0;AF_fin=0.00000;AN_fin=2;AC_mid=0;AF_mid=0.00000;AN_mid=2;AC_sas=0;AF_sas=0.00000;AN_sas=2  

```
We provide a Snakemake file in the repository that may help you create your own EAF file. The Snakefile and associated scripts will help to get started but will likely need to be heavily edited to reflect things like:
  1. Where you are downloading your vcf files from - **edit Snakefile rule `download_chromosomes`**
  2. If the chromosomes are in seperate files (currently what the Snakefile expects) or if they are all in a single file - **edit Snakefile rule `download_chromosomes`**
  3. Which ancestries you are parsing from the vcf - **edit vcf_parsing/main.cpp** AND **code/gnomadV4_EAF_Calculation.R**
  4. Where you want to write output files - **edit Snakefile `workdir` AND `here_path` AND path used in rule `parse_chromosomes` AND path used in rule `calculate_EAF**

To use your calculated EAF in a PhyloFrame run please make sure it is in the `data-raw/` directory with the file name: `gnomad.exomes.all.tsv`. Or change where the file is being read from in the function `load_EAF()` in the `load_large_data.R` script. 
