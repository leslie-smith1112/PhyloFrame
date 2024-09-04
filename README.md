
# PhyloFrame
PhyloFrame uses population level information to correct for unequal ancestral representation in genomic training data. It utilizes population resources such as Gnomad, functional interaction networks like HumanBase, and transcriptomic patient data.  

## Installation

You can install the development version of PhyloFrame like so:

```{r example}
library(githubinstall)
githubinstall("PhyloFrame")

```

## Example: Replication of paper results

This is a basic example which shows you how to recreate a run from the paper as well as to run PhyloFrame on a single datatset.

### 1. Set up data
PhyloFrame utilizes two large data sources to run: Functional interaction networks and enhanced allele frequencies created from Gnomadv4.1. 

#### Functional Interaction Network
Functional interaction networks used in PhyloFrame are from Humanbase and can be downloaded here: https://hb.flatironinstitute.org/download. We use the mammary epithelium, thyroid, and uterine endometrium networks in the associated paper. We use the Full network in our analysis, however any network with the format below can be used.  In this file the columns are listed on Hummanbase as follows: [entrez gene id 1][entrez gene id 2][posterior prob., with known edges set to 1][posterior prob.] We use the [posterior prob., with known edges set to 1] connection in our analysis and drop the fourth column.

```{network file format}    
5983	1663	1	0.217193
5983	1645	1	0.0558941
5983	2255	1	0.0264535
5983	9401	1	0.513181
```
Once the network is downloaded please ensure it is unzipped and in the `data-raw` directory within the project's home directory. Then run the script `run_network_annotation.R` with the name of the network file as an argument. An example run is below:
```
devtools::load_all()
annotate_network(TISSUE_NETWORK)
```

To run a disease, download a network. 

#### Enhanced Allele Frequencies
Previously calculated enhanced allele frequencies (EAF) can be downloaded at this link: . Please keep in mind this is a large file ~28Gb. EAFs in this file were caluclated with 8 ancestries from Gnomadv4.1 exome files. 

```{r example}
library(PhyloFrame)
## basic example code
devtools::load_all()
## For reproduction of PhyloFrame results on multiple TCGA datasets
## Define disease (breast, uterine, thyroid), name of output directory (will be in results/), and whether would would like to create new training batches
devtools::load_all()
main(breast, TCGA_Breast, FALSE)

## To run PhyloFrame on a single dataset, define a list of training samples (samples not in training set are automatically used as test samples), expression matrix with prediction task column names as "subtype", the name of the directory you want your output in (will be in results/USER_DEFINED_RESULTS_DIR, and the name of the training set (this is used to name output files)
## Here we use the BRCA TCGA samples with EUR samples as training data

train_samples <- samples_ancestry$patient[samples_ancestry$consensus_ancestry == "eur"]
single_expr_driver(expression_breast,train_samples, "breast", "TCGA_BRCA", "EUR")
```


## Expected Output 
As PhyloFrame trains on many small, independent, single ancestry dataset reproduction of PhyloFrame results creates directories seperating each model's ancestry training data and the model within that ancestry. For example within the user defined results directory, the EUR model 1 results will be in: `model_runs/phyloFrame/eur/model_1`. Within this directory test results on other ancestry samples, as well as the model and its gene signature. Associated benchmark model results will be in `model_runs/benchmark/eur/model_1`.

The results shown in this paper were run on a single core and utilized 92GB. The model runs in approximatley 2 hours and 40 minutes.  

For single runs (runs trained on a single datasets not seperated into smaller dataset batches) results will be in the user defined results directory in directories `phyloFrame/` and `benchmark/`. Wihtin the directory will be the test set results as well as model and model gene signature.

The single model was run on a single core and utilizes 87GB. The model runs in approximatley 1 hour and 40 minutes. 


