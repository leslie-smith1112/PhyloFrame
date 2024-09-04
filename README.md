
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

# Set up data


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


