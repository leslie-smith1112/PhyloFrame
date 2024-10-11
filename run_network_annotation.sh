#!/bin/bash
#SBATCH --job-name=network_annotation    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=leslie.smith1@ufl.edu     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150gb
#SBATCH --time=120:05:00              # Time limit hrs:min:sec
#SBATCH --output=network_annotation%j.stdout   # Standard output and error log
#SBATCH --error=network_annotation%j.err # Error log
ml R

#define links from humanbase
hb_links=("https://s3-us-west-2.amazonaws.com/humanbase/networks/blood.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/adrenal_gland.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/urinary_bladder.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/brain.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/uterine_cervix.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/global.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/bone_marrow.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/colon.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/esophagus.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/kidney.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/liver.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/lung.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/lymph_node.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/pancreas.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/peripheral_nervous_system.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/prostate_gland.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/gastrointestinal_tract.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/skin.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/stomach.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/testis.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/uterus.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/eye.gz")

prefix=/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/PhyloFrame/data-raw/
for i in ${hb_links[@]}
do
	#get the tissue network file
	wget ${i} -P ${prefix}
	#extract the tissue network name from the file
	tissue_gz=( $(echo ${i} | awk '{split($0,chopped,"/"); print chopped[6]}') )
	#unzip the file so we can work on it in R
	gzip -d ${prefix}${tissue_gz}
	#get the name of the file without .gz ending
	tissue=( $(echo ${tissue_gz} | sed 's/.gz//g') )
	#annotate the network in R - output in ./data-raw/network_downloads_annotated
	Rscript ./run_network_annotation.R ${tissue}

done



