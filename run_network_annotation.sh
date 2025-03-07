#!/bin/bash
#SBATCH --job-name=network_annotation    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email>     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --time=120:05:00              # Time limit hrs:min:sec
#SBATCH --output=network_annotation%j.stdout   # Standard output and error log
#SBATCH --error=network_annotation%j.err # Error log
ml R

#define links from humanbase
hb_links=("https://s3-us-west-2.amazonaws.com/humanbase/networks/mammary_epithelium.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/thyroid_gland.gz" "https://s3-us-west-2.amazonaws.com/humanbase/networks/uterine_endometrium.gz")
# set full path to data raw from working directory
prefix=$(pwd)/data-raw/
for i in ${hb_links[@]}
do
	#get the tissue network file - raw network file put in data-raw/
	wget ${i} -P ${prefix}
	#extract the tissue network name from the file
	tissue_gz=( $(echo ${i} | awk '{split($0,chopped,"/"); print chopped[6]}') )
	#unzip the file so we can work on it in R
	gzip -d ${prefix}${tissue_gz}
	#get the name of the file without .gz ending
	tissue=( $(echo ${tissue_gz} | sed 's/.gz//g') )
	# run_network_annotation.R calls function from R/network_annotation.R
	#annotate the network in R - output in data-raw/
	Rscript ./run_network_annotation.R ${tissue}

done



