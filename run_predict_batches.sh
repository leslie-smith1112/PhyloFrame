#!/bin/bash
#SBATCH --job-name=predict_batches    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email>     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=120:05:00              # Time limit hrs:min:sec
#SBATCH --output=predict_batches_%j.stdout   # Standard output and error log
#SBATCH --error=predict_batches_%j.err # Error log

ml R

echo -e "\e[1;32mGetting list of batch associations for each sample. \e[0m"
# Get the batches for all samples - needed for prep_batches()
cd data-raw #needs to run in the directory where disease sample are
pwd
sample_dir=("brca" "ucec" "thca")
ancestries=("admix" "afr" "eas" "eur" "mixed")
for disease in ${sample_dir[@]}; do #for loop for disaese directories
  dz_samples="${disease}_samples" #directory that holds disease samples
  for ancestry in ${ancestries[@]}; do
      sample_files=( $(ls "${dz_samples}/${ancestry}") )
      #echo ${sample_files[@]}
          for sample in ${sample_files[@]}; do
              name="${sample%.*}" # get rid of .tsv ending
              toadd=( $(echo ${name} | sed "s/samples/batch/")) # replace samples with batch to get: ex: batch3
              #add ancestry and batch # to a new sample file
              sed "s/$/\t${ancestry}\t${toadd}/" ${dz_samples}/${ancestry}/${sample} > ${dz_samples}/${ancestry}/${name}_batch.tsv
          done
  done
  #append disease as a 4th column
  cat *_samples/*/*_batch.tsv | awk -v var="$disease" '{print $0, "\t" var}' >> all_disease_sample_batches.tsv
done


cd .. # go back up to main directory

echo -e "\e[1;32mFinding predicted subtype for each sample and matching to batch.\e[0m"
# Get predictions for all samples - needed for prep_batches()
# Set Driectories

brca_dir="TCGA_Breast_Gnomad4_tester"
echo -e "\e[31mUsing defined breast directory: ${brca_dir}\e[0m"
ucec_dir="TCGA_Uterine_Gnomad4_tester"
echo -e "\e[31mUsing defined uterine directory: ${ucec_dir}\e[0m"
thca_dir="TCGA_Thyroid_Gnomad4_tester"
echo -e "\e[31mUsing defined thyroid directory: ${thca_dir}\e[0m"
TCGA=TRUE #keep true by default

#whether this directory is a TCGA run directory (result file names are different in validation set)
#TCGA=$2
#echo ${master_dir}
#for both phyloframe and benchmark
    #for every ancestry model
        #get list of all models for ancestry
        #append ancestry model, model number, and version to results file
        #write to temp edited file
    #cat all edited files to a main results file containing all ancestries and models
#benchmark and phyloFrame models will be matched together and with metadata in R

versions=("phyloFrame" "benchmark")
ancestries=("admix" "afr" "eas" "eur" "mixed")
disease_dir=("results/$brca_dir/model_runs" "results/$ucec_dir/model_runs" "results/$thca_dir/model_runs")

for master_dir in ${disease_dir[@]}; do
  echo -e "sample_id\tsubtype\t.pred_class\t.pred_subtype1\t.pred_subtype2\ttrain_data\tmodel_num\ttest_data\tversion" > ${master_dir}/final_results.tsv

  for version in ${versions[@]}; do
    for ancestry in ${ancestries[@]}; do
      models=( $(ls ${master_dir}/${version}/${ancestry} | grep model) ) #get list of models
      #echo ${models[@]}
        for model in ${models[@]}; do
          # For the validation set - where we dont have multiple batches
          if [[ ${TCGA} == FALSE ]]; then # results in all ancestry batch tests for the multi are the same
            sed "s/$/\t${ancestry}\t${model}\t${version}/" ${master_dir}/${version}/${ancestry}/${model}/multi_Multi_study_validation_results.tsv >> ${master_dir}/final_results.tsv
          else
            mod_num=( $(echo ${model} | awk '{split($0,a,"_"); print a[2]}' ) ) #get model number
          for anc_test in ${ancestries[@]}; do
            file_name=${anc_test}_${mod_num}_results.tsv #get the sample predictions for this ancestry -
            #echo "12|23|11" | awk '{split($0,a,"|"); print a[3],a[2],a[1]}'
            if test -f ${master_dir}/${version}/${ancestry}/${model}/${file_name}; then
              sed "s/$/\t${ancestry}\t${model}\t${anc_test}_test\t${version}/" ${master_dir}/${version}/${ancestry}/${model}/${file_name} >> ${master_dir}/final_results.tsv
            fi
          done
          fi

        done
    done
  done
done

echo -e "\e[1;32mCalculating predictions metrics by batch.\e[0m"
Rscript ./run_predict_batches.R $brca_dir $ucec_dir $thca_dir


