#!/bin/bash
#SBATCH --job-name=EAF_Creation
#SBATCH --account=kgraim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leslie.smith1@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=290:00:00
#SBATCH --output=logs/%j_eaf.stdout
#SBATCH --error=logs/%j_eaf.stderr

pwd; hostname; date
module load snakemake/7.32.4

snakemake --cluster "sbatch -A {cluster.account} -q {cluster.qos}  -c {cluster.cpus-per-task} -N {cluster.Nodes} -t {cluster.runtime} --mem {cluster.mem} -J {cluster.jobname} \
 --mail-user={cluster.mail} --output {cluster.out} --error {cluster.err}" --cluster-config config/cluster_config.json --jobs 5 --latency-wait 20 --rerun-incomplete --use-envmodules --keep-going
