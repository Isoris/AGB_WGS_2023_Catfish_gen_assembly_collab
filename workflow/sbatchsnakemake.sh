#!/bin/bash
#SBATCH -p compute
#SBATCH -n 1 -c 24
#SBATCH -t 24:00:00
#SBATCH -J Snakemake_quentin
#SBATCH -A proj5034

# Initialize Mamba/Conda
eval "$(conda shell.bash hook)"
source /tarafs/scratch/proj5034-AGBKU/catfish/01-PROGRAMS/miniforge3/etc/profile.d/mamba.sh

# Activate the environment
mamba activate snakemake

# Run snakemake 

snakemake --cluster "sbatch -A CLUSTER_ACCOUNT -t CLUSTER_TIME -p CLUSTER_PARTITION -N CLUSTER_NODES" --jobs NUM_JOBS_TO_SUBMIT

conda deactivate
