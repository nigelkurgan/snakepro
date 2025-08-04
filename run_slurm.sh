#!/bin/bash

### slurm job options

#SBATCH --job-name=quant_%j    # job name
#SBATCH --mail-type=END,FAIL    # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jwn351@ku.dk    # email address to receive the notification    
#SBATCH -c 1    # number of requested cores
#SBATCH --mem=2gb    # total requested RAM
#SBATCH --time=1-00:00:00               # max. running time of the job, format in D-HH:MM:SS
#SBATCH --output=logs/quant_%j.log  # standard output and error log, '%j' gives the job ID 

if [[ $# -eq 0 ]] ; then
    target='all'
else
    target="$@"
fi

echo "Targets: $target"


module load miniconda/24.5.0
eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake                                                            \
    --slurm                                                          \
    --default-resources slurm_partition=standardqueue                \
    --rerun-triggers mtime                                           \
    --jobs 32                                                        \
    --cluster-config "config/slurm.yml"                              \
    --configfile config/config.yml                                   \
    --keep-going      \
    --latency-wait 60                                                \
    $target


