#!/ibin/bash

if [[ $# -eq 0 ]] ; then
    target='all'
else
    target="$@"
fi

echo "Targets: $target"


module load miniconda/24.5.0
module load snakemake/7.30.1
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
