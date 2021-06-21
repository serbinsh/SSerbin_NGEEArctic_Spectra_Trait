#!/bin/bash
#######################################################################################
#SBATCH --time=68:00:00
#SBATCH --job-name=prospect-slurm
#SBATCH --out="prospect_job-%j.out"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=6GB
#SBATCH --mail-type=ALL
#######################################################################################

#######################################################################################
### SLURM usage example
# sbatch -w node03 -c 3 --partition compute --job-name=prospect --mail-user=sserbin@bnl.gov batch/invert_prospect_slurm.sh
#######################################################################################

script_name=R_scripts/PROSPECT/Invert_leaf_refl_spectra_PROSPECT-D.R

echo "Running "${script_name}
echo "Starting at: " `date`
echo "Job submitted to the ${SLURM_JOB_PARTITION} partition on ${SLURM_CLUSTER_NAME}"
echo "Job name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}"
echo "There are ${SLURM_CPUS_ON_NODE} CPUs on compute node $(hostname)"

echo "Run started on: " `date`

Rscript --vanilla ${script_name}
