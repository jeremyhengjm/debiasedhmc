#!/bin/bash
#SBATCH -J hmc_tuning                           # Job name
#SBATCH -n 1                                    # Number of cores
#SBATCH -N 1                                    # All cores on one machine
#SBATCH -t 0-02:30                              # Runtime in D-HH:MM
#SBATCH -p shared                               # Partition to submit to
#SBATCH --mem-per-cpu=500M                    # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END                         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jjmheng@fas.harvard.edu     # Email
#SBATCH --array=1-13

## LOAD SOFTWARE ENV ##
source new-modules.sh
module load R/3.4.2-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
input=hmc.tuning.R
cd /n/home07/jeremyhengjm/apps/R/debiasedhmc/inst/logistic

srun R CMD BATCH --no-save $input out/$input.$SLURM_ARRAY_TASK_ID.out


