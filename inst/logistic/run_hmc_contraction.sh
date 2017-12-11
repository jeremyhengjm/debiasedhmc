#!/bin/bash
#SBATCH -J contract66                             # Job name
#SBATCH -n 1                                    # Number of cores
#SBATCH -N 1                                    # All cores on one machine
#SBATCH -t 0-01:00                              # Runtime in D-HH:MM
#SBATCH -p shared                               # Partition to submit to
#SBATCH --mem-per-cpu=500M                    # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END                         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jjmheng@fas.harvard.edu   # Email
#SBATCH --array=1-10

## LOAD SOFTWARE ENV ##
source new-modules.sh
module load R/3.4.2-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
input=hmc.contraction.R
cd /n/home07/jeremyhengjm/apps/R/debiasedhmc/inst/logistic

srun R CMD BATCH "--args 1" --no-save $input out/$input.$SLURM_ARRAY_TASK_ID.out

# dimension = 66 requires 0.4 hours
# dimension = 130 requires 0.7 hours
# dimension = 258 requires 1.6 hours
# dimension = 514 requires 3.2 hours
# dimension = 1026 requires 7 hours
