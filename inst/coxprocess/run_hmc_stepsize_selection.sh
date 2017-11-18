#!/bin/bash
#SBATCH -n 1                                    # Number of cores
#SBATCH -N 1                                     # Ensure that all cores are on one machine
#SBATCH -t 0-100:00                              # Runtime in D-HH:MM
#SBATCH -p serial_requeue                                 # Partition to submit to
#SBATCH --mem-per-cpu=10000M                     # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END                          # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jjmheng@fas.harvard.edu   # Email

## LOAD SOFTWARE ENV ##
source new-modules.sh
module load R/3.4.2-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
input=hmc.stepsize.selection.R
# cd /n/home07/jeremyhengjm/apps/R/debiasedhmc/inst/coxprocess

R CMD BATCH $input out/$input.out
