#!/bin/bash
#SBATCH -J para1                                # Job name
#SBATCH -n 1                                    # Number of cores
#SBATCH -N 1                                    # All cores on one machine
#SBATCH -t 0-02:30                              # Runtime in D-HH:MM (nsteps = 1, 10, 20, 30 takes 01:30, 02:30, 03:30, 05:00)
#SBATCH -p shared                               # Partition to submit to
#SBATCH --mem-per-cpu=3000M                     # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END                         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jjmheng@fas.harvard.edu     # Email
#SBATCH --array=1-50

## LOAD SOFTWARE ENV ##
source new-modules.sh
module load R/3.4.2-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
input=rm.hmc.meetingtimes.R
cd /n/home07/jeremyhengjm/apps/R/debiasedhmc/inst/coxprocess

srun R CMD BATCH "--args 1" --no-save $input out/$input.$SLURM_ARRAY_TASK_ID.out

#sleep 1                                          # pause to be kind to the scheduler
# parameter: 1  2  3  4  5  6  7  8  9 10
# nsteps:   10 20 30 10 20 30 10 20 10 20
# parameter: 11 12 13 14 15 16 17 18 19 20 21 22 23
# nsteps:    10 20 10 10 10 10 10 10  1  1  1  1  1
