rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)

# load model
# load("germancredit.RData")
args <- commandArgs(TRUE)
idimension <- as.integer(args[1])
filename <- switch(idimension, "simulated_nsamples1000_dimension66.RData",
                               "simulated_nsamples1000_dimension130.RData",
                               "simulated_nsamples1000_dimension258.RData",
                               "simulated_nsamples1000_dimension514.RData",
                               "simulated_nsamples1000_dimension1026.RData")
load(file = filename)
cat("Dimension =", dimension, "\n")

# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

# compute distance
compute_distance <- function(cchain){
  niterations <- dim(cchain$samples2)[1]
  return(sqrt(sum((cchain$samples1[1+niterations, ] - cchain$samples2[niterations, ])^2)))
}

# no. of repetitions and mcmc iterations
nreps <- 5
max_iterations <- 1000

# vary stepsize and no. of steps
# grid_stepsize <- seq(0.01, 0.04, by = 0.0025) # german credit dataset
grid_stepsize <- switch(idimension, seq(0.01, 0.12, length.out = 10),
                                    seq(0.01, 0.10, length.out = 10),
                                    seq(0.01, 0.10, length.out = 10),
                                    seq(0.01, 0.06, length.out = 10),
                                    seq(0.01, 0.06, length.out = 10) ) # simulated dataset
ngrid_stepsize <- length(grid_stepsize)
stepsize <- grid_stepsize[igrid]
grid_nsteps <- c(1, 10, 20, 30)
# c(0.0185, 0.0865, 0.16, 0.23) seconds per mcmc iteration for germancredit dataset
# c(0.005, 0.025, 0.04, 0.06) seconds per mcmc iteration for simulated dataset with d = 66
# c(0.01, 0.045, 0.075, 0.11) seconds per mcmc iteration for simulated dataset with d = 130
# c(0.02, 0.10, 0.18, 0.26) seconds per mcmc iteration for simulated dataset with d = 258
# c(0.048, 0.20, 0.40, 0.50) seconds per mcmc iteration for simulated dataset with d = 514
# c(0.09, 0.40, 0.80, 1.2) seconds per mcmc iteration for simulated dataset with d = 1026
ngrid_nsteps <- length(grid_nsteps)

# pre-allocate
distance <- matrix(nrow = nreps, ncol = ngrid_nsteps)

for (istep in 1:ngrid_nsteps){
  # define hmc kernel
  nsteps <- grid_nsteps[istep]
  hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

  # run hmc
  for (irep in 1:nreps){
    cchains <- coupled_chains(logtarget, hmc$kernel, hmc$coupled_kernel, rinit, max_iterations = max_iterations)
    distance[irep, istep] <- compute_distance(cchains)
    cat("No. of steps", nsteps, "Repetition:", irep, "/", nreps, "\n")
    filename <- paste("output.hmc.dimension", idimension, ".contraction", igrid, ".RData", sep = "")
    save(grid_stepsize, ngrid_stepsize, grid_nsteps, ngrid_nsteps,
         distance, file = filename, safe = F)
  }

}






