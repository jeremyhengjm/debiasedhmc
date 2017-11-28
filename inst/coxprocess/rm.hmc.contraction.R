rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)

# load cox process model
load("coxprocess_with_metric.RData")

# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

# compute distance
compute_distance <- function(cchain){
  niterations <- nrow(cchain$samples2)
  return(sqrt(sum((cchain$samples1[1+niterations, ] - cchain$samples2[niterations, ])^2)))
}

# no. of repetitions and mcmc iterations
nreps <- 5
max_iterations <- 1000

# vary stepsize and no. of steps
grid_stepsize <- seq(0.05, 0.45, by = 0.02)
ngrid_stepsize <- length(grid_stepsize)
stepsize <- grid_stepsize[igrid]
grid_nsteps <- c(1, 10, 20, 30) # c(0.4, 1.2, 2.4, 3.6) seconds per mcmc iteration
ngrid_nsteps <- length(grid_nsteps)

# pre-allocate
distance <- matrix(nrow = nreps, ncol = ngrid_nsteps)
filename <- paste("output.rm.hmc.contraction", igrid, ".RData", sep = "")

for (istep in 1:ngrid_nsteps){
  # define hmc kernel
  nsteps <- grid_nsteps[istep]
  hmc <- get_rm_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension, metric)

  # run hmc
  for (irep in 1:nreps){
    cchains <- coupled_chains(hmc$kernel, hmc$coupled_kernel, rinit, max_iterations = max_iterations)
    distance[irep, istep] <- compute_distance(cchains)
    cat("No. of steps", nsteps, "Repetition:", irep, "/", nreps, "\n")
    save(grid_stepsize, ngrid_stepsize, grid_nsteps, ngrid_nsteps,
         distance, file = filename, safe = F)
  }

}





