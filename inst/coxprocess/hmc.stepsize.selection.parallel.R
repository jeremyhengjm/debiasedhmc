rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)

# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

# load cox process model
load("coxprocess.RData")

# Hamiltonian function
hamiltonian <- function(x, v) -logtarget(x) + sum(v^2) / 2

# stepsize selection
nsteps <- 100
grid_stepsize <- 1e-2
# grid_stepsize <- seq(1e-2, 1e-1, length.out = 10)
nstepsizes <- length(grid_stepsize)
nreps <- 2
hamiltonian_error <- rep(0, nreps)

tic("Runtime:")
  stepsize <- grid_stepsize[igrid]
  for (irep in 1:nreps){
    cat("Repetition:", irep, "/", nreps, "\n")
    # store hamiltonian error along trajectory
    h_error <- rep(0, nsteps)

    # initialize
    current_x <- rinit()
    current_v <- rnorm(dimension)
    initial_hamiltonian <- hamiltonian(current_x, current_v)

    # perform leapfrog integration
    for (istep in 1:nsteps){
      current_v <- current_v + stepsize * gradlogtarget(current_x) / 2
      current_x <- current_x + stepsize * current_v
      current_v <- current_v + stepsize * gradlogtarget(current_x) / 2
      h_error[istep] <- abs( hamiltonian(current_x, current_v) - initial_hamiltonian )
    }

    # store maximum hamiltonian error along trajectory
    hamiltonian_error[irep] <- max(h_error)
  }


toc()
filename <- paste("output.hmc.stepsize.selection.", igrid, ".RData", sep = "")

# save(nstepsizes, grid_stepsize, nreps, hamiltonian_error, file = filename, safe = F)





