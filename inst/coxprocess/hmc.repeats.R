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

# no. of repetitions and mcmc iterations
nreps <- 10
K <- 1000

# specify stepsize and no. of steps
stepsize <- 0.1
nsteps <- 10

# define hmc kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# define mixture kernels
omega <- 1 / 20
Sigma_proposal <- 1e-10 * diag(1, dimension, dimension)
mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension)

# mixture kernels
mixture_kernel <- function(chain_state, iteration){
  if (runif(1) < omega){
    return(mh$kernel(chain_state, iteration))
  } else {
    return(hmc$kernel(chain_state, iteration))
  }
}

mixture_coupled_kernel <- function(chain_state1, chain_state2, iteration){
  if (runif(1) < omega){
    return(mh$coupled_kernel(chain_state1, chain_state2, iteration))
  } else {
    return(hmc$coupled_kernel(chain_state1, chain_state2, iteration))
  }
}

# compute estimates
mean_estimates <- matrix(nrow = nreps, ncol = dimension)
var_estimates <- matrix(nrow = nreps, ncol = dimension)
for(irep in 1:nreps){
  cchains <- coupled_chains(mixture_kernel, mixture_coupled_kernel, rinit, K = K)
  mean_estimates[irep, ] <- H_bar(cchains, h = function(x) x, k, K)
  var_estimates[irep, ] <- H_bar(cchains, h = function(x) x^2, k, K)
  cat("Repetition:", irep, "/", nreps, "\n")
}

filename <- paste("output.hmc.meetingtimes", igrid, ".RData", sep = "")
save(nreps, K, stepsize, nsteps,
     mean_estimates, var_estimates, file = filename, safe = F)
