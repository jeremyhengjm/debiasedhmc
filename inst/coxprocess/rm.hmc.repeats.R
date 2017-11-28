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

# no. of repetitions and mcmc iterations
nreps <- 10
k <- 71
K <- 710

# specify stepsize and no. of steps
stepsize <- 0.13
nsteps <- 10

# define hmc kernel
hmc <- get_rm_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension, metric)

# define mixture kernels
omega <- 1 / 20
Sigma_std <- 1e-3
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension)
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
runtimes <- rep(0, nreps)
meetingtime <- rep(0, nreps)
mean_estimates <- matrix(nrow = nreps, ncol = dimension)
var_estimates <- matrix(nrow = nreps, ncol = dimension)
for(irep in 1:nreps){
  tic()
  cchains <- coupled_chains(mixture_kernel, mixture_coupled_kernel, rinit, K = K)
  timing <- toc()
  runtime <- timing$toc - timing$tic
  runtimes[irep] <- runtime
  meetingtime[irep] <- cchains$meetingtime
  mean_estimates[irep, ] <- H_bar(cchains, h = function(x) x, k, K)
  var_estimates[irep, ] <- H_bar(cchains, h = function(x) x^2, k, K)
  cat("Repetition:", irep, "/", nreps, "\n")
}

filename <- paste("output.rm.hmc.repeat", igrid, ".RData", sep = "")
save(nreps, k, K, stepsize, nsteps, runtimes, meetingtime,
     mean_estimates, var_estimates, file = filename, safe = F)
