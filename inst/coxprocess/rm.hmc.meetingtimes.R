rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)

# load cox process model
load("coxprocess_with_metric.RData")

# specify stepsize and no. of steps
load("rm.hmc.contraction.parameters.RData")
args <- commandArgs(TRUE)
iparameter <- as.integer(args[1])
cat("iparameter:", iparameter, "\n")
stepsize <- contraction.df$stepsize[iparameter]
nsteps <- contraction.df$nsteps[iparameter] # L = 10, 20, 30 takes 1.2, 2.4, 3.6 secs per HMC iteration
# cat("Expected runtime:", 2.4 * 2 * max_iterations * nreps / 3600, "hours \n")

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
  return(sapply(1:niterations, function(index) sqrt(sum((cchain$samples1[1+index, ] - cchain$samples2[index, ])^2))))
}

# no. of repetitions and mcmc iterations
nreps <- 2
m <- 1000
max_iterations <- 1000 # L = 1, 10, 20, 30 takes 1.3, 3, 5, 7 seconds per iteration for coupled chain

# define rm-hmc kernel
hmc <- get_rm_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension, metric)

# define mixture kernels
omega <- 1 / 20
Sigma_std <- 1e-3
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension) # try 1e-8, 1e-10 also
mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension)

# mixture kernels
mixture_kernel <- function(chain_state, current_pdf, iteration){
  if (runif(1) < omega){
    return(mh$kernel(chain_state, current_pdf, iteration))
  } else {
    return(hmc$kernel(chain_state, current_pdf, iteration))
  }
}

mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration){
  if (runif(1) < omega){
    return(mh$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration))
  } else {
    return(hmc$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration))
  }
}

distance <- list()
meetingtime <- rep(0, nreps)
for (irep in 1:nreps){
  cat("Repetition:", irep, "/", nreps, "\n")
  tic()
  cchains <- coupled_chains(logtarget, mixture_kernel, mixture_coupled_kernel, rinit, m = m, max_iterations = max_iterations)
  toc()
  distance[[irep]] <- compute_distance(cchains)
  meetingtime[irep] <- cchains$meetingtime
}

filename <- paste("output.rm.hmc.parameter", iparameter, ".rep", igrid, ".RData", sep = "")
save(nreps, max_iterations, stepsize, nsteps, Sigma_std,
     distance, meetingtime, file = filename, safe = F)
