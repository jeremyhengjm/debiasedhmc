rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)

# load model
load("germancredit.RData")

# specify stepsize and no. of steps
stepsize <- 0.0125
nsteps <- 10
# load("hmc.germancredit.contraction.parameters.RData")
# args <- commandArgs(TRUE)
# iparameter <- as.integer(args[1])
# cat("iparameter:", iparameter, "\n")
# stepsize <- contraction.df$stepsize[iparameter]
# nsteps <- contraction.df$nsteps[iparameter] # L = 10, 20, 30 takes c(0.0185, 0.0865, 0.16, 0.23) seconds per mcmc iteration


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
  return(sapply(1:niterations, function(index) sqrt(sum((cchain$samples1[1+index, ] - cchain$samples2[index, ])^2))))
}

# no. of repetitions and mcmc iterations
nreps <- 10
m <- 1000
max_iterations <- 1000

# define hmc kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# define mixture kernels
omega <- 1 / 20
Sigma_std <- 1e-7
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension) # try 1e-3, 1e-5, 1e-7
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

filename <- paste("output.hmc.germancredit.sigma1e-7.", igrid, ".RData", sep = "")
save(nreps, max_iterations, stepsize, nsteps, Sigma_std,
     distance, meetingtime, file = filename, safe = F)
