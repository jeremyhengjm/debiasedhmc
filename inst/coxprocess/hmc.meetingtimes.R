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

# compute distance
compute_distance <- function(cchain){
  niterations <- nrow(cchain$samples2)
  return(sapply(1:niterations, function(index) sqrt(sum((cchain$samples1[1+index, ] - cchain$samples2[index, ])^2))))
}

# no. of repetitions and mcmc iterations
nreps <- 5
K <- 1000
max_iterations <- 1000

# specify stepsize and no. of steps
stepsize <- 0.25
nsteps <- 10 # 2.2 sec per iteration of coupled chain

# define hmc kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# define mixture kernels
omega <- 1 / 20
Sigma_proposal <- 1e-10 * diag(1, dimension, dimension)
tic()
mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension)
toc()
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

distance <- list()
meetingtime <- rep(0, nreps)
for (irep in 1:nreps){
  cat("Repetition:", irep, "/", nreps, "\n")
  tic()
  cchains <- coupled_chains(mixture_kernel, mixture_coupled_kernel, rinit, K = K, max_iterations = max_iterations)
  toc()
  distance[[irep]] <- compute_distance(cchains)
  meetingtime[irep] <- cchains$meetingtime
}

filename <- paste("output.hmc.meetingtimes", igrid, ".RData", sep = "")
save(nreps, max_iterations, stepsize, nsteps,
     distance, meetingtime, file = filename, safe = F)




