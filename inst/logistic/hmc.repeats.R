rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)

# load model
load("germancredit.RData")

# initialize with crude normal approximation
load("germancredit.normalapprox.RData")
rinit <- function() as.numeric(fast_rmvnorm(1, approx_mean, diag(approx_var, dimension, dimension)))

# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

# no. of repetitions and mcmc iterations
nreps <- 10
k <- 330
m <- 10 * k

# specify stepsize and no. of steps
stepsize <- 0.0125
nsteps <- 10

# define hmc kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# define mixture kernels
omega <- 1 / 20
Sigma_std <- 1e-3
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension)
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

# compute estimates
runtimes <- rep(0, nreps)
meetingtime <- rep(0, nreps)
unbiased_estimates <- list()
# unbiased_estimates <- matrix(nrow = nreps, ncol = 2 * dimension) # first and second moment
# mcmc_estimates <- matrix(nrow = nreps, ncol = 2 * dimension)

for(irep in 1:nreps){
  tic()
  # estimation_output <- unbiased_estimator(logtarget, mixture_kernel, mixture_coupled_kernel, rinit,
  # h = function(x) c(x, x^2), k = k, m = m)
  cchains <- coupled_chains(logtarget, mixture_kernel, mixture_coupled_kernel, rinit, m = m)
  timing <- toc()
  runtime <- timing$toc - timing$tic
  runtimes[irep] <- runtime

  meetingtime[irep] <- cchains$meetingtime
  H_bar_estimates <- matrix(nrow = 2 * dimension, ncol = m)
  for(t in 1:m){
    H_bar_estimates[, t] <- H_bar(cchains, h = function(x) c(x, x^2), k = t, m = t)
  }
  unbiased_estimates[[irep]] <- H_bar_estimates

  # meetingtime[irep] <- estimation_output$meetingtime
  # unbiased_estimates[irep, ] <- estimation_output$uestimator
  # mcmc_estimates[irep, ] <- estimation_output$mcmcestimator
  cat("Repetition:", irep, "/", nreps, "\n")
}

filename <- paste("output.hmc.germancredit.repeat", igrid, ".RData", sep = "")
save(nreps, k, m, stepsize, nsteps, runtimes, meetingtime,
     unbiased_estimates, file = filename, safe = F)
     # unbiased_estimates, mcmc_estimates, file = filename, safe = F)


