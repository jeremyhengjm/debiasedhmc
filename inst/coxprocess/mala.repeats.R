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
adaptive_tol <- 0.1 # tune this?

# define mala kernel
mala <- get_mala_kernel(logtarget, gradlogtarget, stepsize, dimension, adaptive_tol)

# compute estimates
mean_estimates <- matrix(nrow = nreps, ncol = dimension)
var_estimates <- matrix(nrow = nreps, ncol = dimension)
for(irep in 1:nreps){
  cchains <- coupled_chains(mala$kernel, mala$coupled_kernel, rinit, K = K)
  mean_estimates[irep, ] <- H_bar(cchains, h = function(x) x, k, K)
  var_estimates[irep, ] <- H_bar(cchains, h = function(x) x^2, k, K)
  cat("Repetition:", irep, "/", nreps, "\n")
}

filename <- paste("output.mala.meetingtimes", igrid, ".RData", sep = "")
save(nreps, K, stepsize, adaptive_tol,
     mean_estimates, var_estimates, file = filename, safe = F)
