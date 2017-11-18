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
adaptive_tol <- 0.1 # tune this?

# define mala kernel
mala <- get_mala_kernel(logtarget, gradlogtarget, stepsize, dimension, adaptive_tol)

# compute meeting times
distance <- list()
meetingtime <- rep(0, nreps)
for (irep in 1:nreps){
  cat("Repetition:", irep, "/", nreps, "\n")
  tic()
  cchains <- coupled_chains(mala$kernel, mala$coupled_kernel, rinit, K = K, max_iterations = max_iterations)
  toc()
  distance[[irep]] <- compute_distance(cchains)
  meetingtime[irep] <- cchains$meetingtime
}

filename <- paste("output.mala.meetingtimes", igrid, ".RData", sep = "")
save(nreps, max_iterations, stepsize, adaptive_tol,
     distance, meetingtime, file = filename, safe = F)




