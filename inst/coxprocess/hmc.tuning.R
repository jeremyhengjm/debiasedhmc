rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)
library(coda)

# load cox process model
load("coxprocess.RData")

# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

# compute MCMC variance
compute_variance <- function(trajectory){
  # variance of means
  total_var <- sum(spectrum0.ar(trajectory)$spec)
  # variance of second moments
  total_var <- total_var + sum(spectrum0.ar(trajectory^2)$spec)

  # ess of covariances (too expensive)
  # sample_mean <- colMeans(trajectory)
  # for (i in 1:dimension){
  #   for (j in 1:dimension){
  #     total_ess <- total_ess + effectiveSize( (trajectory[, i] - sample_mean[i]) *
  #                                             (trajectory[, j] - sample_mean[j]) )
  #   }
  # }
  return(total_var)
}

# no. of mcmc iterations
nmcmc <- 11000
burnin <- 1001

# vary stepsize and no. of steps
grid_stepsize <- seq(0.05, 0.45, by = 0.02)
ngrid_stepsize <- length(grid_stepsize)
stepsize <- grid_stepsize[igrid]
grid_nsteps <- c(1, 10, 20, 30) # c(0.4, 1.2, 2.4, 3.6) seconds per mcmc iteration
ngrid_nsteps <- length(grid_nsteps)

# pre-allocate
variance <- rep(0, ngrid_nsteps)
acceptprob <- rep(0, ngrid_nsteps)
runtimes <- rep(0, ngrid_nsteps)
filename <- paste("output.hmc.tuning", igrid, ".RData", sep = "")

for (istep in 1:ngrid_nsteps){
  # define hmc kernel
  nsteps <- grid_nsteps[istep]
  hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

  # run hmc
  chain <- matrix(nrow = nmcmc, ncol = dimension)
  current_x <- rinit() # initialize
  current_pdf <- logtarget(current_x)
  accept <- 0
  tic()
  for (imcmc in 1:nmcmc){
    current_state <- hmc$kernel(current_x, current_pdf, imcmc)
    current_x <- current_state$chain_state
    current_pdf <- current_state$current_pdf
    chain[imcmc, ] <- current_x
    accept <- accept + current_state$accept
    cat("No. of steps", nsteps, "Iteration:", imcmc, "/", nmcmc, "\n")
  }
  timing <- toc()
  runtime <- timing$toc - timing$tic
  variance[istep] <- compute_variance(chain[burnin:nmcmc, ])
  acceptprob[istep] <- accept / nmcmc
  runtimes[istep] <- runtime
  save(grid_stepsize, ngrid_stepsize, grid_nsteps, ngrid_nsteps,
       variance, acceptprob, runtimes, file = filename, safe = F)

}

