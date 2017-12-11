rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(coda)

# load model
load("germancredit.RData")

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
grid_stepsize <- seq(0.01, 0.1, by = 0.01)
ngrid_stepsize <- length(grid_stepsize)
nsteps <- 1 # 0.0185 seconds per mcmc iteration

# pre-allocate
variance <- rep(0, ngrid_stepsize)
acceptprob <- rep(0, ngrid_stepsize)
runtimes <- rep(0, ngrid_stepsize)
filename <- paste("output.mala.germancredit.tuning.RData")

for (istepsize in 1:ngrid_stepsize){
  # define hmc kernel
  stepsize <- grid_stepsize[istepsize]
  hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

  # run hmc
  chain <- matrix(0, nrow = nmcmc, ncol = dimension)
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
    cat("No. of steps", stepsize, "Iteration:", imcmc, "/", nmcmc, "\n")
  }
  timing <- toc()
  runtime <- timing$toc - timing$tic
  variance[istepsize] <- compute_variance(chain[burnin:nmcmc, ])
  acceptprob[istepsize] <- accept / nmcmc
  runtimes[istepsize] <- runtime
  save(grid_stepsize, ngrid_stepsize,
       variance, acceptprob, runtimes, file = filename, safe = F)

}

