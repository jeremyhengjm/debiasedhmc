rm(list = ls())
library(tictoc)

# load pine saplings dataset
library(spatstat)
data(finpines)
data_x <- (finpines$x + 5) / 10 # normalize data to unit square
data_y <- (finpines$y + 8) / 10
plot(x = data_x, y = data_y, type = "p")

ngrid <- 64
grid <- seq(from = 0, to = 1, length.out = ngrid+1)
dimension <- ngrid^2
data_counts <- rep(0, dimension)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    logical_y <- (data_x > grid[i]) * (data_x < grid[i+1])
    logical_x <- (data_y > grid[j]) * (data_y < grid[j+1])
    data_counts[(i-1)*ngrid + j] <- sum(logical_y * logical_x)
  }
}

# prior distribution
parameter_sigmasq <- 1.91
parameter_mu <- log(126) - 0.5 * parameter_sigmasq
parameter_beta <- 1 / 33
parameter_area <- 1 / dimension

prior_mean <- rep(parameter_mu, dimension)
prior_cov <- matrix(nrow = dimension, ncol = dimension)
for (m in 1:dimension){
  for (n in 1:dimension){
    index_m <- c( floor((m-1) / ngrid) + 1, ((m-1) %% ngrid) + 1 )
    index_n <- c( floor((n-1) / ngrid) + 1, ((n-1) %% ngrid) + 1 )
    prior_cov[m,n] <- parameter_sigmasq * exp(- sqrt(sum((index_m - index_n)^2)) / (ngrid * parameter_beta) )
  }
}
prior_precision <- solve(prior_cov)
prior_precision_chol <- t(chol(prior_precision))
prior <- list()
prior$logdensity <- function(x){
  return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), prior_mean, prior_precision_chol))
}
prior$gradlogdensity <- function(x){
  return(gradlognormal(x, prior_mean, prior_precision))
}

# likelihood function
likelihood <- list()
likelihood$logdensity <- function(x){
  return(coxprocess_loglikelihood(matrix(x, nrow=1), data_counts, parameter_area))
}
likelihood$gradlogdensity <- function(x){
  return( data_counts - parameter_area * exp(x) )
}

# posterior distribution
logtarget <- function(x) prior$logdensity(x) + likelihood$logdensity(x)
gradlogtarget <- function(x) prior$gradlogdensity(x) + likelihood$gradlogdensity(x)

# initial distribution
rinit <- function() fast_rmvnorm(1, prior_mean, prior_cov)

# Hamiltonian function
hamiltonian <- function(x, v) -logtarget(x) + sum(v^2) / 2

# stepsize selection
nsteps <- 100
# grid_stepsize <- seq(1e-2, 1e-1, length.out = 10)
grid_stepsize <- c(1e-2, 1e-1)
nstepsizes <- length(grid_stepsize)
# nreps <- 100
nreps <- 2 # test
hamiltonian_error <- matrix(nrow = nreps, ncol = nstepsizes)
tic("Runtime:")
for (igrid in 1:nstepsizes){
  igrid
  cat("Grid:", igrid, "/", nstepsizes, "\n")
  stepsize <- grid_stepsize[igrid]

  for (irep in 1:nreps){
    cat("Repetition:", irep, "/", nreps, "\n")
    # store hamiltonian error along trajectory
    h_error <- rep(0, nsteps)

    # initialize
    current_x <- rinit()
    current_v <- rnorm(dimension)
    initial_hamiltonian <- hamiltonian(current_x, current_v)

    # perform leapfrog integration
    for (istep in 1:nsteps){
      current_v <- current_v + stepsize * gradlogtarget(current_x) / 2
      current_x <- current_x + stepsize * current_v
      current_v <- current_v + stepsize * gradlogtarget(current_x) / 2
      h_error[istep] <- abs( hamiltonian(current_x, current_v) - initial_hamiltonian )
    }

    # store maximum hamiltonian error along trajectory
    hamiltonian_error[irep, igrid] <- max(h_error)
  }

}
toc()

# save(nstepsizes, grid_stepsize, nreps, hamiltonian_error, file = "inst/coxprocess/hmc.stepsize.selection.RData")





