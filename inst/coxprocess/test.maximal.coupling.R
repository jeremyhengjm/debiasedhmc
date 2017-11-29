rm(list = ls())
# load cox process model
load("inst/coxprocess/coxprocess.RData")

# specify stepsize and no. of steps
load("inst/coxprocess/hmc.contraction.parameters.RData")
iparameter <- 1
stepsize <- contraction.df$stepsize[iparameter]
nsteps <- contraction.df$nsteps[iparameter]

# define hmc kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# define mixture kernels
omega <- 1 / 20

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

# check maximal coupling
Sigma_std <- 1e-3
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension)
mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension)
threshold <- 1e-10
nreps <- 10
meet <- 0
tic()
for (irep in 1:nreps){
  current_x1 <- rinit()
  current_x2 <- current_x1 + rep(threshold / sqrt(dimension) , dimension) # chains are distance threshold apart

  current_pdf1 <- logtarget(current_x1)
  current_pdf2 <- logtarget(current_x2)
  output_mh <- mh$coupled_kernel(current_x1, current_x2, current_pdf1, current_pdf2, 0)
  current_x1 <- output_mh$chain_state1
  current_x2 <- output_mh$chain_state2
  meet <- meet + all(current_x1 == current_x2)

}
toc()
cat("Meeting probability:", meet / nreps)
meet / nreps
