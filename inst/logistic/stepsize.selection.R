rm(list = ls())
library(tictoc)
library(debiasedhmc)

# load("inst/logistic/germancredit.RData")
# load("inst/logistic/simulated_nsamples1000_dimension66.RData")
# load("inst/logistic/simulated_nsamples1000_dimension130.RData")
# load("inst/logistic/simulated_nsamples1000_dimension258.RData")
load("inst/logistic/simulated_nsamples1000_dimension514.RData")
# load("inst/logistic/simulated_nsamples1000_dimension1026.RData")

nmcmc <- 50
stepsize <- 0.5
nsteps <- 30

hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

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
  # cat("No. of steps", nsteps, "Iteration:", imcmc, "/", nmcmc, "\n")
}
timing <- toc()
runtime <- timing$toc - timing$tic
cat("Accept prob:", accept / nmcmc , "\n")
cat("Time per iteration:", runtime / nmcmc , "\n")
# acf(chain[, 1])
