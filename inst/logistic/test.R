rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)


# load model
load("inst/logistic/germancredit.RData")

# specify stepsize and no. of steps
stepsize <- 0.0125
nsteps <- 10

# compute distance
compute_distance <- function(cchain){
  niterations <- dim(cchain$samples2)[1]
  return(sapply(1:niterations, function(index) sqrt(sum((cchain$samples1[1+index, ] - cchain$samples2[index, ])^2))))
}

# no. of repetitions and mcmc iterations
m <- 500
max_iterations <- 500

# define hmc kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# define mixture kernels
omega <- 1 / 20
Sigma_std <- 1e-3
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension) # try 1e-3, 1e-5, 1e-7
mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension)

# mixture kernels
mixture_kernel <- function(chain_state, current_pdf, iteration){
  if (runif(1) < omega){
    output <- mh$kernel(chain_state, current_pdf, iteration)
    output$rwmh <- TRUE
  } else {
    output <- hmc$kernel(chain_state, current_pdf, iteration)
    output$rwmh <- FALSE
  }
  return(output)
}

mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration){
  if (runif(1) < omega){
    output <- mh$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration)
    output$rwmh <- TRUE
  } else {
    output <- hmc$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration)
    output$rwmh <- FALSE
    output$overlap <- FALSE
  }
  return(output)
}

tic()
cchains <- coupled_chains(logtarget, mixture_kernel, mixture_coupled_kernel, rinit, m = m, max_iterations = max_iterations)
toc()
distance <- log10(compute_distance(cchains))
meetingtime <- cchains$meetingtime

time_iteration <- 1:max_iterations
distance.df <- data.frame(iteration = time_iteration, distance = distance)
rwmh.df <- data.frame(iteration = time_iteration[cchains$logical_rwmh], distance = distance[cchains$logical_rwmh])
overlap.df <- data.frame(iteration = time_iteration[cchains$logical_overlap], distance = distance[cchains$logical_overlap])
logical_reject <- !(cchains$logical_accept1 & cchains$logical_accept2) # either chains rejecting
reject.df <- data.frame(iteration = time_iteration[logical_reject], distance = distance[logical_reject])

g <- ggplot() +
  geom_point(data = distance.df, aes(x = iteration, y = distance)) +
  geom_point(data = rwmh.df, aes(x = iteration, y = distance), color = "red") +
  geom_point(data = overlap.df, aes(x = iteration, y = distance), color = "green") +
  geom_point(data = reject.df, aes(x = iteration, y = distance), color = "blue") +
  xlab('iteration') + ylab('distance')  #+ xlim(0, max_iterations) + scale_y_log10()
g

# times where rwmh moves + rejected
time_iteration[cchains$logical_rwmh & logical_reject]

# times where rwmh moves + overlap
time_iteration[cchains$logical_overlap]
