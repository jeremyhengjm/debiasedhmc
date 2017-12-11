# estimate variance reduction weights
rm(list = ls())
library(tictoc)

parameter_k <- 330 # median = 208, 90% quantile = 330
parameter_m <- 10 * parameter_k # factor = 1, 5, 10

ntestfunct <- 604
njobs <- 100
nreps <- 10
nrepeats <- njobs * nreps
weight <- matrix(nrow = ntestfunct, ncol = parameter_m-parameter_k+1)

for (itestfunc in 1:ntestfunct){
  estimates <- matrix(nrow = nrepeats, ncol = parameter_m-parameter_k+1)
  for (ijob in 1:njobs){
    filename <- paste("output/output.hmc.germancredit.repeat", ijob, ".RData", sep = "")
    load(file = filename)
    for (irep in 1:nreps){
      # cat("Test function:", itestfunc, "/", ntestfunct,
          # "Repetition:", (ijob-1)*nreps+irep, "/", nrepeats, "\n")
      estimates[(ijob-1)*nreps+irep, ] <- unbiased_estimates[[irep]][itestfunc, parameter_k:parameter_m]
    }
  }
  tic()
  cat("Test function:", itestfunc, "/", ntestfunct, "\n")
  cat("Condition number:", rcond(cbind(rbind(cov(estimates), rep(1, parameter_m-parameter_k+1)),
              c(rep(1, parameter_m-parameter_k+1), 0))), "\n")
  # weight[itestfunc, ] <- solve(cbind(rbind(cov(estimates), rep(1, parameter_m-parameter_k+1)),
                                     # c(rep(1, parameter_m-parameter_k+1), 0)),
                               # c(rep(0, parameter_m-parameter_k+1), 1))[1:(parameter_m-parameter_k+1)]
  toc()
  # save(weight, file = "variance.reduction.germancredit.RData")

}








