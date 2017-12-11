# process repeat results for various k,m
rm(list=ls())

parameter_k <- 1 # median = 208, 90% quantile = 330
parameter_m <- 1 * parameter_k # factor = 1, 5, 10

ntestfunct <- 604
njobs <- 100
nreps <- 10
nrepeats <- njobs * nreps
estimates <- matrix(nrow = nrepeats, ncol = ntestfunct)
cost <- rep(0, nrepeats)

for (ijob in 1:njobs){
  filename <- paste("output/output.hmc.germancredit.repeat", ijob, ".RData", sep = "")
  load(file = filename)
  for (irep in 1:nreps){
    cat((ijob-1)*nreps+irep, "\n")
    if (parameter_k == parameter_m){
      estimates[(ijob-1)*nreps+irep, ] <- unbiased_estimates[[irep]][, parameter_k]
    } else {
      estimates[(ijob-1)*nreps+irep, ] <- rowMeans(unbiased_estimates[[irep]][, parameter_k:parameter_m])
    }
     # time averaged estimator
    cost[(ijob-1)*nreps+irep] <- 2 * meetingtime[irep] + max(1, parameter_m + 1 - meetingtime[irep])
  }
}

load("output/optimal.hmc.germancredit.RData") # load HMC inefficiencies
mean_cost <- mean(cost)
variance <- sum(apply(estimates, 2, var))
inefficiency <- mean_cost * variance
cat("cost =", mean_cost, "\n")
cat("variance =", variance, "\n")
cat("inefficiency =", inefficiency, "\n")
cat("relative inefficiency =", inefficiency / suboptimal_hmc$asymptotic_variance, "\n")
cat("relative HMC inefficiency =", suboptimal_hmc$asymptotic_variance / optimal_hmc$asymptotic_variance, "\n")

