# get parameters with desired contraction
load(file = "inst/coxprocess/output.hmc.contraction.RData")
index_contraction <- (results.df$distance < -10 )
contraction.df <- unique(results.df[index_contraction, c("stepsize", "nsteps")])

# get corresponding ESS
load(file = "inst/coxprocess/output.hmc.tuning.RData")
for (i in 1:nrow(contraction.df)){
  stepsize <- contraction.df[i, "stepsize"]
  nsteps <- contraction.df[i, "nsteps"]
  index <- (results.df$stepsize == stepsize) & (results.df$nsteps == nsteps)
  contraction.df$ess[i] <- results.df$ess[index]
}

# deterministic policy
index_max <- (contraction.df$ess == max(contraction.df$ess))
contraction.df[index_max, c("stepsize", "nsteps")]

# random policy by assigning weights proportional to ESS
contraction.df$weights = contraction.df$ess / sum(contraction.df$ess)
contraction.df
