# process output files
rm(list=ls())
grid_stepsize <- seq(0.01, 0.04, by = 0.0025)
ngrid_stepsize <- length(grid_stepsize)
grid_nsteps <- c(10, 20, 30)
ngrid_nsteps <- length(grid_nsteps)
results.df <- data.frame()
for (igrid in 1:ngrid_stepsize){
  filename <- paste("inst/logistic/output/output.hmc.germancredit.tuning", igrid, ".RData", sep = "")
  load(file = filename)
  results.df <- rbind(results.df, data.frame(stepsize = rep(grid_stepsize[igrid], ngrid_nsteps),
                                             nsteps = grid_nsteps,
                                             asymptotic_variance = variance,
                                             runtime = runtimes,
                                             acceptprob = acceptprob))
  # file.remove(file = filename)
}
# load("inst/logistic/output/output.mala.germancredit.tuning.RData")
# results.df <- rbind(results.df, data.frame(stepsize = grid_stepsize,
                                           # nsteps = rep(1, ngrid_stepsize),
                                           # asymptotic_variance = variance,
                                           # runtime = runtimes,
                                           # acceptprob = acceptprob))

save(results.df, file = "inst/logistic/output/output.hmc.germancredit.tuning.RData")

rm(list=ls())
library(ggplot2)
setmytheme()
# order by asymptotic efficiency
load("inst/logistic/output/output.hmc.germancredit.tuning.RData")
asymptotic_efficiency <- results.df$nsteps * results.df$asymptotic_variance # leapfrog steps as cost
# asymptotic_efficiency <- results.df$runtime * results.df$asymptotic_variance # runtime as cost
sorted_efficiency <- sort(asymptotic_efficiency, index.return = TRUE)
results.df[sorted_efficiency$ix[1:5], ]

# Optimal marginal HMC kernel
cat("Optimal marginal HMC kernel")
index_min <- sorted_efficiency$ix[1]
optimal_hmc <- results.df[index_min, ] # optimal hmc parameters
suboptimal_hmc <- results.df[4, ] # parameters used for contraction
save(suboptimal_hmc, optimal_hmc, file = "inst/logistic/output/optimal.hmc.germancredit.RData")

# variance plot (no colour)
g <- ggplot(results.df, aes(x = stepsize, y = asymptotic_variance)) +
  geom_point(aes(shape = factor(nsteps))) +
  scale_shape_discrete(name = "steps") +
  geom_line(aes(linetype = factor(nsteps))) +
  scale_linetype_discrete(name = "steps") +
  labs(x = "stepsize", y = "variance")
g
ggsave(filename = "inst/logistic/plots/logistic.hmc.germancredit.tuning.pdf", plot = g, width = 9, height = 6)


# variance plot (in colour)
g <- ggplot(results.df, aes(x = stepsize, y = asymptotic_variance)) +
  geom_line(aes(colour = factor(nsteps))) +
  geom_point(aes(colour = factor(nsteps))) +
  scale_colour_discrete(name = "steps") +
  labs(x = "stepsize", y = "variance")
g
ggsave(filename = "inst/logistic/plots/logistic.hmc.germancredit.tuning.pdf", plot = g, width = 9, height = 6)

# acceptance probability (no colour)
ggplot(results.df, aes(x = stepsize, y = acceptprob)) +
  geom_point(aes(shape = factor(nsteps))) +
  scale_shape_discrete(name = "steps") +
  geom_line(aes(linetype = factor(nsteps))) +
  scale_linetype_discrete(name = "steps") +
  labs(x = "stepsize", y = "acceptance probability")


# acceptance probability (in colour)
ggplot(results.df, aes(x = stepsize, y = acceptprob)) +
  geom_line(aes(colour = factor(nsteps))) +
  geom_point(aes(colour = factor(nsteps))) +
  scale_colour_discrete(name = "steps") +
  labs(x = "stepsize", y = "acceptance probability")

