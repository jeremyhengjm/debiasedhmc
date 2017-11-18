# process output files
rm(list=ls())
library(ggplot2)
grid_stepsize <- seq(0.05, 0.45, by = 0.02)
# grid_stepsize <- seq(0.05, 0.25, length.out = 11)
ngrid_stepsize <- length(grid_stepsize)
grid_nsteps <- c(1, 10, 20, 30) # c(0.4, 1.2, 2.4, 3.6) seconds per mcmc iteration
ngrid_nsteps <- length(grid_nsteps)
results.df <- data.frame()
for (igrid in 1:ngrid_stepsize){
  filename <- paste("inst/coxprocess/output.hmc.tuning", igrid, ".RData", sep = "")
  load(file = filename)
  results.df <- rbind(results.df, data.frame(stepsize = rep(grid_stepsize[igrid], ngrid_nsteps),
                                     nsteps = grid_nsteps,
                                     ess = ess,
                                     acceptprob = acceptprob))
  # file.remove(file = filename)
}
results.df$weights <- results.df$ess / sum(results.df$ess) # random policy by assigning weights proportional to ESS
save(grid_stepsize, ngrid_stepsize, grid_nsteps, ngrid_nsteps,
     results.df, file = "inst/coxprocess/output.hmc.tuning.RData")

# ess per second plot (no colour)
load("inst/coxprocess/output.hmc.tuning.RData")
ggplot(results.df, aes(x = stepsize, y = ess)) +
  geom_point(aes(shape = factor(nsteps))) +
  scale_shape_discrete(name = "steps") +
  geom_line(aes(linetype = factor(nsteps))) +
  scale_linetype_discrete(name = "steps") +
  labs(x = "stepsize", y = "effective sample size / sec")
g
ggsave(filename = "coxprocess.hmc.tuning.pdf", plot = g, width = 9, height = 6)


# ess per second plot (in colour)
g <- ggplot(results.df, aes(x = stepsize, y = ess)) +
  geom_line(aes(colour = factor(nsteps))) +
  geom_point(aes(colour = factor(nsteps))) +
  scale_colour_discrete(name = "steps") +
  labs(x = "stepsize", y = "effective sample size / sec")
g
ggsave(filename = "inst/coxprocess/coxprocess.hmc.tuning.pdf", plot = g, width = 9, height = 6)

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

