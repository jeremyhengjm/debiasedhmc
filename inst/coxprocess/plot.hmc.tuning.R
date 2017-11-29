# process output files
rm(list=ls())
library(ggplot2)
dimension <- 64^2
grid_stepsize <- seq(0.05, 0.45, by = 0.02)
ngrid_stepsize <- length(grid_stepsize)
grid_nsteps <- c(1, 10, 20, 30)
ngrid_nsteps <- length(grid_nsteps)
results.df <- data.frame()
for (igrid in 1:ngrid_stepsize){
  filename <- paste("inst/coxprocess/output/output.hmc.tuning", igrid, ".RData", sep = "")
  load(file = filename)
  results.df <- rbind(results.df, data.frame(stepsize = rep(grid_stepsize[igrid], ngrid_nsteps),
                                     nsteps = grid_nsteps,
                                     var = variance * grid_nsteps / dimension,
                                     acceptprob = acceptprob))
  # file.remove(file = filename)
}
save(grid_stepsize, ngrid_stepsize, grid_nsteps, ngrid_nsteps,
     results.df, file = "inst/coxprocess/output/output.hmc.tuning.RData")

# Optimal marginal HMC kernel
load("inst/coxprocess/output/output.hmc.tuning.RData")
index_min <- (results.df$var == min(results.df$var))
cat("Optimal marginal HMC kernel")
results.df[index_min, ]

# variance plot (no colour)
ggplot(results.df, aes(x = stepsize, y = var)) +
  geom_point(aes(shape = factor(nsteps))) +
  scale_shape_discrete(name = "steps") +
  geom_line(aes(linetype = factor(nsteps))) +
  scale_linetype_discrete(name = "steps") +
  labs(x = "stepsize", y = "variance")
g
ggsave(filename = "inst/coxprocess/plots/coxprocess.hmc.tuning.pdf", plot = g, width = 9, height = 6)


# variance plot (in colour)
g <- ggplot(results.df, aes(x = stepsize, y = var)) +
  geom_line(aes(colour = factor(nsteps))) +
  geom_point(aes(colour = factor(nsteps))) +
  scale_colour_discrete(name = "steps") +
  labs(x = "stepsize", y = "variance")
g
ggsave(filename = "inst/coxprocess/plots/coxprocess.hmc.tuning.pdf", plot = g, width = 9, height = 6)

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

