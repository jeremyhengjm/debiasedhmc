# process output files
rm(list=ls())
library(ggplot2)
nreps <- 5
grid_stepsize <- seq(0.05, 0.45, by = 0.02)
ngrid_stepsize <- length(grid_stepsize)
grid_nsteps <- c(1, 10, 20, 30) # c(0.4, 1.2, 2.4, 3.6) seconds per mcmc iteration
ngrid_nsteps <- length(grid_nsteps)
results.df <- data.frame()
for (igrid in 1:ngrid_stepsize){
  filename <- paste("inst/coxprocess/output.hmc.contraction", igrid, ".RData", sep = "")
  load(file = filename)
  # ngrid_nsteps <- 1
  for (istep in 1:ngrid_nsteps){
  results.df <- rbind(results.df, data.frame(stepsize = rep(grid_stepsize[igrid], nreps),
                                             nsteps = rep(grid_nsteps[istep], nreps),
                                             time = rep(grid_stepsize[igrid] * grid_nsteps[istep], nreps),
                                             distance = log(distance[, istep])))
  }
  # file.remove(file = filename)
}
save(nreps, grid_stepsize, ngrid_stepsize, grid_nsteps, ngrid_nsteps,
     results.df, file = "inst/coxprocess/output.hmc.contraction.RData")

# distance after 1000 iterations (no colour)
load(file = "inst/coxprocess/output.hmc.contraction.RData")
g <- ggplot(results.df, aes(x = stepsize, y = distance)) +
  geom_point(aes(shape = factor(nsteps))) +
  scale_shape_discrete(name = "steps") +
  labs(x = "stepsize", y = "log distance after 1000 iterations")
g
ggsave(filename = "inst/coxprocess/coxprocess.hmc.contraction.pdf", plot = g, width = 9, height = 6)

# distance after 1000 iterations (in colour)
g <- ggplot(results.df, aes(x = stepsize, y = distance)) +
  geom_point(aes(colour = factor(nsteps))) +
  scale_colour_discrete(name = "steps") +
  labs(x = "stepsize", y = "log distance after 1000 iterations")
g
ggsave(filename = "inst/coxprocess/coxprocess.hmc.contraction.pdf", plot = g, width = 9, height = 6)

# distance against integration time
ggplot(results.df, aes(x = time, y = distance)) +
  geom_point(aes(colour = factor(nsteps))) +
  scale_colour_discrete(name = "steps") +
  labs(x = "integration time", y = "log distance after 1000 iterations")
