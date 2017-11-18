# process output files
nstepsizes <- 10
mean_error.df <- data.frame()
hamiltonian_error.df <- data.frame()
for (igrid in 1:nstepsizes){
  filename <- paste("inst/coxprocess/output.hmc.stepsize.selection.", igrid, ".RData", sep = "")
  load(file = filename)
  mean_error.df <- rbind(mean_error.df, data.frame(stepsize = grid_stepsize[igrid], error = mean(hamiltonian_error)))
  hamiltonian_error.df <- rbind(hamiltonian_error.df, data.frame(stepsize = rep(grid_stepsize[igrid], nreps),
                                                                 error = hamiltonian_error))
  file.remove(file = filename)
}
save(nstepsizes, grid_stepsize, nreps, mean_error.df, hamiltonian_error.df, file = "inst/coxprocess/output.hmc.stepsize.selection.RData", safe = F)

# plot results
rm(list = ls())
library(ggplot2)
load(file = "inst/coxprocess/output.hmc.stepsize.selection.RData")
threshold.df <- data.frame(stepsize = c(grid_stepsize[1]-0.005, grid_stepsize[nstepsizes]+0.005), error = rep(abs(log(0.8)),2))
ggplot() + geom_point(data = hamiltonian_error.df, aes(x = stepsize, y = error)) +
  geom_line(data = threshold.df, aes(x = stepsize, y = error), colour = "red") +
  labs(x = "stepsize", y = "Hamiltonian error")

ggplot() + geom_point(data = mean_error.df, aes(x = stepsize, y = error)) +
  geom_line(data = threshold.df, aes(x = stepsize, y = error), colour = "red") +
  labs(x = "stepsize", y = "Average Hamiltonian error")
