rm(list = ls())
library(ggplot2)
load("inst/coxprocess/hmc.stepsize.selection.RData")
hamiltonian_error.df <- data.frame()
for (igrid in 1:nstepsizes){
  hamiltonian_error.df <- rbind(hamiltonian_error.df, data.frame(stepsize = rep(grid_stepsize[igrid], nreps),
                                                                 error = hamiltonian_error[, igrid] ) )
}
ggplot() + geom_point(data = hamiltonian_error.df, aes(x = stepsize, y = error))
