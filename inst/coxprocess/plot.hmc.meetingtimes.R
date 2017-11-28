# process output files
rm(list=ls())
library(ggplot2)
njobs <- 50
nreps <- 2
nrepeats <- njobs * nreps

# load all parameter configurations
load("inst/coxprocess/hmc.contraction.parameters.RData")

# process output files
for (iparameter in 2:nrow(contraction.df)){
  # current configuration
  stepsize <- contraction.df$stepsize[iparameter]
  nsteps <- contraction.df$nsteps[iparameter]

  # pre-allocate
  meetingtimes <- rep(0, nrepeats)
  distance.df <- data.frame()
  for (ijob in 1:njobs){
    filename <- paste("inst/coxprocess/output.hmc.parameter", iparameter, ".rep", ijob, ".RData", sep = "")
    load(file = filename)

    for (irep in 1:nreps){
    meetingtimes[(ijob-1)*nreps+irep] <- meetingtime[irep]
    distance.df <- rbind(distance.df, data.frame(iteration = 1:max_iterations,
                                                 repetition = rep((ijob-1)*nreps+irep, max_iterations),
                                                 distance = distance[[irep]]))
    }
    file.remove(file = filename)
  }
  combined_filename <- paste("inst/coxprocess/output.hmc.parameter", iparameter, ".RData", sep = "")
  save(nrepeats, stepsize, nsteps, meetingtimes, distance.df, file = combined_filename)
}

# plotting
load("inst/coxprocess/hmc.contraction.parameters.RData")
nrepeats <- 100
ghist <- list()
gdistance <- list()
mean_cost <- rep(0, nrow(contraction.df))
store_meetingtimes <- matrix(nrow = nrow(contraction.df), ncol = nrepeats)
for (iparameter in 1:nrow(contraction.df)){
  # current configuration
  combined_filename <- paste("inst/coxprocess/output.hmc.parameter", iparameter, ".RData", sep = "")
  load(file = combined_filename)
  # cat("stepsize =", contraction.df$stepsize[iparameter], ",", "steps =", contraction.df$nsteps[iparameter], "\n")
  glabel <- paste("stepsize = ", contraction.df$stepsize[iparameter], ", ",
                  "steps = ", contraction.df$nsteps[iparameter], sep = "")

  # meeting times
  store_meetingtimes[iparameter, ] <- meetingtimes
  meetingtimes.df <- data.frame(meetingtime = meetingtimes)
  mean_cost[iparameter] <- mean(meetingtimes * contraction.df$nsteps[iparameter])

  g <- ggplot(data = meetingtimes.df, aes(x = meetingtime)) +
        geom_histogram(aes(y = ..density..), binwidth = 50) +
        xlab("meeting times") + xlim(0, 1000) +
        ggtitle(label = glabel)
  ghist[[iparameter]] <- g

  # distance over iterations
  g <- ggplot(distance.df, aes(x = iteration, y = distance, group = repetition)) + geom_line(alpha = 0.25) +
        xlab('iteration') + ylab('distance') + xlim(0, 1000) +
        scale_y_log10() + #breaks = c(1e-5, 1e-3, 1e-1, 10), limits = c(1e-6, 3))
        ggtitle(label = glabel)
  gdistance[[iparameter]] <- g

  # choose k and m parameters
  contraction.df$k[iparameter] <- floor(quantile(meetingtimes, 0.9))
  contraction.df$m[iparameter] <- 10 * contraction.df$k[iparameter]
}
contraction.df
rowMeans(store_meetingtimes)
mean_cost
index_min <- (mean_cost == min(mean_cost))
contraction.df[index_min, ]
summary(store_meetingtimes[index_min, ])

# plot meeting times
grid.arrange(arrangeGrob(ghist[[1]], ghist[[2]], ghist[[3]],
                         ghist[[4]], ghist[[5]], ghist[[6]],
                         nrow = 3, ncol = 2))
# g_filename <- paste("inst/coxprocess/hmc.meetingtimes.pdf")
# ggsave(filename = g_filename, plot = g, width = 9, height = 6)

# plot distance over iterations
grid.arrange(arrangeGrob(gdistance[[1]], gdistance[[2]], gdistance[[3]],
                              gdistance[[4]], gdistance[[5]], gdistance[[6]],
                              nrow = 3, ncol = 2))
# g_filename <- paste("inst/coxprocess/hmc.distance.pdf")
# ggsave(filename = g_filename, plot = g, width = 9, height = 6)


