# process output files
rm(list=ls())
library(ggplot2)
library(gridExtra)
nrepeats <- 100

# load all parameter configurations
load("inst/logistic/hmc.germancredit.contraction.parameters.RData")

# plotting
nrepeats <- 100
ghist <- list()
gdistance <- list()
mean_cost <- rep(0, nrow(contraction.df))
store_meetingtimes <- matrix(nrow = nrow(contraction.df), ncol = nrepeats)
for (iparameter in 1:nrow(contraction.df)){
  # current configuration
  combined_filename <- paste("inst/logistic/output/output.hmc.germancredit.parameter", iparameter, ".RData", sep = "")
  load(file = combined_filename)
  glabel <- paste("stepsize = ", contraction.df$stepsize[iparameter], ", ",
                  "steps = ", contraction.df$nsteps[iparameter], sep = "")

  # meeting times
  meetingtimes <- meetingtime
  store_meetingtimes[iparameter, ] <- meetingtimes
  meetingtimes.df <- data.frame(meetingtime = meetingtimes)
  mean_cost[iparameter] <- mean(meetingtimes * contraction.df$nsteps[iparameter])

  g <- ggplot(data = meetingtimes.df, aes(x = meetingtime)) +
    geom_histogram(aes(y = ..density..), binwidth = 50) +
    xlab("meeting times") + xlim(0, 1000) # + ggtitle(label = glabel)
  ghist[[iparameter]] <- g

  distance.df <- data.frame()
  # distance over iterations
  for (irep in 1:nrepeats){
    distance.df <- rbind(distance.df, data.frame(iteration = 1:max_iterations,
                                                 repetition = rep(irep, max_iterations),
                                                 distance = distance[[irep]]))
  }
  g <- ggplot(distance.df, aes(x = iteration, y = distance, group = repetition)) +
    geom_line(alpha = 0.25) +
    xlab('iteration') + ylab('distance') + xlim(0, 1000) +
    scale_y_log10() # + ggtitle(label = glabel)
  gdistance[[iparameter]] <- g

  # choose k and m parameters
  contraction.df$k[iparameter] <- floor(quantile(meetingtimes, 0.9))
  contraction.df$m[iparameter] <- 10 * contraction.df$k[iparameter]
}
contraction.df
rowMeans(store_meetingtimes)
mean_cost
index_min <- 1:nrow(contraction.df)
logical_min <- (mean_cost == min(mean_cost))
index_min <- index_min[logical_min]
contraction.df[index_min, ]
summary(store_meetingtimes[index_min, ])

# plot meeting times and distance for best configuration
ghist[[index_min]]
gdistance[[index_min]]

# plot meeting times
grid.arrange(arrangeGrob(ghist[[1]], ghist[[2]], ghist[[3]],
                         ghist[[4]], ghist[[5]], ghist[[6]],
                         ghist[[7]], ghist[[8]], ghist[[9]], ghist[[10]],
                         nrow = 5, ncol = 2))
# g_filename <- paste("inst/logistic/plots/hmc.meetingtimes.pdf")
# ggsave(filename = g_filename, plot = g, width = 9, height = 6)

# plot distance over iterations
grid.arrange(arrangeGrob(gdistance[[1]], gdistance[[2]], gdistance[[3]],
                         gdistance[[4]], gdistance[[5]], gdistance[[6]],
                         gdistance[[7]], gdistance[[8]], gdistance[[9]], gdistance[[10]],
                         nrow = 5, ncol = 2))
# g_filename <- paste("inst/logistic/plots/hmc.distance.pdf")
# ggsave(filename = g_filename, plot = g, width = 9, height = 6)


