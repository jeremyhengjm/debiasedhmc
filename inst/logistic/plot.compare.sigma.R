# process output files
rm(list=ls())
njobs <- 10
nreps <- 10
nrepeats <- njobs * nreps

# process output files
# pre-allocate
meetingtimes <- rep(0, nrepeats)
distance.df <- data.frame()
store_distance <- matrix(nrow = nrepeats, ncol = 1000)
for (ijob in 1:njobs){
  filename <- paste("inst/logistic/output/output.hmc.germancredit.betterinit", ijob, ".RData", sep = "")
  load(file = filename)

  for (irep in 1:nreps){
    meetingtimes[(ijob-1)*nreps+irep] <- meetingtime[irep]
    distance.df <- rbind(distance.df, data.frame(iteration = 1:max_iterations,
                                                 repetition = rep((ijob-1)*nreps+irep, max_iterations),
                                                 distance = distance[[irep]]))
    store_distance[(ijob-1)*nreps+irep, ] <- distance[[irep]]
  }
  file.remove(file = filename)
}
combined_filename <- paste("inst/logistic/output/output.hmc.germancredit.betterinit.RData")
save(nrepeats, max_iterations, stepsize, nsteps, Sigma_std,
     meetingtimes, distance.df, file = combined_filename)

# plot result with better initialization
rm(list=ls())
library(ggplot2)
library(gridExtra)
load("inst/logistic/output/output.hmc.germancredit.betterinit.RData")
summary(meetingtimes)
floor(quantile(meetingtimes, 0.9))
meetingtimes.df <- data.frame(meetingtime = meetingtimes)
ggplot(data = meetingtimes.df, aes(x = meetingtime)) +
  geom_histogram(aes(y = ..density..), binwidth = 50) +
  xlab("meeting times") + xlim(0, 1000)
sum(meetingtimes > 500)
ggplot(distance.df, aes(x = iteration, y = distance, group = repetition)) +
  geom_line(alpha = 0.25) +
  xlab('iteration') + ylab('distance') + xlim(0, 1000) +
  scale_y_log10(breaks = 10^seq(2, -10, by = -2))

# plotting results for different sigma
rm(list=ls())
library(ggplot2)
library(gridExtra)
load("inst/logistic/output/output.hmc.germancredit.sigma1e-1.RData")
summary(meetingtimes)
meetingtimes.df <- data.frame(meetingtime = meetingtimes)
ghist1 <- ggplot(data = meetingtimes.df, aes(x = meetingtime)) +
  geom_histogram(aes(y = ..density..), binwidth = 50) +
  xlab("meeting times") + xlim(0, 1000)
sum(meetingtimes > 500)
gdist1 <- ggplot(distance.df, aes(x = iteration, y = distance, group = repetition)) +
  geom_line(alpha = 0.25) +
  xlab('iteration') + ylab('distance') + xlim(0, 1000) +
  scale_y_log10(breaks = 10^seq(2, -10, by = -2))

load("inst/logistic/output/output.hmc.germancredit.sigma1e-3.RData")
summary(meetingtimes)
meetingtimes.df <- data.frame(meetingtime = meetingtimes)
ghist3 <- ggplot(data = meetingtimes.df, aes(x = meetingtime)) +
  geom_histogram(aes(y = ..density..), binwidth = 50) +
  xlab("meeting times") + xlim(0, 1000)
sum(meetingtimes > 500)
gdist3 <- ggplot(distance.df, aes(x = iteration, y = distance, group = repetition)) +
  geom_line(alpha = 0.25) +
  xlab('iteration') + ylab('distance') + xlim(0, 1000) +
  scale_y_log10(breaks = 10^seq(2, -10, by = -2))

load("inst/logistic/output/output.hmc.germancredit.sigma1e-5.RData")
summary(meetingtimes)
meetingtimes.df <- data.frame(meetingtime = meetingtimes)
ghist5 <- ggplot(data = meetingtimes.df, aes(x = meetingtime)) +
  geom_histogram(aes(y = ..density..), binwidth = 50) +
  xlab("meeting times") + xlim(0, 1000)
sum(meetingtimes > 500)
gdist5 <- ggplot(distance.df, aes(x = iteration, y = distance, group = repetition)) +
  geom_line(alpha = 0.25) +
  xlab('iteration') + ylab('distance') + xlim(0, 1000) +
  scale_y_log10(breaks = 10^seq(2, -10, by = -2))

load("inst/logistic/output/output.hmc.germancredit.sigma1e-7.RData")
summary(meetingtimes)
meetingtimes.df <- data.frame(meetingtime = meetingtimes)
ghist7 <- ggplot(data = meetingtimes.df, aes(x = meetingtime)) +
  geom_histogram(aes(y = ..density..), binwidth = 50) +
  xlab("meeting times") + xlim(0, 1000)
sum(meetingtimes > 500)
gdist7 <- ggplot(distance.df, aes(x = iteration, y = distance, group = repetition)) +
  geom_line(alpha = 0.25) +
  xlab('iteration') + ylab('distance') + xlim(0, 1000) +
  scale_y_log10(breaks = 10^seq(2, -10, by = -2))

grid.arrange(arrangeGrob(ghist3, ghist5, ghist7, nrow = 1, ncol = 3))
grid.arrange(arrangeGrob(gdist3, gdist5, gdist7, nrow = 1, ncol = 3))

