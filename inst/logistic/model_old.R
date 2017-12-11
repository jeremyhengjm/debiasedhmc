rm(list=ls())
library(debiasedhmc)
library(tictoc)
library(expm)

# simulate dataset
# dimension <- 300
# nsamples <- 1000
# true_parameter <- rep(1, dimension)
# logistic_function <- function(x) 1 / (1 + exp(-x))
# design_matrix <- matrix(nrow = nsamples, ncol = dimension)
# response <- rep(0, nsamples)
# for (isample in 1:nsamples){
#   # sample from Rademacher distribution
#   covariate <- sample(c(-1,1), size = dimension, replace = TRUE, prob = rep(0.5, 2))
#   # normalize
#   covariate <- covariate / sqrt(sum(covariate^2))
#   design_matrix[isample, ] <- covariate
#   response[isample] <- ( runif(1) < logistic_function(sum(true_parameter * covariate)) )
# }

rm(list=ls())
# german credit dataset
germancredit <- read.table(system.file("", "germancredit.txt", package = "debiasedhmc"))
design_matrix <- scale(germancredit[, 1:24])
response <- germancredit[, 25] - 1
dimension <- ncol(design_matrix)
nsamples <- nrow(design_matrix)
interaction_terms <- matrix(nrow = nsamples, ncol = dimension*(dimension-1) / 2)
index <- 1
for (j in 1:(dimension-1)){
  for (jprime in (j+1):dimension){
    interaction_terms[, index] <- design_matrix[, j] * design_matrix[, jprime]
    index <- index + 1
  }
}
design_matrix <- cbind(rep(1, nsamples), design_matrix, scale(interaction_terms))
colnames(design_matrix) <- NULL
dimension <- ncol(design_matrix)

# prior distribution
prior_mean <- rep(0, dimension)
prior_sigmasq <- 1
# prior_precision <- (3 * dimension / pi^2) * (t(design_matrix) %*% design_matrix) / nsamples
prior_precision <- diag(prior_sigmasq, dimension, dimension)
prior_cov <- solve(prior_precision)
prior_precision_chol <- t(chol(prior_precision))
prior <- list()
prior$logdensity <- function(x){
  return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), prior_mean, prior_precision_chol))
}
prior$gradlogdensity <- function(x){
  return(gradlognormal(x, prior_mean, prior_precision))
}

# pre-conditioning
metric <- 0.25 * (t(design_matrix) %*% design_matrix) + prior_precision
(eigen(gram_matrix)$values > 0)
preconditioner <- solve(sqrtm(gram_matrix))


# likelihood function
response_design <- response %*% design_matrix # of length dimension
likelihood <- list()
likelihood$logdensity <- function(x){
  return(logistic_likelihood(x, design_matrix, response))
}
likelihood$gradlogdensity <- function(x){
  return(logistic_gradloglikelihood(x, design_matrix, response_design))
}

# posterior distribution
logtarget <- function(x) prior$logdensity(x) + likelihood$logdensity(x)
gradlogtarget <- function(x) prior$gradlogdensity(x) + likelihood$gradlogdensity(x)

# pre-conditioned posterior distribution
preconditioned_logtarget <- function(x) logtarget(as.numeric(preconditioner %*% x))
preconditioned_gradlogtarget <- function(x) as.numeric(gradlogtarget(as.numeric(preconditioner %*% x)) %*% preconditioner)

# initial distribution
rinit <- function() as.numeric(fast_rmvnorm(1, prior_mean, prior_cov))
filename <- paste("inst/logistic/logistic_nsamples", nsamples, "_dimension", dimension, ".RData", sep = "")

save.image(file = filename)

