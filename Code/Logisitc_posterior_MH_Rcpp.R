library(LaplacesDemon)
library(invgamma)
library(MASS)
library(tidyr)
library(dplyr)
library(ggplot2)
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
library(rstan)

# setwd("C:/Users/tay15002/Dropbox/Data Science/Final Project/Code")
setwd("C:/Users/Qiangsuper/Dropbox/Data Science/Final Project/Code")
# setwd("~/Dropbox/Data Science/Final Project/Code")

sourceCpp('RCpp_code.cpp')

################################################################################
### 0 - Simulate Data 
################################################################################
set.seed(520)

sim_dat <- function(N){
  d <- matrix(data = NA, nrow = N, ncol = 3)
  
  d[, 1] <- rnorm(n = N, mean = 0, sd = 1)
  d[, 2] <- rnorm(n = N, mean = 0.5, sd = 1)
  d[, 3] <- abs(d[, 1])
  
  y<-rbinom(n = N, size = 1, prob = invlogit(1 + 3 *d[, 1] + 3 *d[, 2] - 3 *d[, 3]))
  
  X<-as.matrix(cbind(1,d[,1:3])) # model matrix
  
  colnames(X) <- c('Intercept', 'X1', 'X2', 'X3')

  return(list('X'=X, 'Y'=y))  
}

d <- sim_dat(N=1000)
X <- d$X
Y <- d$Y

################################################################################
### 1 - functions to sample from conditional posterior distributions
################################################################################

# unnormalized log posterior of beta vector
log_posterior<-function(beta, X, Y){

    # calculate likelihood
  xb <- X %*% beta
  # xb <- ifelse(xb>10, 10, ifelse( xb< (-10) ,-10, xb))
  p_i <- invlogit(xb)
  
  likelihood <- sum(log(dbern(Y, p_i)))
  
  # calculate prior 
  prior <- log(dmvn(x = beta, mu = rep(0,p), Sigma = (10^2)*diag(p)))
  
  log_cond_post <- likelihood + prior
  return(log_cond_post)
}

# Metropolis-Hastings Sampler using log_posterior(),
# which is the log posterior coded in R.
sample_mh<-function(X, Y, iteration, burnin, jump_v){
  
  # create shells
  p <- ncol(X)
  BETA <- matrix(NA, nrow = iteration, ncol = p)
  accept_rate <- numeric(length = iteration)
  
  # starting values
  BETA[1,] <- rep(1, p)
  
  for(i in 2:iteration){
    beta_s <- BETA[i-1, ]
    sd_proposal <- jump_v*diag(p)
    
    # draw from proposal distribution
    beta_star <- mvrnorm(n = 1, beta_s, Sigma = sd_proposal)
    
    # calculate ratio of conditional posterior densities
    r_num <- log_posterior(beta_star, X, Y ) 
              + log(dmvn(x = beta_s, mu = beta_star, Sigma = sd_proposal))
    
    r_denom <- log_posterior(beta_s, X, Y ) 
              + log(dmvn(x = beta_star, mu = beta_s, Sigma = sd_proposal))
    
    # calculate acceptance probability
    r <- exp(r_num - r_denom)
    rmin<-min(r,1)
    
    # accept or reject proposal
    if( rmin > runif(n = 1, min = 0, max = 1) ){ 
      BETA[i, ] <- beta_star
    }else{
      BETA[i, ] <- beta_s
    }
    accept_rate[i] <- rmin
    
  }
  
  colnames(BETA) <- colnames(X)
  colnames(BETA)[1] <- 'intercept'
  
  return(list('BETA' = BETA, 'accept_rate' = accept_rate) )
}

# Metropolis-Hastings Sampler using log_post(),
# which is the log posterior coded in C++.

sample_mh_cpp <-function(X, Y, iteration, jump_v){
  # create shells
  p <- ncol(X)
  BETA <- matrix(NA, nrow = iteration, ncol = p)
  accept_rate <- numeric(length = iteration)
  
  # starting values
  BETA[1,] <- rep(1, p)
  
  for(i in 2:iteration){
    beta_s <- BETA[i-1, ]
    
    # draw from proposal distribution
    beta_star <- mvrnorm(n = 1, beta_s, Sigma = sqrt(jump_v*diag(p)))
    
    # calculate ratio of conditional posterior densities
    r_num <- log_post(beta_star, Y, X, beta_s, jump_v)
    r_denom <- log_post(beta_s, Y, X, beta_star, jump_v)
    
    # calculate acceptance probability
    r <- exp(r_num - r_denom)
    rmin<-min(r,1)
    
    # accept or reject proposal
    if( rmin > runif(n = 1, min = 0, max = 1) ){ 
      BETA[i, ] <- beta_star
    }else{
      BETA[i, ] <- beta_s
    }
    accept_rate[i] <- rmin
    
  }
  colnames(BETA) <- colnames(X)
  colnames(BETA)[1] <- 'intercept'
  return(list('BETA' = BETA, 'accept_rate' = accept_rate) )
  
  
}

plot.together <- function(result) {
  par(mfrow=c(2,2))
  plot(result$BETA[, 1], type='l',
       xlab='MH Iteration', ylab='Posterior Draw', main='Intercept')

  plot(result$BETA[, 2], type='l',
       xlab='MH Iteration', ylab='Posterior Draw', main='beta_1')

  plot(result$BETA[, 3], type='l',
       xlab='MH Iteration', ylab='Posterior Draw', main='beta_2')

  plot(result$BETA[, 4], type='l',
       xlab='MH Iteration', ylab='Posterior Draw', main='beta_3')

  par(mfrow=c(1,1))
}

################################################################################
### 2 - Test the Samplers
################################################################################

burnin <- 25000
iteration <- 50000
jump_v <- 0.01
p <- ncol(X)

# pure R
res_mh <- sample_mh(X = X, Y = Y, iteration = iteration, burnin = burnin, jump_v = jump_v)

plot.together(result = res_mh)

# RCpp
res_mh_cpp <- sample_mh_cpp(X, Y, iteration = iteration, jump_v = jump_v)

plot.together(result = res_mh_cpp)

# pure Cpp
res_pure_cpp <- for_log_post(X = X, Y = Y, iteration = iteration, jump_v = jump_v)

plot.together(result = res_pure_cpp)

# RStan
res_rstan <- stan(file = 'posterior_in_stan.stan', 
                  data = list(x = X, y = Y, N = nrow(X), p = ncol(Y)), iter = iteration, chains = 1)

res_rstan <- as.data.frame(res_rstan)

par(mfrow=c(2,2))
plot(res_rstan$`beta[1]`, type='l',
     xlab='MH Iteration', ylab='Posterior Draw', main='Intercept')

plot(res_rstan$`beta[2]`, type='l',
     xlab='MH Iteration', ylab='Posterior Draw', main='beta_1')

plot(res_rstan$`beta[3]`, type='l',
     xlab='MH Iteration', ylab='Posterior Draw', main='beta_2')

plot(res_rstan$`beta[4]`, type='l',
     xlab='MH Iteration', ylab='Posterior Draw', main='beta_3')

par(mfrow=c(1,1))

################################################################################
### 3 - Benchmarks
################################################################################
iteration <- 20000

ss <- seq(from = 50, by = 50, to = 500)

rel_time <- matrix(data = NA, nrow = length(ss), ncol = 3)

r_time <- numeric(length = length(ss))
rcpp_time <- numeric(length = length(ss))
cpp_time <- numeric(length = length(ss))
stan_time <- numeric(length = length(ss))


for(i in 1:length(ss) ){
  d <- sim_dat(N = ss[i])
  X <- d$X
  Y <- d$Y
  
  bench<-microbenchmark(R_MH = sample_mh(X, Y, iteration = iteration, jump_v = jump_v),
                        Cpp_MH = sample_mh_cpp(X, Y, iteration = iteration, jump_v = jump_v),
                        Pure_cpp = for_log_post(X = X, Y = Y, iteration = iteration, jump_v = jump_v),
                        Stan = stan(file = 'posterior_in_stan.stan', 
                                    data = list(x = X, y = Y, N = nrow(X), p = ncol(Y)), iter = iteration, chains = 1),
                        times = 2)
  bench_sum <- summary(bench)
  r_time[i] <- bench_sum$mean[bench_sum$expr=='R_MH']
  rcpp_time[i] <- bench_sum$mean[bench_sum$expr=='Cpp_MH']
  cpp_time[i] <- bench_sum$mean[bench_sum$expr=='Pure_cpp']
  stan_time[i] <- bench_sum$mean[bench_sum$expr=='Stan']
  rel_time[i, 1] <- r_time[i]/rcpp_time[i]
  rel_time[i, 2] <- r_time[i]/cpp_time[i]
  rel_time[i, 3] <- r_time[i]/stan_time[i]
}

time_summary <- data.frame(cbind(ss, r_time, rcpp_time, cpp_time, stan_time))
colnames(time_summary) <- c('Sample Size', 'R', 'RCpp', 'Cpp', 'Stan')

rel_time_summary <- data.frame(cbind(ss, rel_time))
colnames(rel_time_summary) <- c('Sample Size', 'R/Rcpp', 'R/Cpp', 'R/Stan')

ggplot(time_summary, aes(time_summary$`Sample Size`))+
  geom_line(aes(y = time_summary$`R`, colour = 'R'))+ 
  geom_line(aes(y = time_summary$`RCpp`, colour = 'RCpp'))+
  geom_line(aes(y = time_summary$`Stan`, colour = 'Stan'))+
  geom_line(aes(y = time_summary$`Cpp`, colour = 'Cpp'))+
  xlab('Data Sample Size')+ 
  ylab('Runtime Time')+
  ggtitle('Running Time Summary')

ggplot(rel_time_summary, aes(rel_time_summary$`Sample Size`))+
  geom_line(aes(y = rel_time_summary$`R/Rcpp`, colour = 'R/Rcpp'))+ 
  geom_line(aes(y = rel_time_summary$`R/Cpp`, colour = 'R/Cpp'))+
  geom_line(aes(y = rel_time_summary$`R/Stan`, colour = 'R/Stan'))+
  xlab('Data Sample Size')+ 
  ylab('Relative Runtime Time (R as baseline)')+
  ggtitle('Relative Running Time against Sample Size')

d <- sim_dat(N = 50)
X <- d$X
Y <- d$Y

benchmark_result <- microbenchmark(
  sample_mh(X, Y, iteration = iteration, jump_v = jump_v),
  sample_mh_cpp(X, Y, iteration = iteration, jump_v = jump_v),
  for_log_post(X = X, Y = Y, iteration = iteration, jump_v = jump_v),
  stan(file = 'posterior_in_stan.stan', 
       data = list(x = X, y = Y, N = nrow(X), p = ncol(Y)), iter = iteration, chains = 1)
)

autoplot(benchmark_result) + scale_x_discrete(labels = c('Pure R', 'RCpp', 'Pure Cpp', 'Stan'))

save.image(file = "project.RData")
