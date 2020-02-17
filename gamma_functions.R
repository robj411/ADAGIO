incperiod_shape <- 3.11
incperiod_rate <- 0.32

vacc_shape <- 1
vacc_rate <- 3

hist(rgamma(1000,shape=vacc_shape,rate=vacc_rate))

hist(rgamma(1000,shape=incperiod_shape,rate=incperiod_rate)+rgamma(1000,shape=vacc_shape,rate=vacc_rate))

incperiod_scale <- 1/incperiod_rate
vacc_scale <- 1/vacc_rate

beta1 <- vacc_scale
beta2 <- incperiod_scale
p <- 1 - beta1/beta2

r <- incperiod_shape

rho <- incperiod_shape + vacc_shape

G <- function(x)dgamma(x,rho,beta1)

library(CharFun)