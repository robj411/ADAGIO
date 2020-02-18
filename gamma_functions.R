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
y <- 1:30
t <- p*y/beta1

hypergeom1F1(t,r,rho)*G(y)*(1-p)^r
plotGraf(function(t)
  hypergeom1F1(1i * t, r, rho), t,
  title = "")

a  = 1 / 2
b  = 1 / 2
t  = seq(-50, 50, length.out = 2 ^ 11)
plotGraf(function(t)
  hypergeom1F1(1i * t, a, b + b), t,
  title = "CF of the Beta distribution with alpha = 1/2, beta = 1/2")


library(fAsianOptions)
y <- 1:30
t <- p*y/beta1

hist(rgamma(1000,shape=incperiod_shape,rate=incperiod_rate)+rgamma(1000,shape=vacc_shape,rate=vacc_rate))
lines(y,(Re(hypergeom1F1(-t,r,rho)*G(y)*(1-p)^r))*2e9)

x = seq(0, 16, length = 200)
plot(x = x, y = kummerM(x, -4.5, 1), type = "l", ylim = c(-25,125),main = "Figure 13.2:  M(-4.5, 1, x)")
lines(x = c(0, 16), y = c(0, 0), col = 2)




mu <- incperiod_shape*incperiod_scale + vacc_shape*vacc_scale
sig2 <- incperiod_shape*incperiod_scale^2 + vacc_shape*vacc_scale^2
alpha <- mu^2/sig2
beta <- sig2/mu

hist(rgamma(1000,shape=incperiod_shape,rate=incperiod_rate)+rgamma(1000,shape=vacc_shape,rate=vacc_rate))
lines(y,dgamma(y,shape=alpha,scale=beta)*5e3)




probability_by_lag_given_removal <- readRDS('probability_by_lag_given_removal_8080.Rds')
probability_after_day_0_given_removal <- readRDS('probability_after_day_0_given_removal_8080.Rds')

probability_by_lag_given_removal_mat <- sapply(probability_by_lag_given_removal,function(x)x[,1])

probability_by_lag_given_removal_mat[1:20,1:10]

probability_by_lag_given_removal_mat <<- sapply(1:20,function(x)pgamma(1:80,shape=incperiod_shape,rate=incperiod_rate)-pgamma((1-x):(80-x),shape=incperiod_shape,rate=incperiod_rate))


recruitment_time <- 30
probability_after_day_0_given_removal <<- lapply(1:20,function(x){
  sapply(1:recruitment_time,function(j){
    poss_inc_val_start <- 1:80
    poss_inc_val_stop <- poss_inc_val_start - x
    denom <- pgamma(poss_inc_val_start,shape=incperiod_shape,rate=incperiod_rate)-pgamma(poss_inc_val_stop,shape=incperiod_shape,rate=incperiod_rate)
    if(j >= recruitment_time){
      poss_inc_val_start <- 1:80
    }else if(recruitment_time > x+j){
      poss_inc_val_start <- 1:80 - x
    }else{
      poss_inc_val_start <- 1:80 + j - recruitment_time
  }
    raw <- pgamma(poss_inc_val_start,shape=incperiod_shape,rate=incperiod_rate)-pgamma(poss_inc_val_stop,shape=incperiod_shape,rate=incperiod_rate)
    raw/denom
  })
  })


