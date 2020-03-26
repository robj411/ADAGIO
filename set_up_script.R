rm(list=ls())
setwd('~/overflow_dropbox/ADAGIO/')
source('Code/functions_network.R')
library(igraph)
library(truncnorm)
library(infotheo)
library(xtable)
library(RColorBrewer)
library(plotrix)
library(profvis)
library(funique)
library(doParallel)
library(foreach)
library(survival)
library(coxme)
library(pracma)
library(dplyr)
library(rBeta2009)

## build network ###############################################################

source('build_network.R')

## functions #######################################################

source('network_functions.R')
source('evaluation_functions.R')

probability_by_lag <<- readRDS(paste0('storage/probability_by_lag_8080.Rds'))
probability_after_day_0 <<- readRDS(paste0('storage/probability_after_day_0_8080.Rds'))
#probability_by_lag_given_removal <<- readRDS(paste0('probability_by_lag_given_removal_8080.Rds'))
#probability_after_day_0_given_removal <<- readRDS(paste0('probability_after_day_0_given_removal_8080.Rds'))
ref_recruit_day <<- 30

## set up #######################################################

# Per-time-step hazard of infection for a susceptible nodes from an infectious
# neighbour
beta_base <<- 0.0065
high_risk_scalar <<- 2.17
# fraction of beta_base applied to neighbours ("contacts of contacts")
neighbour_scalar <<- 0.39
# Gamma-distribution parameters of incubation and infectious period and wait times
incperiod_shape <<- 3.11
incperiod_rate <<- 0.32
infperiod_shape <<- 1.13
infperiod_rate <<- 0.226
hosp_shape_index <<- 2
hosp_rate_index <<- 0.5
hosp_shape <<- 2
hosp_rate <<- 1.5
recruit_shape <<- 5.4
recruit_rate <<- 0.47
hosp_mean_index <<- 3.85
hosp_sd_index <<- 2.76
hosp_mean <<- 2.8
hosp_sd <<- 1.5
vacc_mean <<- 1.5
vacc_sd <<- 1
recruit_mean <<- 10.32
recruit_sd <<- 4.79
enrollment_rate <<- 0.5
direct_VE <- 0.0

g <<- new_g

g_name <<- as.numeric(as.vector(V(g)$name))
vertices <- V(g)
cluster_size <- hosp_times <- recruit_times <- c()
results_list <- list()


###########################################################################

gamma_integral <- function(z,s2,recruitment_time) {
  sapply(z,function(xp) {
    dgamma(xp,shape=incperiod_shape,rate=incperiod_rate)*
      pgamma(s2-recruitment_time-xp,shape=vacc_shape,rate=vacc_rate)
  })
}

set_variables_from_gamma_distributions <- function(){
  vacc_shape <<- 1
  vacc_rate <<- 1
  
  incperiod_scale <- 1/incperiod_rate
  vacc_scale <- 1/vacc_rate
  
  mu <- incperiod_shape*incperiod_scale + vacc_shape*vacc_scale
  sig2 <- incperiod_shape*incperiod_scale^2 + vacc_shape*vacc_scale^2
  alpha <- mu^2/sig2
  beta <- sig2/mu
  
  inc_plus_vacc_shape <<- alpha
  inc_plus_vacc_rate <<- 1/beta
  zero <<- 50
  pgamma_vector <<- pgamma(1:100,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
  dgamma_vector <<- dgamma(1:100,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
  pgamma_inc_vector <<- c(rep(0,zero),pgamma(1:100,shape=incperiod_shape,rate=incperiod_rate))
  pgamma_vacc_vector <<- c(rep(0,zero),pgamma(1:100,shape=vacc_shape,rate=vacc_rate))
  dgamma_inc_vector <<- c(rep(0,zero),dgamma(1:100,shape=incperiod_shape,rate=incperiod_rate))
  
  probability_by_lag_given_removal_mat <<- sapply(1:20,function(x)pgamma_inc_vector[zero+1:80]-pgamma_inc_vector[zero+(1-x):(80-x)])
  
  
  recruitment_time <<- 30
  #system.time(
  probability_after_day_0_given_removal <<- lapply(1:20,function(x){
    sapply(1:80,function(j){
      poss_inc_val_start <- 1:80
      poss_inc_val_stop <- poss_inc_val_start - x
      denom <- pgamma_inc_vector[zero+poss_inc_val_start]-pgamma_inc_vector[zero+poss_inc_val_stop]
      subtract <- max(0, min(recruitment_time - j, x) )
      #if(j >= recruitment_time){
      #  poss_inc_val_start <- 1:80
      #}else if(recruitment_time > x+j){
      #  poss_inc_val_start <- 1:80 - x
      #}else{
      #  poss_inc_val_start <- 1:80 + j - recruitment_time
      #}
      #num <- pgamma(poss_inc_val_start,shape=incperiod_shape,rate=incperiod_rate)-pgamma(poss_inc_val_stop,shape=incperiod_shape,rate=incperiod_rate)
      num <- sapply(poss_inc_val_start,function(ii)integrate(gamma_integral , ii-x, ii - subtract,s2=ii+j,recruitment_time=recruitment_time)$value)
      pmin(num/denom,1)
    })
  }
  )
  #)
}
set_variables_from_gamma_distributions()

x <- 20
j <- 43
i <- 3
s1 <- j
s2 <- i+j
s1p <- x+s1
subtract <- max(0, min(recruitment_time - j, x) )
inc_period <- rgamma(10000,shape=incperiod_shape,rate=incperiod_rate)
vacc_period <- rgamma(10000,shape=vacc_shape,rate=vacc_rate)
denom <- sapply(i,function(ii)sum(inc_period < ii & inc_period > ii-x))
num <- sapply(i,function(ii)sum(inc_period < ii - subtract - vacc_period & inc_period > ii-x))
c(num,denom)/10000
print(num/denom)
denom <- pgamma(i,shape=incperiod_shape,rate=incperiod_rate)-pgamma(i-x,shape=incperiod_shape,rate=incperiod_rate)
# int from 0 to x with dgamma
num <- sapply(i,function(ii)integrate(gamma_integral , ii-x, ii - subtract,ii+j,recruitment_time)$value)
c(num,denom)
print(num/denom)

s1 <- 20
s1p <- 31
s2 <- 35
i <- s2 - s1
j <- s1
x <- s1p - s1
inc_period <- rgamma(10000,shape=incperiod_shape,rate=incperiod_rate)
denom <- sapply(i,function(ii)sum(inc_period < ii & inc_period > ii-x))
subtract <- max(0, min(recruitment_time - j, x) )
num <- sapply(i,function(ii)sum(inc_period < ii - subtract & inc_period > ii-x))
print(num/denom)
denom <- pgamma(i,shape=incperiod_shape,rate=incperiod_rate)-pgamma(i-x,shape=incperiod_shape,rate=incperiod_rate)
num <- pgamma(i - subtract,shape=incperiod_shape,rate=incperiod_rate)-pgamma(i-x,shape=incperiod_shape,rate=incperiod_rate)
print(num/denom)


