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

## build network ###############################################################

source('build_network.R')

## functions #######################################################

source('network_functions.R')
source('evaluation_functions.R')

probability_by_lag <<- readRDS(paste0('probability_by_lag_8080.Rds'))
probability_after_day_0 <<- readRDS(paste0('probability_after_day_0_8080.Rds'))
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
direct_VE <- 0.0

g <<- new_g

g_name <<- as.numeric(as.vector(V(g)$name))
vertices <- V(g)
cluster_size <- hosp_times <- recruit_times <- c()
results_list <- list()


###########################################################################

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


