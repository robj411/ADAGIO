rm(list=ls())
setwd('~/overflow_dropbox/ADAGIO/COVID19/')
#source('../Code/functions_network.R')
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
library(readxl)
library(rBeta2009)
set.seed(0)
## build network ###############################################################

source('build_network.R')
eligible_first_person <- sapply(contact_of_contact_list,length)>10

## functions #######################################################

source('../network_functions.R')
source('../evaluation_functions.R')

ref_recruit_day <<- 30
observed <<- 0.8
eval_day <<- 25
target_weight <<- 24

covid_spread_wrapper <- function(i_nodes_info,s_nodes,v_nodes,e_nodes_info,direct_VE){
  # to contacts
  # e infects house and work and anyone - only enodes infected one day ago or more, and only enodes with one day left incubating
  ##!! a subset of i_nodes are nonsymptomatic and therefore continue to infect contacts. these should be a fixed list, not sampled randomly every time.
  current_infectious <- c(i_nodes_info[c(runif(nrow(i_nodes_info))<observed),1],e_nodes_info[e_nodes_info[,2]>=e_nodes_info[,3],1])
  if(length(current_infectious)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_list,beta_scalar=nonrandom_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=random_list,beta_scalar=random_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
  }
  # i infects house
  current_infectious <- c(i_nodes_info[,1])
  if(length(current_infectious)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=household_list,beta_scalar=1)
  }
  return(e_nodes_info)
}

## set up #######################################################

# Per-time-step hazard of infection for a susceptible nodes from an infectious
# neighbour
beta_base <<- 0.01
# Gamma-distribution parameters of incubation and infectious period and wait times
# hist(rgamma(1000,shape=1.43,rate=0.549)+2)
infperiod_shape <<- 1.43
infperiod_rate <<- 0.549
infperiod_const <<- 2
## assume there is no difference between infectious with and without symptoms - all I
#mn <- 6.5 # =shape/rate # 5.2
#sd <- 2.6 # =sqrt(shape/rate^2) # 2.8
# hist(rgamma(1000,shape=13.3,rate=4.16)+2)
incperiod_rate <<- 4.16 # mn/sd^2 # 
incperiod_shape <<- 13.3 # incperiod_rate*mn # 
incperiod_const <<- 2
#hosp_shape_index <<- 2
#hosp_rate_index <<- 0.5
#hosp_shape <<- 2
#hosp_rate <<- 0.5
recruit_shape <<- 5.4
recruit_rate <<- 0.47
#hosp_mean_index <<- 3.85
#hosp_sd_index <<- 2.76
#hosp_mean <<- 3.85
#hosp_sd <<- 2.76
recruit_mean <<- 10.32
recruit_sd <<- 4.79
enrollment_rate <<- 0.7
nonrandom_scalar <<- 1
random_scalar <<- 1/10
length(E(new_g))
length(E(random_g))
direct_VE <- 0.0

enrolled_per_contact <<- enrollment_rate*mean(sapply(contact_of_contact_list,length)[eligible_first_person])

g <<- new_g

g_name <<- as.numeric(as.vector(V(g)$name))
vertices <- V(g)
cluster_size <- hosp_times <- recruit_times <- c()
results_list <- list()


###########################################################################

set_variables_from_gamma_distributions <- function(){
  vacc_shape <<- 3
  vacc_rate <<- 1
  
  infperiod_scale <- 1/infperiod_rate
  incperiod_scale <- 1/incperiod_rate
  vacc_scale <- 1/vacc_rate
  
  mu <- incperiod_shape*incperiod_scale + vacc_shape*vacc_scale
  sig2 <- incperiod_shape*incperiod_scale^2 + vacc_shape*vacc_scale^2
  alpha <- mu^2/sig2
  beta <- sig2/mu
  inc_plus_vacc_shape <<- alpha
  inc_plus_vacc_rate <<- 1/beta
  
  zero <<- 50
  pgamma_vector <<- pgamma(1:100-incperiod_const,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
  dgamma_vector <<- dgamma(1:100-incperiod_const,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
  pgamma_inc_vector <<- c(rep(0,zero),pgamma(1:100-incperiod_const,shape=incperiod_shape,rate=incperiod_rate))
  pgamma_vacc_vector <<- c(rep(0,zero),pgamma(1:100,shape=vacc_shape,rate=vacc_rate))
  dgamma_inc_vector <<- c(rep(0,zero),dgamma(1:100-incperiod_const,shape=incperiod_shape,rate=incperiod_rate))
  
  recruitment_time <<- 30
}
set_variables_from_gamma_distributions()



