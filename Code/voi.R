{rm(list=ls())
  library(Rlab)
  library(lme4)
  library(coxme)
  library(Matrix)
  library(igraph)
  library(deSolve)
  library(NetSurv)
  library(ggplot2)
  library(caTools)
  library(reshape2)
  library(survival)
  library(frailtypack)
  library(rootSolve)
  library(fields)
  library(foreach)
  library(doParallel)
  library(rethinking)
  library(latex2exp)
  setwd('~/overflow_dropbox/ADAGIO/Code')
  source('functions.R')
}
# August 29th, 2017
# Supplementary Material to accompany "Competing effects of indirect 
# protection and clustering on the power of cluster-randomized 
# controlled vaccine trials"

# Script to simulate an individually-randomized and cluster-randomized
# trial a given number of times, and return vaccine effect estimate and
# estimated power for a given set of parameters
# Also contains code to run the sample size calculations based on final
# size equations
{
  # Number of simulated trials
  nsim <- 100
  
  # Population structure parameters:
  # Average size of one community
  ave_community_size <- 100
  # Range of community sizes (sizes are uniformly distributed on this range)
  community_size_range <- 40
  # Number of communities
  num_communities <- 80
  # Probability of an edge between two nodes in the same community
  rate_within <- 0.15
  # Probability of an edge between two nodes in different communities
  rate_between <- 0
  
  # Disease characteristics:
  # Per-time-step hazard of infection for a susceptible nodes from an infectious
  # neighbour
  beta <- 0.01
  # Expected number of importations to the population over two years
  num_introductions <- 50
  # Leaky multiplicative efficacy of vaccine
  direct_VE <- 0.6
  
  disease_dynamics <- list(beta=beta,
                           num_introductions=num_introductions)
  
  trial_design <- list_trial_parameters(bTrial=2,
                                               adaptation_day = 40,
                                               adaptation='Ney') # iRingNey-40
  
  
  ## get infection trajectory for source population
  num_timesteps <- trial_design$trial_startday + trial_design$trial_length + trial_design$enrollment_period - 1
  times <- seq(0,num_timesteps,1)
  infected_trajectory <- get_infected_trajectory(times)
  
  ## set up for results
  registerDoParallel(cores=18)
  simnum <- sday <- tr <- 1
}

# Gamma-distribution parameters of incubation and infectious period
incperiod_shape <- 3.11
incperiod_rate <- 0.32
infperiod_shape <- 1.13
infperiod_rate <- 0.226
ave_inc_period <- incperiod_shape/incperiod_rate
ave_inf_period <- infperiod_shape/infperiod_rate
nsim <- 100
{par(mfrow=c(2,2),mar=c(5,4,4,2))
hist(rgamma(nsim,shape=3.11,rate=0.32))
hist(rgamma(nsim,shape=1.13,rate=0.226))
inc_shape <- rgamma(nsim,shape=ave_inc_period*2,rate=2)
inf_shape <- rgamma(nsim,shape=ave_inf_period*2,rate=2)
hist(rlnorm(nsim,log(inc_shape),0.2))
hist(rlnorm(nsim,log(inf_shape),0.2))
}
x11(); plot(0:30,dgamma(0:30,shape=inc_shape[1]*0.311,rate=0.311)); for(i in 2:10) lines(0:30,dgamma(0:30,shape=inc_shape[i]*0.311,rate=0.311))
x11(); plot(0:30,dgamma(0:30,shape=inf_shape[1]*0.311,rate=0.311)); for(i in 2:10) lines(0:30,dgamma(0:30,shape=inf_shape[i]*0.311,rate=0.311))
pvals <- c()
ves <-c()
for (simnum in 1:nsim) {
  
  disease_dynamics <- list(beta=beta,
                           num_introductions=num_introductions,
                           incperiod_shape=inc_shape[simnum]*incperiod_rate,
                           incperiod_rate=incperiod_rate,
                           infperiod_shape=inf_shape[simnum]*infperiod_rate,
                           infperiod_rate=infperiod_rate,
                           ave_inc_period=inc_shape[simnum])
  
  g <<- make_network(ave_community_size, community_size_range, num_communities,rate_within, rate_between)
  
  #trial_outcomes <- foreach(tr = trial_indicies) %dopar% {
  out <- core_trial_script(trial_design=trial_design,g)
  pvals[simnum] <- out$pval
  ves[simnum] <- out$VaccineEfficacy
  #}
  cat("Simulation ",simnum,"\n")
}

plot(inc_shape[1:nsim],pvals)
plot(inc_shape[1:nsim],ves)
plot(inf_shape[1:nsim],pvals)
plot(inf_shape[1:nsim],ves)
