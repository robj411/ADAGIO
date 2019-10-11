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
                                               follow_up = 40) # iRingFR-40
  
  
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

nsim <- 1000
pvals <- matrix(0,nrow=nsim,ncol=2)
ves <-matrix(0,nrow=nsim,ncol=2)
for (simnum in 1:nsim) {
  
  #disease_dynamics <- list(beta=beta,
  #                         num_introductions=num_introductions,
  #                         incperiod_shape=inc_shape[simnum]*incperiod_rate,
  #                         incperiod_rate=incperiod_rate,
  #                         infperiod_shape=inf_shape[simnum]*infperiod_rate,
  #                         infperiod_rate=infperiod_rate,
  #                         ave_inc_period=inc_shape[simnum])
  disease_dynamics <- list(beta=beta,
                           num_introductions=num_introductions,
                           incperiod_shape=incperiod_shape,
                           incperiod_rate=incperiod_rate,
                           infperiod_shape=infperiod_shape,
                           infperiod_rate=infperiod_rate,
                           ave_inc_period=ave_inc_period)
  g <<- make_network(ave_community_size, community_size_range, num_communities,rate_within, rate_between)
  
  #trial_outcomes <- foreach(tr = trial_indicies) %dopar% {
  trial_startday <- trial_design$trial_startday#100#50 + (sday-1)*100
  trial_length <- trial_design$trial_length#500 - trial_startday
  bTrial <- trial_design$bTrial
  bCluster <- trial_design$bCluster
  #adaptation <- trial_designs[[tr]]$adaptation
  #vaccination_gap <- trial_designs[[tr]]$vaccination_gap
  follow_up <- trial_design$follow_up
  #adaptation_day <- trial_designs[[tr]]$adaptation_day
  #profvis(
  list[results,trial_nodes,trajectories,allocation_rates]<-
    network_epidemic(g,disease_dynamics,direct_VE,infected_trajectory,trial_design)
  #)
  
  results_to_join <- results
  colnames(results_to_join)[1] <- 'Node'
  trial_results <- left_join(trial_nodes,results_to_join[,c(1,2)],by='Node')
  trial_results$infected <- !is.na(trial_results$DayInfected) & trial_results$DayInfected - trial_results$DayVaccinated < trial_design$follow_up
  trial_results$DayInfected[is.na(trial_results$DayInfected)] <- trial_length + trial_startday
  trial_results$incubation <- trial_results$DayInfected - trial_results$DayVaccinated
  trial_results$weight <- pgamma(trial_results$incubation,shape=incperiod_shape,rate=incperiod_rate)
  
  fail0 <- sum(subset(trial_results,TrialStatus==0&infected==T)$weight)
  fail1 <- sum(subset(trial_results,TrialStatus==1&infected==T)$weight)
  success0 <- sum(subset(trial_results,TrialStatus==0&infected==F)$weight)
  success1 <- sum(subset(trial_results,TrialStatus==1&infected==F)$weight)
  n0 <- fail0 + success0
  n1 <- fail1 + success1
  
  p0 <- success0/n0
  p1 <- success1/n1
  sigma0 <- p0 * ( 1 - p0 ) /n0
  sigma1 <- p1 * ( 1 - p1 ) /n1
  zval <- (p1-p0)/(sqrt(sigma0+sigma1))
  pval_binary_mle <- dnorm(zval)
  
  VE_pointest_binary_mle <- 1 - (fail1/n1)/(fail0/n0)
  
  list[VE,pval,events_vacc,events_cont,analysed_trialsize] <- 
    analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,bCluster,follow_up,trial_design$revisit)
  pvals[simnum,] <- c(pval,pval_binary_mle)
  ves[simnum,] <- c(VE,VE_pointest_binary_mle)
  #}
  cat("Simulation ",simnum,"\n")
}
apply(pvals,2,function(x)sum(x<0.05))
apply(ves,2,function(x)c(mean(x),sd(x)))

plot(inc_shape[1:nsim],pvals)
plot(inc_shape[1:nsim],ves)
plot(inf_shape[1:nsim],pvals)
plot(inf_shape[1:nsim],ves)
