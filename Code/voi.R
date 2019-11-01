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
  library(dplyr)
  library(Rcpp)
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
  ave_community_size <- 1000
  # Range of community sizes (sizes are uniformly distributed on this range)
  community_size_range <- 40
  # Number of communities
  num_communities <- 1
  # Probability of an edge between two nodes in the same community
  rate_within <- 0.015
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

extF <- -log(1-num_introductions/(ave_community_size*num_communities))/trapz(times,infected_trajectory) * infected_trajectory
bg <- 1 - exp(-extF)  # 1 - exp(-(1-direct_VE)*extF)
prethreshold <- function(x,day1) {dgamma(day1-x,shape=incperiod_shape,rate=incperiod_rate)*approx(x=0:(trial_length + trial_startday),y=bg,xout=x)$y }
postthreshold <- function(x,day0,day1) {
  dgamma(day1-x,shape=incperiod_shape,rate=incperiod_rate)*
    (approx(x=0:(trial_length + trial_startday),y=bg,xout=x)$y+(1-exp(-beta))*dgamma(x-day0,shape=infperiod_shape,rate=infperiod_rate)) }
likelihood <- function(x,day0,day1) {
  if(x<day0){
    dgamma(day1-x,shape=incperiod_shape,rate=incperiod_rate)*approx(x=0:(trial_length + trial_startday),y=bg,xout=x)$y
  }else{
    dgamma(day1-x,shape=incperiod_shape,rate=incperiod_rate)*
    (approx(x=0:(trial_length + trial_startday),y=bg,xout=x)$y+(1-exp(-beta))*dgamma(x-day0,shape=infperiod_shape,rate=infperiod_rate)) 
  }
}

nsim <- 1000
pvals <- matrix(0,nrow=nsim,ncol=2)
ves <-matrix(0,nrow=nsim,ncol=2)
results_list <- list()
for(i in 1:2){
  direct_VE <- c(0,0.6)[i]
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
    
    g_name <<- V(g)$name
    g_matrix <<- as.matrix(igraph::as_adjacency_matrix(g))
    g_community <<- V(g)$community
    
    #profvis(
    list[results,trial_nodes,trajectories,allocation_rates]<-
      network_epidemic(disease_dynamics,direct_VE,infected_trajectory,trial_design)
    #)
    cat("Simulation ",i,simnum,"\n")
    
    results_to_join <- results
    colnames(results_to_join)[1] <- 'Node'
    trial_results <- left_join(trial_nodes,results_to_join[,c(1,2)],by='Node')
    results_list[[simnum]] <- trial_results
  }
  saveRDS(results_list,paste0('1000simulations',direct_VE,'.Rds'))
}

  
pvals <- list()
ves <-list()
for(fu in 1:3){
  pvals[[fu]] <- list()
  ves[[fu]] <-list()
  trial_design$follow_up <- c(21,40,400)[fu]
  for(i in 1:2){
    pvals[[fu]][[i]] <- matrix(0,nrow=nsim,ncol=5)
    ves[[fu]][[i]] <-matrix(0,nrow=nsim,ncol=5)
    direct_VE <- c(0,0.6)[i]
    results_list <- readRDS(paste0('1000simulations',direct_VE,'.Rds'))
    for(simnum in 1:length(results_list)){
      
      trial_results <- results_list[[simnum]] 
      
      trial_results$infected <- !is.na(trial_results$DayInfected) & trial_results$DayInfected - trial_results$DayVaccinated < trial_design$follow_up
      trial_results$DayInfected[is.na(trial_results$DayInfected)] <- trial_length + trial_startday
      trial_results$threshold <- trial_results$DayInfected - trial_results$DayVaccinated
      trial_results$weight <- 1
      if(sum(trial_results$infected)==0){
        pvals[[fu]][[i]][simnum,] <- 1
        ves[[fu]][[i]][simnum,] <- 0
      }else{
        for(j in 1:3){
          if(j==1){
            trial_results$weight[trial_results$threshold<10] <- 0
          }else if (j==2){
            trial_results$weight <- pgamma(trial_results$threshold,shape=incperiod_shape,rate=incperiod_rate)
          }else if(j==3){
            trial_results$weight <- pgamma(trial_results$threshold,shape=incperiod_shape,rate=incperiod_rate)
            pos_results <- which(trial_results$infected == T)
            pos_set <- trial_results[pos_results,]
            prob_before <- apply(pos_set,1,function(x) integrate(prethreshold,x[6],lower=0,upper=x[5])$value)
            prob_after <- apply(pos_set,1,function(x) integrate(postthreshold,x[5],x[6],lower=x[5],upper=x[6])$value)
            pos_set$weight2 <- prob_after/(prob_after+prob_before)
            trial_results$weight[pos_results] <- prob_after/(prob_after+prob_before)
          }
          
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
          
          pvals[[fu]][[i]][simnum,j] <- pval_binary_mle
          ves[[fu]][[i]][simnum,j] <- VE_pointest_binary_mle
          if(j>1){
            survmodel<-coxph(Surv(DayVaccinated,DayInfected, infected) ~ TrialStatus,#+tt(DayVaccinated),
                             data=trial_results,
                             weights=weight,
                             tt=function(x,t,...) dgamma(t-x,shape=infperiod_shape,rate=infperiod_rate))
            vaccEffEst <- 1-exp(survmodel$coefficient[1] + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var[1])))
            zval <- survmodel$coefficient[1]/sqrt(survmodel$var[1])
            pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
            
            pvals[[fu]][[i]][simnum,j+2] <- pval
            ves[[fu]][[i]][simnum,j+2] <- vaccEffEst[1]
          }
        }
      }
    }
  }
}
lapply(pvals,function(y)sapply(y,function(z) apply(z,2,function(x)sum(x<0.05,na.rm=T))))
apply(ves,2,function(x)c(mean(x,na.rm=T),sd(x,na.rm=T)))
col_labels <- c('Exclusion binary','Incubation-weighted binary','Infectious-weighted binary','Incubation-weighted TTE','Infectious-weighted TTE')
for(fu in 1:3){
  for(i in 1:5){
    cat(paste0(col_labels[i],' & ',c(21,40,400)[fu],' & ',round(sum(pvals[[fu]][[2]][,i]<0.05,na.rm=T)/nsim,3),
               ' & ',round(sum(pvals[[fu]][[1]][,i]<0.05,na.rm=T)/nsim,3),
               ' & ',round(mean(ves[[fu]][[2]][,i],na.rm=T),3),' (',round(sd(ves[[fu]][[2]][,i],na.rm=T),3),') \\\\ \n'))
  }
  cat('\\hline \n')
}

day0 <- 99
day1 <- 100
pb <- integrate(prethreshold,day1,lower=0,upper=day0)$value
pa <- integrate(postthreshold,day0,day1,lower=day0,upper=day1)$value
c(pb,pa)/(pb+pa)


plot(inc_shape[1:nsim],pvals)
plot(inc_shape[1:nsim],ves)
plot(inf_shape[1:nsim],pvals)
plot(inf_shape[1:nsim],ves)
