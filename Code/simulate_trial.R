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
# Gamma-distribution parameters of incubation and infectious period
incperiod_shape <- 3.11
incperiod_rate <- 0.32
infperiod_shape <- 1.13
infperiod_rate <- 0.226
ave_inc_period <- ceiling(incperiod_shape/incperiod_rate)

disease_dynamics <- list(beta=beta,
                         num_introductions=num_introductions,
                         incperiod_shape=incperiod_shape,
                         incperiod_rate=incperiod_rate,
                         infperiod_shape=infperiod_shape,
                         infperiod_rate=infperiod_rate,
                         ave_inc_period=ave_inc_period)

# Calculate R0
R0 <- (1 - (infperiod_rate/(infperiod_rate+beta))^infperiod_shape) *
  (((ave_community_size-1)*(1-rate_within)*rate_within + 
      (num_communities-1)*ave_community_size*(1-rate_between)*rate_between + 
      ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)^2)/
     ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)-1)
cat("R0: ",R0,"\n")
pdf('R0.pdf');
par(mar=c(5,5,2,2));plot(seq(0.0025,0.4,by=0.0025),(1 - (infperiod_rate/(infperiod_rate+seq(0.0025,0.4,by=0.0025)))^infperiod_shape) *
                           (((ave_community_size-1)*(1-rate_within)*rate_within + 
                               (num_communities-1)*ave_community_size*(1-rate_between)*rate_between + 
                               ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)^2)/
                              ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)-1),
                         ylab=TeX('$R_0$'),xlab=TeX('$\\beta$'),frame=F,typ='l',lwd=2,col='navyblue',cex.axis=1.5,cex.lab=1.5)
abline(h=1,lwd=2,lty=2,col='grey')
abline(v=0.01,lwd=2,lty=2,col='grey')
text(x=0.01,y=8,labels=TeX('$\\beta$=0.01'),pos=4,cex=1.5,col='navyblue')
text(x=0.3,y=1,labels=TeX('$R_0$=1'),pos=3,cex=1.5,col='navyblue')
dev.off() 

# Initialize vectors/data frames to store results
pvals<-data.frame('iRCT'=rep(NA,nsim),
                  'DiRCT'=rep(NA,nsim),
                  'FAiRCT'=rep(NA,nsim),
                  'TSiRCT'=rep(NA,nsim),
                  'cRCT_gaussian_coxme'=rep(NA,nsim),
                  'cRCT_gaussian_coxph'=rep(NA,nsim),
                  'cRCT_gamma_coxph'=rep(NA,nsim),
                  'cRCT_gee'=rep(NA,nsim))
VaccineEfficacy<-data.frame('iRCT'=rep(NA,nsim),
                            'DiRCT'=rep(NA,nsim),
                            'FAiRCT'=rep(NA,nsim),
                            'TSiRCT'=rep(NA,nsim),
                            'cRCT_gaussian_coxme'=rep(NA,nsim),
                            'cRCT_gaussian_coxph'=rep(NA,nsim),
                            'cRCT_gamma_coxph'=rep(NA,nsim),
                            'cRCT_gee'=rep(NA,nsim))
mle_VaccineEfficacy <- mle_pvals <- data.frame('iRCT'=rep(NA,nsim),
                                               'cRCT_gee'=rep(NA,nsim),
                                               'DiRCT'=rep(NA,nsim),
                                               'FAiRCT'=rep(NA,nsim),
                                               'TSiRCT'=rep(NA,nsim))
ICCs <- rep(NA,nsim)
deffs <- rep(NA,nsim)
props_zeros <- rep(NA,nsim)

list_trial_parameters <- function(# First day of trial enrollment, relative to start of epidemic
                                  trial_startday=100,
                                  # Number of days over which subjects are vaccinated
                                  vaccination_gap=1,
                                  # As defined for primary endpoint
                                  follow_up=40,
                                  revisit=0,
                                  trial_length=400,
                                  bCluster=0,
                                  bTrial=1,
                                  reevaluate=0,
                                  adaptation_day=0,
                                  adaptation='',
                                  # Number of days over which subjects are enrolled
                                  enrollment_period = 1,
                                  enrollment_gap = 1,
                                  # Target community enrollment proportion
                                  cluster_coverage = 0.75){
  list(trial_startday=trial_startday,
       vaccination_gap=vaccination_gap,
       follow_up=follow_up,
       revisit=revisit,
       trial_length=trial_length,
       bCluster=bCluster,
       bTrial=bTrial,
       reevaluate=reevaluate,
       adaptation_day=adaptation_day,
       adaptation=adaptation,
       enrollment_period = enrollment_period,
       enrollment_gap = enrollment_gap,
       # Number of clusters targeted for enrollment
       # Must be less than or equal to the number of communities
       num_enrolled_per_day = floor(num_communities/enrollment_period),
       cluster_coverage=cluster_coverage,
       name=paste0(ifelse(bCluster==1,'c','i'),ifelse(vaccination_gap==trial_length,'RCT',paste0(ifelse(bTrial==2,'Ring',''),ifelse(adaptation=='','FR',adaptation))),'-',follow_up))
}

trial_designs <- list()
trial_designs[[1]] <- list_trial_parameters(vaccination_gap=400,
                                            follow_up=400)
trial_designs[[2]] <- list_trial_parameters(vaccination_gap=400,
                                            follow_up=400,
                                            bCluster=1)
trial_designs[[3]] <- list_trial_parameters(reevaluate=1)
trial_designs[[4]] <- list_trial_parameters(adaptation_day = 40,
                                            adaptation='FA')
trial_designs[[5]] <- list_trial_parameters(adaptation_day = 40,
                                            adaptation='TS')
trial_designs[[6]] <- list_trial_parameters(adaptation_day = 40,
                                            adaptation='TST')
trial_designs[[7]] <- list_trial_parameters(bTrial=2,
                                            reevaluate=1)
trial_designs[[8]] <- list_trial_parameters(bTrial=2,
                                            adaptation_day = 40,
                                            adaptation='TS')
trial_designs[[9]] <- list_trial_parameters(bTrial=2,
                                            adaptation_day = 40,
                                            adaptation='TST')
trial_designs[[10]] <- list_trial_parameters(revisit=1,
                                             #follow_up=400,
                                            adaptation_day = 40,
                                            adaptation='TS')
trial_designs[[11]] <- list_trial_parameters(revisit=1,
                                             #follow_up=400,
                                            adaptation_day = 40,
                                            adaptation='TST')
trial_designs[[12]] <- list_trial_parameters(bTrial=2,
                                            adaptation_day = 40,
                                            revisit=1,
                                            follow_up=400,
                                            adaptation='TS')
trial_designs[[13]] <- list_trial_parameters(bTrial=2,
                                            adaptation_day = 40,
                                            revisit=1,
                                            #follow_up=400,
                                            adaptation='TST')
trial_designs[[14]] <- list_trial_parameters(revisit=1,
                                             #follow_up=400,
                                             adaptation_day = 40,
                                            adaptation='FA')

## get infection trajectory for source population
num_timesteps <- max(sapply(trial_designs,function(x) x$trial_startday + x$trial_length + x$enrollment_period - 1))
times <- seq(0,num_timesteps,1)
infected_trajectory <- get_infected_trajectory(times)

## set up for results
trials <- length(trial_designs)
extra_trials <- sum(sapply(trial_designs,function(x) x$reevaluate))
numevents_cont <- numevents <- numevents_vacc <- num_vacc <-  num_enrolled <- matrix(NA,nrow=trials+extra_trials,ncol=nsim)
trajectory_list <- list()
registerDoParallel(cores=10)
simnum <- sday <- 1
trial_outcomes <- list()
}

print(system.time(for(direct_VE in c(0,0.6)){ # sday in c(5:1)){ #
  for (simnum in 1:nsim) {
    if(direct_VE==0.6) trajectory_list[[simnum]] <- list() 
    g<-make_network(ave_community_size, community_size_range, num_communities,rate_within, rate_between)
    
    trial_outcomes[[sday]] <- foreach(tr = 1:trials) %dopar% {
      trial_startday <- trial_designs[[tr]]$trial_startday#100#50 + (sday-1)*100
      trial_length <- trial_designs[[tr]]$trial_length#500 - trial_startday
      bTrial <- trial_designs[[tr]]$bTrial
      bCluster <- trial_designs[[tr]]$bCluster
      #adaptation <- trial_designs[[tr]]$adaptation
      #vaccination_gap <- trial_designs[[tr]]$vaccination_gap
      follow_up <- trial_designs[[tr]]$follow_up
      #adaptation_day <- trial_designs[[tr]]$adaptation_day
      #profvis(
      list[results,trial_nodes,trajectories]<-
        network_epidemic(g,disease_dynamics,direct_VE,infected_trajectory,trial_designs[[tr]])
      #)
      
      list[VE,pval,events_vacc,events_cont,analysed_trialsize] <- 
        analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,bCluster,follow_up,trial_designs[[tr]]$revisit)
      VE <- VE[1]
      # if(tr!=2){
      ## add analysis for ring-end
        if(trial_designs[[tr]]$reevaluate==1){
          pvals$DiRCT[simnum] <- pval
          VaccineEfficacy$DiRCT[simnum] <- VE[1]
          # duplicate results with different end point
          list[VE2,pval2,events_vacc,events_cont,analysed_trialsize] <- 
            analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,bCluster,follow_up,revisit=1)
          pval <- c(pval2,pval)
          VE <- c(VE2[1],VE[1])
        }
      # }else{
        # list[VE_gaussian_coxme,pval_gaussian_coxme,VE_gaussian_coxph,pval_gaussian_coxph,
        #      VE_gamma_coxph,pval_gamma_coxph,VE_gee,pval_gee,events_vacc,events_cont,analysed_trialsize,
        #      ICC,deff,prop_zeros] <- analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,bCluster,vaccination_gap)
        # 
        # pvals$cRCT_gaussian_coxme[simnum] <- pval_gaussian_coxme
        # VaccineEfficacy$cRCT_gaussian_coxme[simnum] <- VE_gaussian_coxme[1]
        # 
        # pvals$cRCT_gaussian_coxph[simnum] <- pval_gaussian_coxph
        # VaccineEfficacy$cRCT_gaussian_coxph[simnum] <- VE_gaussian_coxph[1]
        # 
        # pvals$cRCT_gamma_coxph[simnum] <- pval_gamma_coxph
        # VaccineEfficacy$cRCT_gamma_coxph[simnum] <- VE_gamma_coxph[1]
        # 
        # pvals$cRCT_gee[simnum] <- pval_gee
        # VaccineEfficacy$cRCT_gee[simnum] <- VE_gee[1]
        # 
        # ICCs[simnum] <- ICC
        # deffs[simnum] <- deff
        # props_zeros[simnum] <- prop_zeros
        # 
        # pval <- pval_gee
        # VE <- VE_gee[1]
      #}
      num_enrolled <- nrow(trial_nodes)#analysed_trialsize
      print(num_enrolled)
      num_vacc <- sum(trial_nodes$TrialStatus==1)
      numevents <- nrow(results)
      list(num_enrolled=num_enrolled,num_vacc=num_vacc,events_vacc=events_vacc,events_cont=events_cont,numevents=numevents,pval=pval,VaccineEfficacy=VE,
           trajectories=trajectories,vaccinationDays=trial_nodes$DayVaccinated)
    }
    index <- 0
    for(tr in 1:trials){
      
      ## add analysis for ring-end
      index <- index + 1
      if(trial_designs[[tr]]$reevaluate==1) index <- c(index,index+1)
      num_enrolled[index,simnum] <- trial_outcomes[[sday]][[tr]]$num_enrolled
      num_vacc[index,simnum] <- trial_outcomes[[sday]][[tr]]$num_vacc
      numevents_vacc[index,simnum] <- trial_outcomes[[sday]][[tr]]$events_vacc
      numevents_cont[index,simnum] <- trial_outcomes[[sday]][[tr]]$events_cont
      numevents[index,simnum] <- trial_outcomes[[sday]][[tr]]$numevents
      mle_pvals[simnum,index] <- trial_outcomes[[sday]][[tr]]$pval
      mle_VaccineEfficacy[simnum,index] <- trial_outcomes[[sday]][[tr]]$VaccineEfficacy
      if(direct_VE==0.6) trajectory_list[[simnum]][[tr]] <- trial_outcomes[[sday]][[tr]]$trajectories
      index <- max(index)
    }
    cat("Simulation ",simnum,"\n")
    
  }
  assign(paste0('num_enrolled',direct_VE),num_enrolled)
  assign(paste0('mle_pvals',direct_VE),mle_pvals)
  #assign(paste0('mle_pvals',sday),mle_pvals)
}))

N <- 50000
y <- c(S=N-1,E1=0,E2=0,E3=0,I1=1,I2=0,I3=0,R=0)
parms <- c(betahat=0.94,a1=0.19,a2=0.6,atau=27.79,sigma=0.14,gamma=0.33)
out <- as.data.frame(lsoda(y,times,source_population_model,parms))
sources <- list()
sources[[1]] <- out$S
sources[[2]] <- out$E1 + out$E2 + out$E3
sources[[3]] <- out$I1 + out$I2 + out$I3
sources[[4]] <- out$R
pop <- c('S','E','I','R')
cols <- c('darkorange','navyblue','hotpink','turquoise')
pdf('source.pdf'); par(mfrow=c(2,2),mar=c(5,5,2,0.5))
for(k in 1:4){
  plot(1:length(sources[[1]]),sources[[k]],typ='l',col=cols[k],frame=F,lwd=2,xlab='Day',ylab=pop[k],cex.axis=1.5,cex.lab=1.5)
#  for(vl in seq(0,max(vaccination_gaps*(vaccination_periods-1)),by=max(vaccination_gaps))+trial_startday) abline(v=vl,lwd=2,col='grey',lty=2)
}
dev.off()

#rowlabels <- c('iRCT','cRCT','FR-end','FR-40','FA-40','TS-40','Ring-end','Ring-40','Ring-TS','TS-end','Ring-TS-end','FA-end')
rowlabels <- sapply(trial_designs,function(x) x$name)
for(i in length(trial_designs):1) 
  if(trial_designs[[i]]$reevaluate==1) {
    end_time <- trial_designs[[i]]$trial_length
    previous_end_time <- trial_designs[[i]]$follow_up
    new_name <- gsub(previous_end_time,end_time,trial_designs[[i]]$name)
    rowlabels <- c(rowlabels[1:(i-1)],new_name,rowlabels[i:length(rowlabels)])
  }

# plot_inf <- sources[[3]][50:450]
# power_plot <- matrix(0,nrow=5,ncol=7)
# for(i in 1:5) power_plot[i,] <- apply(get(paste0('mle_pvals',i)),2,function(x)sum(x<0.05,na.rm=T)/sum(!is.na(x)))
# cols <- rainbow(7)
# pdf('powerplot.pdf'); par(mar=c(5,5,2,2),mfrow=c(1,1))
# plot(50+0:4*100,power_plot[,1],typ='b',lwd=2,col=cols[1],xlab='Start day',ylab='Power',cex.lab=1.5,cex.axis=1.5,frame=F,ylim=range(power_plot))
# for(i in 2:7) lines(50+0:4*100,power_plot[,i],typ='b',lwd=2,col=cols[i])
# lines(50:450,plot_inf/max(plot_inf),lwd=2,col=col.alpha('grey',0.6))
# legend(x=50,y=0.7,legend=rowlabels,col=cols,lwd=3,bty='n')
# dev.off()

sd_tab <- results_tab <- matrix(0,nrow=length(rowlabels),ncol=9)
results_tab[,1] <- rowMeans(numevents )
sd_tab[,1] <- apply(numevents,1,sd)
results_tab[,2] <- rowMeans(num_vacc)
sd_tab[,2] <- apply(num_vacc,1,sd)
results_tab[,3] <- rowMeans(numevents_vacc)
sd_tab[,3] <- apply(numevents_vacc,1,sd)
results_tab[,4] <- rowMeans(numevents/num_vacc)
sd_tab[,4] <- apply(numevents/num_vacc,1,sd)
results_tab[,5] <- rowMeans(num_enrolled)
sd_tab[,5] <- apply(num_enrolled,1,sd)
results_tab[,6] <- rowMeans(num_enrolled0)
sd_tab[,6] <- apply(num_enrolled0,1,sd)
results_tab[,7] <- apply(mle_pvals,2,function(x)sum(x<0.05,na.rm=T)/sum(!is.na(x))) # power
results_tab[,8] <- apply(mle_pvals0,2,function(x)sum(x<0.05,na.rm=T)/sum(!is.na(x))) # type 1 error rate
results_tab[,9] <- apply(mle_VaccineEfficacy,2,function(x)mean(x[is.finite(x)],na.rm=T))
sd_tab[,9] <- apply(mle_VaccineEfficacy,2,function(x)sd(x[is.finite(x)],na.rm=T))
signifs <- c(3,4,3,2,4,4,2,2,2)
for(i in 1:length(rowlabels)){
  cat(rowlabels[i])
  for(j in 1:ncol(results_tab)){
    cat(paste0(' & ',signif(results_tab[i,j],signifs[j])))
    if(!j%in%c(7,8)) cat(paste0(' (',signif(sd_tab[i,j],2),')'))
  }
  cat('\\\\\n')
  if(i%%3==0)cat('\\hline\n')
}

pop <- c('S','New E','New I','New R')
cols <- c('darkorange','navyblue','hotpink','turquoise')
methods <- c('iRCT','cRCT','FR-iRCT','FA-iRCT','TS-iRCT','Ring')
maxmins <- lapply(1:4,function(k)range(sapply(trajectory_list,function(x)sapply(x,function(y)y[[k]]))))
for(i in 1:length(methods)){
  pdf(paste0(methods[i],'.pdf')); par(mfrow=c(2,2),mar=c(5,5,2,0.5))
  for(k in 1:4)
    for(j in 1:5){
      mat <- do.call(rbind,trajectory_list[[j]][[i]])
      total <- colSums(mat)
      props <- mat#/total
      if(j==1) plot(1:ncol(props),props[k,],typ='l',col=col.alpha(cols[k],ifelse(k==1,1,0.25)),frame=F,lwd=2,xlab='Day',ylab=pop[k],cex.axis=1.5,cex.lab=1.5,ylim=maxmins[[k]])
      if(j>1) lines(1:ncol(props),props[k,],lwd=2,col=col.alpha(cols[k],ifelse(k==1,1,0.25)))
      if(j==1){
        abline(v=trial_startday,lwd=2,col='grey',lty=2)
        if(i>2) for(vl in seq(0,max(vaccination_gaps*(vaccination_periods-1)),by=max(vaccination_gaps))+trial_startday) abline(v=vl,lwd=2,col='grey',lty=2)
      }
    }
  dev.off()
}


### plot vaccination and follow-up days
example_routines <- c('iRCT','FR-400','FR-40','Ring')
trial_numbers <- c(1,3,3,7)
trial_end <- 500#trial_startday+trial_length
pdf('enrollment.pdf',width=10); par(mfrow=c(2,2),mar=c(5,5,2,2))
for(i in 1:4){
  vdays <- hist(trial_outcomes[[1]][[trial_numbers[i]]]$vaccinationDays,breaks=seq(0,trial_end,by=20),plot=F)
  fdays <- if(i%in%c(1,2)){
    rep(trial_end,length=length(trial_outcomes[[1]][[trial_numbers[i]]]$vaccinationDays))
  }else if(i%in%c(3,4)){
    trial_outcomes[[1]][[trial_numbers[i]]]$vaccinationDays + 40
  }
  fdays <- fdays[fdays<=trial_end]
  fdaysfreq <- hist(fdays,breaks=seq(0,trial_end,by=20),plot=F)
  toplot <- rbind(fdaysfreq$counts,vdays$counts)
  barplot(toplot,col=c('navyblue','darkorange'),beside=T,xlab='Day',ylab='Number of people',cex.axis=1.5,cex.lab=1.5,main=example_routines[i])
  axis(1,at=seq(0,500,by=100)/7,labels=seq(0,500,by=100),cex.axis=1.5,cex.lab=1.5)
  if(i==1) legend(x=200/7,y=5500,legend = c('Enrolled','Followed up'),fill=c('darkorange','navyblue'),bty='n')
}
dev.off()



# Final size calculations

# Compare iRCT to cRCT
# Vary VE, R0, number of introductions, and/or enrollment percentage
# R0s <- c(seq(1.2,1.7,0.005))
# community_size <- 100
# num_communities<-100
# num_intros_percluster_vec <- 1
# VEs <- seq(0.3,0.6,0.005)
# enrollment_percs<-c(0.6)
# 
# num_params <- length(R0s)*length(num_intros_percluster_vec)*length(VEs)*length(enrollment_percs)
# # Assumption: size of minor outbreaks, when R0>1 but no major outbreak occurs
# AR_nonoutbreak <- 1/community_size
# 
# index <- 1
# 
# results <- data.frame("R0"=rep(NA,num_params),
#                       "num_intros"=rep(NA,num_params),
#                       "VE"=rep(NA,num_params),
#                       "enrollment_perc"=rep(NA,num_params),
#                       "estim_VE"=rep(NA,num_params),
#                       "estim_VE_iRCT"=rep(NA,num_params),
#                       "cRCT_trial_CI"=rep(NA,num_params),
#                       "cRCT_trial_CI_vacc"=rep(NA,num_params),
#                       "cRCT_trial_CI_unvacc"=rep(NA,num_params),
#                       "iRCT_trial_CI"=rep(NA,num_params),
#                       "iRCT_trial_CI_vacc"=rep(NA,num_params),
#                       "iRCT_trial_CI_unvacc"=rep(NA,num_params),
#                       "ICC"=rep(NA,num_params),
#                       "samplesize_cRCT"=rep(NA,num_params),
#                       "samplesize_iRCT"=rep(NA,num_params))
# 
# for (R0 in R0s) { 
#   for (num_intros_percluster in num_intros_percluster_vec) {
#     for (direct_VE in VEs) {
#       for (enrollment_perc in enrollment_percs) {
#         
#         results$R0[index]<-R0
#         results$num_intros[index]<-num_intros_percluster
#         results$VE[index]<-direct_VE
#         results$enrollment_perc[index]<-enrollment_perc
#         
#         R0_vacc <- (1 - enrollment_perc * direct_VE) * R0
#         
#         # Final size function to solve
#         
#         # Calculate final sizes
#         if (R0_vacc>1) {
#           CI_vacc<-multiroot(finalsizes_vacc,c(0.5,0.5),parms=list(R0=R0,enrollment_perc=enrollment_perc,direct_VE=direct_VE))$root[2]
#         } else {
#           # This is true in large populations but is not a good approximation
#           # in small populations and when R0 is close to 1
#           CI_vacc <- 1/(1-R0_vacc)/community_size
#         }
#         if (R0>1) {
#           CI_unvacc<-uniroot(finalsize_unvacc,c(1e-05,1),parms=list(R0=R0))$root
#         } else {
#           # This is true in large populations but is not a good approximation
#           # in small populations and when R0 is close to 1
#           CI_unvacc <- 1/(1-R0)/community_size
#         }
#         
#         # CI of vaccinated individuals in vaccinated clusters, accounting for outbreak probs
#         if (R0_vacc > 1) {
#           p_vacc <- 1-uniroot(outbreakprob_vacc,c(0,1-1e-05),parms=list(R0_vacc=R0_vacc))$root
#           cRCT_CI_vacc <- num_intros_percluster*(CI_vacc*p_vacc + (1-p_vacc)*AR_nonoutbreak)
#           
#           p_unvacc <- CI_unvacc
#           cRCT_CI_unvacc <- num_intros_percluster*(CI_unvacc*p_unvacc+(1-p_unvacc)*AR_nonoutbreak)
#           
#           # CI_O, overall CI in the trial population
#           results$cRCT_trial_CI[index] <- (cRCT_CI_unvacc+cRCT_CI_vacc)/2
#           results$cRCT_trial_CI_vacc[index] <- cRCT_CI_vacc
#           results$cRCT_trial_CI_unvacc[index] <- cRCT_CI_unvacc
#           # Estimated VE for hazard rate analysis
#           results$estim_VE[index]<-1-log(1-cRCT_CI_vacc)/log(1-cRCT_CI_unvacc)
#           
#           # Calculate ICC by calculating within- and between- cluster sum of squares
#           # and applying the ANOVA method
#           ARs <- c(rep(CI_vacc,round(num_communities*num_intros_percluster*p_vacc/2)),
#                    rep(CI_unvacc,round(num_communities*num_intros_percluster*p_unvacc/2)),
#                    rep(AR_nonoutbreak,round(num_communities*
#                                               (1-num_intros_percluster*p_vacc/2-num_intros_percluster*p_unvacc/2))))
#           cluster_size <- enrollment_perc*community_size
#           numevents <- ARs*cluster_size
#           
#           cluster_sizes <- rep(cluster_size,num_communities)
#           N <- sum(cluster_sizes)
#           K <- num_communities
#           n0 <- 1/(K-1) * (N - K*cluster_size^2/N)
#           n01 <- 1/(K-1) * ((K-1)*n0 - K*cluster_size^2/N)
#           MSB <- 1/(K-1) * sum((numevents-mean(numevents))^2/cluster_size)
#           MSW <- 1/(N-K-1) * sum(numevents-numevents^2/cluster_size)
#           ICC <- (MSB - MSW) / (MSB + (n01-1) * MSW)
#           
#           results$ICC[index]<-ICC
#           DE<-1+(enrollment_perc*community_size-1)*ICC
#           results$samplesize_cRCT[index] <- 
#             (1.96+1.28)^2/log(1-results$estim_VE[index])^2 * DE / results$cRCT_trial_CI[index]
#           
#           
#         } else {
#           cRCT_CI_vacc <- num_intros_percluster*CI_vacc
#         }
#         
#         R0_iRCT <- (1 - 0.5 * enrollment_perc * direct_VE) * R0
#         if (R0_iRCT>1) {
#           CIs_iRCT<-multiroot(finalsizes_iRCT,c(0.5,0.5),parms=list(R0=R0,enrollment_perc=enrollment_perc,direct_VE=direct_VE))
#           CI_vacc_iRCT<-CIs_iRCT$root[2]
#           CI_unvacc_iRCT<-CIs_iRCT$root[1]
#           
#           p_vacc_iRCT <- 1-uniroot(outbreakprob_iRCT,c(0,1-1e-05),parms=list(R0_iRCT=R0_iRCT))$root
#           iRCT_CI_vacc <- num_intros_percluster*(CI_vacc_iRCT*p_vacc_iRCT + (1-p_vacc_iRCT)*AR_nonoutbreak)
#           iRCT_CI_unvacc <- num_intros_percluster*(CI_unvacc_iRCT*p_vacc_iRCT + (1-p_vacc_iRCT)*AR_nonoutbreak)
#           
#           results$estim_VE_iRCT[index] <- 1-log(1-iRCT_CI_vacc)/log(1-iRCT_CI_unvacc)
#           results$iRCT_trial_CI[index] <- (iRCT_CI_vacc+iRCT_CI_unvacc)/2
#           results$iRCT_trial_CI_vacc[index] <- iRCT_CI_vacc
#           results$iRCT_trial_CI_unvacc[index] <- iRCT_CI_unvacc
#           results$samplesize_iRCT[index] <- (1.96+1.28)^2/log(1-results$estim_VE_iRCT[index])^2 / results$iRCT_trial_CI[index]
#         }
#         index<-index+1
#       }
#     }
#   }
# }
# 
# results$ssratio <- results$samplesize_cRCT/results$samplesize_iRCT
# 
# # Example plot of sample size ratio, varying VE and R0
# x<-matrix(results$ssratio,nrow = length(VEs),ncol=length(R0s),byrow = FALSE)
# x11();image.plot(VEs,R0s,x,
#            legend.lab="cRCT/iRCT sample size ratio",
#            xlab="Vaccine efficacy",
#            ylab=expression(R[0]))
# 
# 
