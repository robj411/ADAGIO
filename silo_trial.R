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

## build network ###############################################################

source('build_network.R')

## functions #######################################################

source('network_functions.R')
source('evaluation_functions.R')

probability_by_lag <<- readRDS(paste0('probability_by_lag_8080.Rds'))
probability_after_day_0 <<- readRDS(paste0('probability_after_day_0_8080.Rds'))
probability_by_lag_given_removal <<- readRDS(paste0('probability_by_lag_given_removal_8080.Rds'))
probability_after_day_0_given_removal <<- readRDS(paste0('probability_after_day_0_given_removal_8080.Rds'))

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

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
ves <- c(0,0.9)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- c(0,1)
trial_designs <- expand.grid(VE=ves,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
trial_designs <- rbind(trial_designs,trial_designs[trial_designs$adapt=='',])
trial_designs$weight[(nCombAdapt+1):(nComb*(length(adaptations)+1))] <- 'binary'
trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
registerDoParallel(cores=16)
func <- get_efficacious_probabilities
eval_day <- 31
latest_infector_time <- eval_day - 0

for(rnd in 2:1){
  trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
    if(rnd==1){
      func <- get_efficacious_probabilities
    }else{
      func <- get_efficacious_probabilities2
    }
    cluster_flag <- trial_designs$cluster[des]
    direct_VE <- trial_designs$VE[des]
    adaptation <- trial_designs$adapt[des]
    vaccinated_count <- infectious_count <- vaccinated_countb <- enrolled_count <- enrolled_countb <- infectious_countb <- 0
    pval_binary_mle3 <- pval_binary_mle2 <- pval_binary_mle <- ve_est3 <- ve_est2 <- ve_est <- c()
    for(tr in 1:nTrials){
      vaccinees <- trial_participants <- recruit_times <- c()
      vaccinees2 <- trial_participants2 <- c()
      infectious_by_vaccine <- excluded <- matrix(0,nrow=nClusters,ncol=2)
      results_list <- list()
      allocation_ratio <- 0.5
      netwk_list <- list()
      for(iter in 1:nClusters){
        ## select random person to start
        first_infected <- sample(g_name,1)
        inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
        #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
        hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
        inf_time <- min(inf_period,hosp_time)
        netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,inf_time,cluster_flag=cluster_flag,allocation_ratio=allocation_ratio,direct_VE=direct_VE)
        netwk_list[[iter]] <- netwk
        results_list[[iter]] <- netwk[[1]]
        results <- results_list[[iter]]
        infectious_by_vaccine[iter,] <- c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+9),sum(!results$vaccinated&results$inTrial&results$DayInfectious>results$RecruitmentDay+9))
        excluded[iter,] <- c(sum(results$vaccinated&results$DayInfectious<results$RecruitmentDay+10),sum(!results$vaccinated&results$inTrial&results$DayInfectious<results$RecruitmentDay+10))
        recruit_times[iter] <- netwk[[3]]
        vaccinees[iter] <- netwk[[4]]
        trial_participants[iter] <- netwk[[5]]
        
        ##!! weighting non-events
        rec_day <- recruit_times[iter]
        infectious_index <- results$DayInfectious<latest_infector_time+rec_day&(results$DayRemoved>rec_day|is.na(results$DayRemoved))
        infectious_names <- results$InfectedNode[infectious_index]
        infectious_ends <- pmin(results$DayRemoved[infectious_index],latest_infector_time+rec_day)
        infectious_ends[is.na(infectious_ends)] <- latest_infector_time+rec_day
        infectious_starts <- pmax(results$DayInfectious[infectious_index],rec_day)
        if(identical(func,get_efficacious_probabilities2)){
          vaccinees[iter] <- trial_participants[iter] <- 0
          if(length(infectious_names)>0){
            popweights <- rowSums(sapply(1:length(infectious_names),function(i){
              x <- infectious_names[i]
              # prepare contacts
              contacts <- contact_list[[x]]
              c_of_c <- contact_of_contact_list[[x]]
              hr <- c(high_risk_list[[x]],household_list[[x]])
              # prepare trial participants
              vax <- netwk[[6]]
              cont <- netwk[[7]]
              # work out total risk presented by infector
              infector_weight <- sum(pgamma(eval_day-infectious_starts[i]:infectious_ends[i],shape=incperiod_shape,rate=incperiod_rate))
              # remove anyone infectious earlier
              earlier_nodes <- results$InfectedNode[results$DayInfectious<infectious_starts[i]]
              contacts <- contacts[!contacts%in%earlier_nodes]
              c_of_c <- c_of_c[!c_of_c%in%earlier_nodes]
              hr <- hr[!hr%in%earlier_nodes]
              # sum of person days times scalars
              total_vax <- sum(vax%in%contacts) + neighbour_scalar*sum(vax%in%c_of_c) + (high_risk_scalar-1)*sum(vax%in%hr)
              total_cont <- sum(cont%in%contacts) + neighbour_scalar*sum(cont%in%c_of_c) + (high_risk_scalar-1)*sum(cont%in%hr)
              c(total_vax,total_cont)*infector_weight
            }))
            if(length(netwk[[6]])>0)
              vaccinees[iter] <- popweights[1]
            trial_participants[iter] <- popweights[2]
          }
        }
        vaccinees2[iter] <- netwk[[4]]
        trial_participants2[iter] <- netwk[[5]]
        
        ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
        if(adaptation!=''&&iter %% eval_day == 0 && sum(vaccinees)>0){
          allocation_ratio <- response_adapt(results_list,vaccinees,trial_participants,adaptation,func=func)
        }
      }
      
      ves <- c(0.6,1)
      while(abs(ves[1]-ves[2])>0.1){
        trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=ves[1])
        
        tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
        form <- 'Surv(time,outcome)~vaccinated'
        if(cluster_flag==1) form <- paste0(form,'+strata(cluster)')
        survmodel<-coxph(as.formula(form),weights=weight,tte)
        vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
        zval <- survmodel$coefficient/sqrt(survmodel$var)
        pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
        #print(c(vaccEffEst,zval,pval))
        ves[2] <- ves[1]
        if(!is.na(vaccEffEst[1]))
          ves[1] <- vaccEffEst[1]
      }
      pval_binary_mle3[tr]  <- pval#calculate_pval(eval_list[[3]],eval_list[[2]])
      ve_est3[tr]  <- ves[1]#eval_list[[1]]
      
      
      if(rnd==1){
        eval_list <- func(results_list,vaccinees2,trial_participants2)
      }else{
        eval_list <- func(results_list,vaccinees,trial_participants)
      }
      pval_binary_mle2[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
      ve_est2[tr]  <- eval_list[[1]]
      vaccinated_count <- vaccinated_count + sum(vaccinees2)/nTrials
      enrolled_count <- enrolled_count + sum(trial_participants2)/nTrials
      infectious_count <- infectious_count + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
      if(adaptation==''){
        pop_sizes <- c(sum(vaccinees2),sum(trial_participants2) - sum(vaccinees2)) - colSums(excluded)
        pval_binary_mle[tr] <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
        ve_est[tr]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
        vaccinated_countb <- vaccinated_countb + sum(vaccinees2)/nTrials
        enrolled_countb <- enrolled_countb + sum(trial_participants2)/nTrials
        infectious_countb <- infectious_countb + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
      }
      ## ICC without weighting
      #if(cluster_flag==1){
      #  vax <- vaccinees
      #  non_vax <- trial_participants - vax
      #  trial_case <- sapply(results_list,function(x)sum(x$inTrial==T))
      #  vax_case <- sapply(results_list,function(x)sum(x$vaccinated==T))
      #  non_vax_case <- trial_case - vax_case
      #  cid <- rep(1:length(trial_participants),times=trial_participants)
      #  non_cases <- trial_participants - trial_case
      #  y <- unlist(sapply(1:length(trial_case),function(x) c(rep(1,times=trial_case[x]),rep(0,times=non_cases[x]))))
      #icc <- iccbin(cid,y,data=data.frame(cid=factor(cid),y=y),method='aov',ci.type='aov')
      #}
    }
    print(c(des,adaptation))
    power <- sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
    print(c(des,power,sum(pval_binary_mle3<0.05,na.rm=T)/sum(!is.na(pval_binary_mle3))))
    VE_est <- mean(ve_est2,na.rm=T)
    print(c(des,VE_est,mean(ve_est3,na.rm=T)))
    VE_sd <- sd(ve_est2,na.rm=T)
    powerb <- VE_estb <- VE_sdb <- c()
    if(adaptation==''){
      powerb <- sum(pval_binary_mle<0.05,na.rm=T)/sum(!is.na(pval_binary_mle))
      VE_estb <- mean(ve_est,na.rm=T)
      VE_sdb <- sd(ve_est,na.rm=T)
    }
    return(list(power, VE_est, VE_sd,powerb, VE_estb, VE_sdb,vaccinated_count, infectious_count, vaccinated_countb, infectious_countb,enrolled_count,enrolled_countb))
  }
  for(des in 1:nCombAdapt){
    cluster_flag <- trial_designs$cluster[des]
    direct_VE <- trial_designs$VE[des]
    adaptation <- trial_designs$adapt[des]
    trial_designs$vaccinated[des] <- trial_results[[des]][[7]]
    trial_designs$infectious[des] <- trial_results[[des]][[8]]
    trial_designs$enrolled[des] <- trial_results[[des]][[11]]
    if(adaptation==''){
      trial_designs$vaccinated[des+nComb] <- trial_results[[des]][[9]]
      trial_designs$infectious[des+nComb] <- trial_results[[des]][[10]]
      trial_designs$enrolled[des+nComb] <- trial_results[[des]][[12]]
    }
    trial_designs$power[des] <- trial_results[[des]][[1]]
    trial_designs$VE_est[des] <- trial_results[[des]][[2]]
    trial_designs$VE_sd[des] <- trial_results[[des]][[3]]
    if(adaptation==''){
      trial_designs$power[des+nComb] <- trial_results[[des]][[4]]
      trial_designs$VE_est[des+nComb] <- trial_results[[des]][[5]]
      trial_designs$VE_sd[des+nComb] <- trial_results[[des]][[6]]
    }
  }
  subset(trial_designs,VE==0)
  subset(trial_designs,VE>0)
  
  result_table <- subset(trial_designs,VE>0)[,c(2,3,4,5,6,7,8,9,10)]
  result_table$t1e <- subset(trial_designs,VE==0)$power
  result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
  result_table <- result_table[,!colnames(result_table)%in%c('VE_est','VE_sd')]
  result_table$adapt <- as.character(result_table$adapt)
  result_table$adapt[result_table$adapt==''] <- 'None'
  result_table$cluster[result_table$cluster==0] <- 'Individual'
  result_table$cluster[result_table$cluster==1] <- 'Cluster'
  colnames(result_table) <- c('Randomisation','Adaptation','Weighting','Sample size','Infectious','Vaccinated','Power','Type 1 error','VE estimate')
  print(xtable(result_table), include.rownames = FALSE)

}
