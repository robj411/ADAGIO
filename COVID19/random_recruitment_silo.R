source('set_up_script.R')
registerDoParallel(cores=12)

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
vaccine_efficacies <- c(0,0.7)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
ref_recruit_day <- 30
latest_infector_time <- eval_day - 0

## store original lists
true_contact_of_contact_list <- contact_of_contact_list
true_household_list <- household_list
true_contact_list <- contact_list
## get list sizes
total_size <- mean(sapply(1:length(contact_list),function(x){
  c_of_c <- contact_of_contact_list[[x]]
  cont <- contact_list[[x]]
  length(unique(c(cont,c_of_c)))
}))
## overwrite with empty lists
household_list <- contact_of_contact_list <- lapply(g_name,function(x)c())
## fill contact list with random indices of mean total size
contact_list <- lapply(g_name,function(x)sample(g_name,total_size,replace=F))

covid_spread_wrapper <- function(i_nodes_info,s_nodes,v_nodes,e_nodes_info,direct_VE){
  contact_list <- true_contact_list
  household_list <- true_household_list
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

eligible_first_person <- sapply(contact_list,length)>10


trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
  set.seed(des)
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  vaccinated_count <- infectious_count <- rr_list <- list()
  for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- 0
  pval_binary_mle3 <- ve_est3 <- pval_binary_mle2 <- ve_est2 <- pval_binary_mle <- ve_est <- ve_estht <- c()
  exports <- enrolled_count <- c()
  for(tr in 1:nTrials){
    randomisation_ratios <- c()
    people_per_ratio <- c()
    vaccinees <- trial_participants <- c()
    infectious_by_vaccine <- excluded <- c()
    results_list <- list()
    allocation_ratio <- 0.5
    netwk_list <- list()
    weight_break <- 0
    iter <- 0
    while(weight_break<target_weight){
      iter <- iter + 1
    #for(iter in 1:nClusters){
      ## select random person to start
      randomisation_ratios[iter] <- allocation_ratio
      first_infected <- sample(g_name[eligible_first_person],1)
      netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,end_time=eval_day,allocation_ratio=allocation_ratio,direct_VE=direct_VE,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper,from_source=1e-5)
      netwk_list[[iter]] <- netwk
      results_list[[iter]] <- netwk[[1]]
      results <- results_list[[iter]]
      infectious_by_vaccine <- rbind(infectious_by_vaccine,c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+9),sum(!results$vaccinated&results$inTrial&results$DayInfectious>results$RecruitmentDay+9)))
      excluded <- rbind(excluded,c(sum(results$vaccinated&results$DayInfectious<results$RecruitmentDay+10),sum(!results$vaccinated&results$inTrial&results$DayInfectious<results$RecruitmentDay+10)))
      
      vaccinees[iter] <- netwk[[4]]
      trial_participants[iter] <- netwk[[5]]
      
      ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
      if(adaptation!=''&&iter %% eval_day == 0 && sum(vaccinees)>0){
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
        pop_sizes2 <- probs[[2]]
        fails <- probs[[3]]
        allocation_ratio <- response_adapt(fails,pop_sizes2,days=iter,adaptation)
        people_per_ratio <- rbind(people_per_ratio,c(sum(trial_participants),iter,allocation_ratio))
        #if(allocation_ratio==0) break
        weight_break <- sum(probs[[3]])
      }else if(iter >= eval_day && sum(vaccinees)>0){
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
        weight_break <- sum(probs[[3]])
      }
    }
    rr_list[[tr]] <- people_per_ratio
    ## regular test
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1,observed=observed)
    pval_binary_mle2[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]
    ## correct VE test
    #eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,randomisation_ratios=randomisation_ratios,#rbht_norm=0,
    #                                           rbht_norm=ifelse(adaptation=='',1,2),
    #                                           people_per_ratio=people_per_ratio,adaptation=adaptation,contact_network=-1,observed=observed)#adaptation=adapt if rbht_norm=2
    ve_estht[tr]  <- eval_list[[1]]
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees)/nTrials
    enrolled_count[tr] <- sum(trial_participants)
    infectious_count[[1]] <- infectious_count[[1]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
    if(adaptation==''){
      pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
      pval_binary_mle[tr] <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      ve_est[tr]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      vaccinated_count[[2]] <- vaccinated_count[[2]] + sum(vaccinees)/nTrials
      infectious_count[[2]] <- infectious_count[[2]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
    }
    ## if a test was done
    #eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=T,contact_network=-1,observed=observed)
    #pval_binary_mle3[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    #ve_est3[tr]  <- eval_list[[1]]
    ## correcting for trend 
    pval_binary_mle3[tr]  <- NA
    ve_est3[tr]  <- NA
    if(adaptation!='')
      pval_binary_mle3[tr] <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                                 tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio,observed=observed)
    
    ## exports
    exports[tr] <- sum(sapply(results_list,function(x)sum(!x$inCluster)-1))/length(results_list)*100
  }
  power <- VE_est <- VE_sd <- c()
  power[1] <- sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
  VE_est[1] <- mean(ve_est2,na.rm=T)
  VE_sd[1] <- sd(ve_est2,na.rm=T)
  power[3] <- sum(pval_binary_mle2<pval_binary_mle3,na.rm=T)/sum(!is.na(pval_binary_mle3)&!is.na(pval_binary_mle2))
  VE_est[3] <- mean(ve_est3,na.rm=T)
  VE_sd[3] <- sd(ve_est3,na.rm=T)
  VE_est[4] <- mean(ve_estht,na.rm=T)
  VE_sd[4] <- sd(ve_estht,na.rm=T)
  #if(adaptation==''){
  #  power[2] <- sum(pval_binary_mle<0.05,na.rm=T)/sum(!is.na(pval_binary_mle))
  #  VE_est[2] <- mean(ve_est,na.rm=T)
  #  VE_sd[2] <- sd(ve_est,na.rm=T)
  #}
  #print(list(des, power, VE_est, VE_sd,vaccinated_count, infectious_count, enrolled_count,mean(exports)))
  power[2] <- quantile(pval_binary_mle2,0.95) - quantile(pval_binary_mle2,0.05)
  enrolled <- list(mean(enrolled_count),sd(enrolled_count))
  print(des)
  print(enrolled)
  return(list(power, VE_est, VE_sd,vaccinated_count, infectious_count, enrolled,rr_list,mean(exports)))
}
trial_designs$prange <- trial_designs$mee <- trial_designs$powertst <- trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- 
  trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolledsd <- trial_designs$enrolled <- 0
for(des in 1:nCombAdapt){
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  trial_designs$vaccinated[des] <- round(trial_results[[des]][[4]][[1]])
  trial_designs$infectious[des] <- round(trial_results[[des]][[5]][[1]])
  trial_designs$enrolled[des] <- round(trial_results[[des]][[6]][[1]])
  trial_designs$enrolledsd[des] <- round(trial_results[[des]][[6]][[2]])
  trial_designs$power[des] <- trial_results[[des]][[1]][1]
  trial_designs$prange[des] <- trial_results[[des]][[1]][2]
  trial_designs$VE_est[des] <- trial_results[[des]][[2]][1]
  trial_designs$VE_sd[des] <- trial_results[[des]][[3]][1]
  trial_designs$powertst[des] <- trial_results[[des]][[1]][3]
  trial_designs$mee[des] <- trial_results[[des]][[8]]
}
subset(trial_designs,VE==0)
subset(trial_designs,VE>0)
result_table <- subset(trial_designs,VE>0)[,-c(1:2)]
result_table$t1e <- subset(trial_designs,VE==0)$power
result_table$t1etst <- subset(trial_designs,VE==0)$powertst
result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
#result_table$power <- paste0(round(result_table$power,2),' (',round(result_table$prange,2),')')
result_table$enrolled <- paste0(result_table$enrolled,' (',result_table$enrolledsd,')')
result_table$adapt <- as.character(result_table$adapt)
result_table$adapt[result_table$adapt==''] <- 'None'
result_table$nmee <- subset(trial_designs,VE==0)$mee - subset(trial_designs,VE>0)$mee
result_table <- result_table[,!colnames(result_table)%in%c('VE_est','VE_sd','enrolledsd','mee','prange')]
colnames(result_table) <- c('Adaptation','Weighting','Sample size','Infectious','Vaccinated','Power','Power (corrected)',
                            'Type 1 error','Type 1 error (corrected)','VE estimate','NMEE')
print(xtable(result_table,digits=c(0,0,0,0,0,0,2,2,2,2,0,2)), include.rownames = FALSE)

saveRDS(trial_results,'storage/silo_trial_results.Rds')

change_days <- trial_results[[1]][[7]][[1]][,2]
adaptation_days <- c(1:change_days[1],c(sapply(2:length(change_days),function(x)(1+change_days[x-1]):change_days[x])))
get_allocation_vector <- function(x){
  sapply(1:5,function(y){
    vals <- trial_results[[x]][[7]][[y]][,3]
    alloc <- c(rep(0.5,change_days[1]),c(sapply(2:length(change_days),function(x)rep(vals[x-1],(change_days[x]-change_days[x-1])))))
    alloc
  })
}

# first index is the trial; 5:8 is TS designs
# second index is 7, extracting the rr_list
# third index is y over the number of samples
# fourth index ,3 extracts the ratio column
## propensities to stop for each 
print(sapply(5:8,function(x)
  apply(
    t(sapply(1:length(trial_results[[x]][[7]]),function(y){
      trial_results[[x]][[7]][[y]][,3]
    })),2,function(z) sum(z>0.99)/length(trial_results[[x]][[7]]))
))

## overall propensity to stop
print(sapply(5:8,function(x)
  sum(apply(
    t(sapply(1:length(trial_results[[x]][[7]]),function(y){
      trial_results[[x]][[7]][[y]][,3]
    })),1,function(z) any(z>0.99)))/
    length(trial_results[[x]][[7]])
))


## sample size for each 
print(sapply(5:8,function(x){
  sss <- sapply(1:length(trial_results[[x]][[7]]),function(y){
    rr <- trial_results[[x]][[7]][[y]]
    rs <- rr[,3]
    ss <- rr[,1]
    if(rs[1]>0.99) ss[1]
    else if(rs[2]>0.99) ss[2]
    else if(rs[3]>0.99) ss[3]
    else ss[3]/93*100
  })
  c(mean(sss),sd(sss))
}
))

