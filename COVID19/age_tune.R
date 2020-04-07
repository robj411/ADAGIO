source('set_up_script.R')
simulate_contact_network_original <- simulate_contact_network
simulate_contact_network_age <- function(first_infected,individual_recruitment_times=F,end_time=31,start_day=0,from_source=0,cluster_flag=0,allocation_ratio=0.5,
                                     direct_VE=0,base_rate=0,spread_wrapper=ebola_spread_wrapper){
  # set up info to store
  trajectories <- list()
  trajectories$S <- length(vertices) - 1
  trajectories$E <- 0
  trajectories$I <- 0
  trajectories$R <- 0
  e_nodes_info <- matrix(nrow=0,ncol=3)
  i_nodes_info <- matrix(nrow=0,ncol=4)
  e_nodes <- rep(0,length(g_name))
  i_nodes <- rep(0,length(g_name))
  v_nodes <- rep(0,length(g_name))
  c_nodes <- rep(0,length(g_name))
  s_nodes <- rep(1,length(g_name))
  r_nodes <- rep(0,length(g_name))
  t_nodes <- rep(0,length(g_name))
  
  # generate info for index case
  inc_time <- rgamma(length(first_infected),shape=incperiod_shape,rate=incperiod_rate)
  ceil_inc_time <- ceiling(inc_time)
  #i_nodes_info <- rbind(i_nodes_info,c(first_infected,rep(0,length(first_infected)),inf_time,inc_time))
  e_nodes_info <- rbind(e_nodes_info,c(first_infected,0,inc_time))
  s_nodes[first_infected] <- 0
  e_nodes[first_infected] <- 1
  
  #recruitment_time <- round(rgamma(1,shape=recruit_shape,rate=recruit_rate))
  recruitment_time <- ceiling(rtruncnorm(1,a=0,mean=recruit_mean,sd=recruit_sd))
  results <- matrix(nrow=0,ncol=4)#c(first_infected,0,-inc_time,NA),nrow=1)
  numinfectious <- 1
  ##!! add in additional infectious people?
  
  # identify contacts of index case
  contacts <- contact_list[[first_infected]]
  contacts <- contacts[contacts!=first_infected]
  order_infected <- first_infected
  ## identify high-risk people
  ##!! all household members are high risk. 
  if(!exists('high_risk_list')) high_risk_list <- lapply(g_name,function(x)c())
  high_risk <- high_risk_list[[first_infected]]
  ## contacts of contacts
  contacts_of_contacts <- contact_of_contact_list[[first_infected]]
  contacts_of_contacts <- contacts_of_contacts[contacts_of_contacts!=first_infected]
  ## add households of high-risk contacts to contacts of contacts
  if(length(high_risk)>0) 
    for(hr in high_risk)
      contacts_of_contacts <- c(contacts_of_contacts,household_list[[hr]])
  #high_risk <- c(high_risk,household_list[[first_infected]])
  cluster_people <- funique(c(contacts,contacts_of_contacts))
  cluster_people_index <- g_name%in%cluster_people
  
  # enroll trial participants
  n_trial_participants <- rbinom(1,length(cluster_people),enrollment_rate)
  trial_participants <- sample(cluster_people,n_trial_participants,replace=F)
  t_nodes[trial_participants] <- 1
  t_nodes <<- t_nodes
  vaccinees <- c()
  if(cluster_flag==0){
    prob1 <- allocation_ratio
    prob1c <- prob1^tune[demographic_index[trial_participants]]
    prob0c <- (1-prob1)^tune[demographic_index[trial_participants]]
    alloc_probs <- prob1c / (prob1c + prob0c)
    vaccinees <- trial_participants[as.logical(rbinom(length(trial_participants),1,alloc_probs))]
    #nvacc <- round(length(trial_participants)*allocation_ratio)
    #vaccinees <- trial_participants[sample.int(length(trial_participants),nvacc,replace=F)]
  }else{
    if(runif(1)<allocation_ratio)
      vaccinees <- trial_participants
  }
  if(length(vaccinees)>0)
    vaccine_incubation_times <- rgamma(length(vaccinees),shape=vacc_shape,rate=vacc_rate)
  if(individual_recruitment_times==F){
    recruitment_times <- rep(recruitment_time,n_trial_participants) + ceil_inc_time
  }else{
    recruitment_times <- sample(1:recruitment_time,n_trial_participants,replace=T) + ceil_inc_time
  }
  # roll epidemic forward one day at a time
  sim_time <- recruitment_time + end_time + ceil_inc_time
  for(time_step in 1:sim_time){
    ## vaccination given time to develop immunity
    if(length(vaccinees)>0) {
      developed <- vaccine_incubation_times<=time_step-recruitment_times[trial_participants%in%vaccinees]-inc_time
      v_nodes[vaccinees[developed]] <- 1
    }
    
    # update everyone's internal clock
    newinfectious <- newremoved <- c()
    if ((nrow(e_nodes_info)>0)||(nrow(i_nodes_info)>0)) {
      time_diff <- NULL
      if(!individual_recruitment_times) time_diff <- recruitment_time-time_step
      rec_list <- recover(e_nodes_info,i_nodes_info,infperiod_shape,infperiod_rate,cluster_people_index=cluster_people_index,time_diff=time_diff)
      e_nodes_info <- rec_list[[1]]
      i_nodes_info <- rec_list[[2]]
      newremoved <- rec_list[[3]]
      i_nodes[newremoved] <- 0
      r_nodes[newremoved] <- 1 
      newinfectious <- rec_list[[4]]
      e_nodes[newinfectious] <- 0
      i_nodes[newinfectious] <- 1
    }
    # store new cases of infectiousness
    numnewinfectious <- length(newinfectious)
    if (numnewinfectious>0) {
      # Update results
      results <- rbind(results,cbind(newinfectious,time_step,time_step-as.numeric(i_nodes_info[match(newinfectious,i_nodes_info[,1]),4]),NA))
      numinfectious <- numinfectious+numnewinfectious
    }
    if(length(newremoved)>0){
      results[results[,1]%in%newremoved,4] <- time_step
    }
    
    ## spread infection
    e_nodes_info <- spread_wrapper(i_nodes_info,s_nodes,v_nodes,e_nodes_info,direct_VE)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes[e_nodes_info[,1]] <- 1
    order_infected <- c(order_infected,e_nodes_info[,1])
    
    # infect from source
    rate_from_source <- max((start_day + time_step)*from_source + base_rate, 0)
    if(rate_from_source>0){
      e_nodes_info <- infect_from_source(s_nodes,v_nodes,e_nodes_info,direct_VE,incperiod_shape,incperiod_rate,rate_from_source)
      s_nodes[e_nodes_info[,1]] <- 0
      e_nodes[e_nodes_info[,1]] <- 1
      order_infected <- c(order_infected,e_nodes_info[,1])
    }
    
    # store information
    trajectories$S[time_step+1] <- sum(s_nodes)# + sum(v_nodes) + sum(c_nodes)
    trajectories$E[time_step+1] <- trajectories$S[time_step] - trajectories$S[time_step+1]
    trajectories$I[time_step+1] <- numnewinfectious
    trajectories$R[time_step+1] <- sum(r_nodes)-sum(trajectories$R)
    
  }
  
  # store information and format to return
  results <- as.data.frame(results)
  colnames(results) <- c('InfectedNode', 'DayInfectious', 'DayInfected', 'DayRemoved')
  results$inCluster <- results$InfectedNode%in%cluster_people
  results$contact <- results$InfectedNode%in%contacts
  #results$highrisk <- results$InfectedNode%in%high_risk
  results$inTrial <- results$InfectedNode%in%trial_participants
  results$vaccinated <- results$InfectedNode%in%vaccinees
  results$RecruitmentDay <- recruitment_times[match(results$InfectedNode,trial_participants)]
  
  return(list(results,length(cluster_people),recruitment_times,length(vaccinees),length(trial_participants),vaccinees,trial_participants,order_infected))
}


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
trial_designs <- rbind(trial_designs,trial_designs[trial_designs$adapt=='',])
trial_designs$weight[(nCombAdapt+1):(nComb*(length(adaptations)+1))] <- 'binary'
trial_designs$powertst <- trial_designs$VE_esttst <- trial_designs$VE_sdtst <- trial_designs$VE_estht <- trial_designs$VE_sdht <- 
  trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
registerDoParallel(cores=12)
eval_day <- 20
latest_infector_time <- eval_day - 0

trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
  set.seed(des)
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  if(adaptation=='TS'){
    simulate_contact_network <- simulate_contact_network_age
  }else{
    simulate_contact_network <- simulate_contact_network_original
  }
  vaccinated_count <- infectious_count <- enrolled_count <- list()
  for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- enrolled_count[[i]] <- 0
  pval_binary_mle3 <- ve_est3 <- pval_binary_mle2 <- ve_est2 <- pval_binary_mle <- ve_est <- ve_estht <- c()
  rr_list <- list()
  exports <- deaths <- c()
  for(tr in 1:nTrials){
    randomisation_ratios <- c()
    people_per_ratio <- c()
    vaccinees <- trial_participants <- c()
    infectious_by_vaccine <- excluded <- matrix(0,nrow=nClusters,ncol=2)
    results_list <- list()
    allocation_ratio <- 0.5
    netwk_list <- list()
    for(iter in 1:nClusters){
      ## select random person to start
      randomisation_ratios[iter] <- allocation_ratio
      first_infected <- sample(g_name[eligible_first_person],1)
      netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,end_time=eval_day,allocation_ratio=allocation_ratio,direct_VE=direct_VE,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper)
      netwk_list[[iter]] <- netwk
      results_list[[iter]] <- netwk[[1]]
      results <- results_list[[iter]]
      infectious_by_vaccine[iter,] <- c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+9),sum(!results$vaccinated&results$inTrial&results$DayInfectious>results$RecruitmentDay+9))
      excluded[iter,] <- c(sum(results$vaccinated&results$DayInfectious<results$RecruitmentDay+10),sum(!results$vaccinated&results$inTrial&results$DayInfectious<results$RecruitmentDay+10))
      
      vaccinees[iter] <- netwk[[4]]
      trial_participants[iter] <- netwk[[5]]
      
      ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
      if(adaptation!=''&&iter %% eval_day == 0 && sum(vaccinees)>0){
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1)
        pop_sizes2 <- probs[[2]]
        fails <- probs[[3]]
        allocation_ratio <- response_adapt(fails,pop_sizes2,days=iter,adaptation)
        people_per_ratio <- rbind(people_per_ratio,c(sum(trial_participants),iter,allocation_ratio))
        #if(allocation_ratio==0) break
      }
    }
    if(tr<6) rr_list[[tr]] <- people_per_ratio
    ## regular test
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1)
    pval_binary_mle2[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]
    ## correct VE test
    #eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,randomisation_ratios=randomisation_ratios,#rbht_norm=0,
    #                                           rbht_norm=ifelse(adaptation=='',1,2),
    #                                           people_per_ratio=people_per_ratio,adaptation=adaptation,contact_network=-1)#adaptation=adapt if rbht_norm=2
    ve_estht[tr]  <- eval_list[[1]]
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees)/nTrials
    enrolled_count[[1]] <- enrolled_count[[1]] + sum(trial_participants)/nTrials
    infectious_count[[1]] <- infectious_count[[1]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
    if(adaptation==''){
      pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
      pval_binary_mle[tr] <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      ve_est[tr]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      vaccinated_count[[2]] <- vaccinated_count[[2]] + sum(vaccinees)/nTrials
      enrolled_count[[2]] <- enrolled_count[[2]] + sum(trial_participants)/nTrials
      infectious_count[[2]] <- infectious_count[[2]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
    }
    ## if a test was done
    #eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=T,contact_network=-1)
    #pval_binary_mle3[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    #ve_est3[tr]  <- eval_list[[1]]
    ## correcting for trend 
    pval_binary_mle3[tr]  <- NA
    ve_est3[tr]  <- NA
    if(adaptation!='')
      pval_binary_mle3[tr] <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                                    tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio)
    
    ## exports
    exports[tr] <- sum(sapply(results_list,function(x)sum(!x$inCluster)-1))
    infected_nodes <- demographic_index[do.call(rbind,results_list)$InfectedNode]
    deaths[tr] <- sum(sapply(unique(demographic_index),function(x) sum(infected_nodes==x)*grouped_cfr[x]))
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
  return(list(power, VE_est, VE_sd,vaccinated_count, infectious_count, enrolled_count,rr_list,mean(exports),mean(deaths)))
}
trial_designs$mee <- 0
for(des in 1:nCombAdapt){
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  trial_designs$vaccinated[des] <- trial_results[[des]][[4]][[1]]
  trial_designs$infectious[des] <- trial_results[[des]][[5]][[1]]
  trial_designs$enrolled[des] <- trial_results[[des]][[6]][[1]]
  #if(adaptation==''){
  #  trial_designs$vaccinated[des+nComb] <- trial_results[[des]][[4]][[2]]
  #  trial_designs$infectious[des+nComb] <- trial_results[[des]][[5]][[2]]
  #  trial_designs$enrolled[des+nComb] <- trial_results[[des]][[6]][[2]]
  #}
  trial_designs$power[des] <- trial_results[[des]][[1]][1]
  trial_designs$VE_est[des] <- trial_results[[des]][[2]][1]
  trial_designs$VE_sd[des] <- trial_results[[des]][[3]][1]
  trial_designs$powertst[des] <- trial_results[[des]][[1]][3]
  trial_designs$VE_esttst[des] <- trial_results[[des]][[2]][3]
  trial_designs$VE_sdtst[des] <- trial_results[[des]][[3]][3]
  trial_designs$VE_estht[des] <- trial_results[[des]][[2]][4]
  trial_designs$VE_sdht[des] <- trial_results[[des]][[3]][4]
  trial_designs$mee[des] <- trial_results[[des]][[8]]
  trial_designs$deaths[des] <- trial_results[[des]][[9]]
  #if(adaptation==''){
  #  trial_designs$power[des+nComb] <- trial_results[[des]][[1]][2]
  #  trial_designs$VE_est[des+nComb] <- trial_results[[des]][[2]][2]
  #  trial_designs$VE_sd[des+nComb] <- trial_results[[des]][[3]][2]
  #  trial_designs$powertst[des+nComb] <- trial_results[[des]][[1]][3]
  #  trial_designs$VE_esttst[des+nComb] <- trial_results[[des]][[2]][3]
  #  trial_designs$VE_sdtst[des+nComb] <- trial_results[[des]][[3]][3]
  #  trial_designs$VE_estht[des+nComb] <- trial_results[[des]][[2]][4]
  #  trial_designs$VE_sdht[des+nComb] <- trial_results[[des]][[3]][4]
  #  trial_designs$mee[des+nComb] <- trial_results[[des]][[8]]
  #}
}
subset(trial_designs,VE==0)
subset(trial_designs,VE>0)
saveRDS(trial_designs,'storage/silo_trials.Rds')
result_table <- subset(trial_designs,VE>0)[,c(3:15)]
result_table$t1e <- subset(trial_designs,VE==0)$power
result_table$t1etst <- subset(trial_designs,VE==0)$powertst
result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
result_table <- result_table[,!colnames(result_table)%in%c('VE_est','VE_sd')]
result_table$htVE <- paste0(round(result_table$VE_estht,2),' (',round(result_table$VE_sdht,2),')')
result_table <- result_table[,!colnames(result_table)%in%c('VE_estht','VE_sdht')]
#result_table$tstVE <- paste0(round(result_table$VE_esttst,2),' (',round(result_table$VE_sdtst,2),')')
result_table <- result_table[,!colnames(result_table)%in%c('VE_esttst','VE_sdtst')]
result_table$adapt <- as.character(result_table$adapt)
result_table$adapt[result_table$adapt==''] <- 'None'
result_table$nmee <- subset(trial_designs,VE==0)$mee - subset(trial_designs,VE>0)$mee
result_table$averted <- subset(trial_designs,VE==0)$deaths - subset(trial_designs,VE>0)$deaths
#result_table$cluster[result_table$cluster==0] <- 'Individual'
#result_table$cluster[result_table$cluster==1] <- 'Cluster'
colnames(result_table) <- c('Adaptation','Weighting','Sample size','Infectious','Vaccinated','Power','Power (corrected)',
                            'Type 1 error','Type 1 error (corrected)','VE estimate','VE estimate (TH)','NMEE','Deaths averted')
print(xtable(result_table), include.rownames = FALSE)



