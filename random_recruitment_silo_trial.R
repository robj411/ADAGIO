source('set_up_script.R')

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
vaccine_efficacies <- c(0,0.7)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- c(0,1)
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
trial_designs <- rbind(trial_designs,trial_designs[trial_designs$adapt=='',])
trial_designs$weight[(nCombAdapt+1):(nComb*(length(adaptations)+1))] <- 'binary'
trial_designs$ttepower <- trial_designs$tteVE_est <- trial_designs$tteVE_sd <- trial_designs$power <- 
  trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
registerDoParallel(cores=20)
eval_day <- 31
latest_infector_time <- eval_day - 0

## store original lists
true_contact_of_contact_list <- contact_of_contact_list
true_high_risk_list <- high_risk_list
true_household_list <- household_list
true_contact_list <- contact_list
## get list sizes
total_size <- mean(sapply(1:length(high_risk_list),function(x){
  hrs <- high_risk_list[[x]]
  hr_hhs <- unlist(sapply(hrs,function(y)household_list[[y]]))
  c_of_c <- contact_of_contact_list[[x]]
  cont <- contact_list[[x]]
  length(unique(c(cont,c_of_c,hr_hhs)))
  }))
## overwrite with empty lists
high_risk_list <- household_list <- contact_of_contact_list <- lapply(g_name,function(x)c())
## fill contact list with random indices of mean total size
contact_list <- lapply(g_name,function(x)sample(g_name,total_size,replace=F))

ebola_spread_wrapper <- function(i_nodes_info,s_nodes,v_nodes,e_nodes_info,direct_VE){
  contact_of_contact_list <- true_contact_of_contact_list
  contact_list <- true_contact_list
  # to contacts
  current_infectious <- i_nodes_info[,1]
  if(length(current_infectious)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_list,beta_scalar=1)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=hr_and_hh_list,beta_scalar=high_risk_scalar-1)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_of_contact_list,beta_scalar=neighbour_scalar)
  }
  return(e_nodes_info)
}

for(rnd in 1){
  trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
    cluster_flag <- trial_designs$cluster[des]
    direct_VE <- trial_designs$VE[des]
    adaptation <- trial_designs$adapt[des]
    vaccinated_count <- infectious_count <- enrolled_count <- list()
    for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- enrolled_count[[i]] <- 0
    pval_binary_mle3 <- ve_est3 <- pval_binary_mle2 <- ve_est2 <- pval_binary_mle <- ve_est <- c()
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
        netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,allocation_ratio=allocation_ratio,direct_VE=direct_VE)
        netwk_list[[iter]] <- netwk
        results_list[[iter]] <- netwk[[1]]
        results <- results_list[[iter]]
        infectious_by_vaccine[iter,] <- c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+9),sum(!results$vaccinated&results$inTrial&results$DayInfectious>results$RecruitmentDay+9))
        excluded[iter,] <- c(sum(results$vaccinated&results$DayInfectious<results$RecruitmentDay+10),sum(!results$vaccinated&results$inTrial&results$DayInfectious<results$RecruitmentDay+10))
        recruit_times[iter] <- netwk[[3]][1]
        vaccinees[iter] <- netwk[[4]]
        trial_participants[iter] <- netwk[[5]]
        
        ##!! weighting non-events
        rec_day <- recruit_times[iter]
        infectious_index <- results$DayInfectious<latest_infector_time+rec_day&(results$DayRemoved>rec_day|is.na(results$DayRemoved))
        infectious_names <- results$InfectedNode[infectious_index]
        infectious_ends <- pmin(results$DayRemoved[infectious_index],latest_infector_time+rec_day)
        infectious_ends[is.na(infectious_ends)] <- latest_infector_time+rec_day
        infectious_starts <- pmax(results$DayInfectious[infectious_index],rec_day)
        vaccinees2[iter] <- netwk[[4]]
        trial_participants2[iter] <- netwk[[5]]
        
        ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
        if(adaptation!=''&&iter %% eval_day == 0 && sum(vaccinees)>0){
          weights <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network = -1)
          allocation_ratio <- response_adapt(weights[[3]],weights[[2]],days=iter,adaptation=adaptation)
        }
      }
      
      #ph_results <- iterate_ph_model(netwk_list,cluster_flag=cluster_flag,pre_randomisation=F)
      
      #pval_binary_mle3[tr]  <- ph_results[1]
      #ve_est3[tr]  <- ph_results[2]
      
      
      eval_list <- get_efficacious_probabilities(results_list,vaccinees=vaccinees2,trial_participants=trial_participants2,contact_network = -1)
      pval_binary_mle2[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
      ve_est2[tr]  <- eval_list[[1]]
      vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees2)/nTrials
      enrolled_count[[1]] <- enrolled_count[[1]] + sum(trial_participants2)/nTrials
      infectious_count[[1]] <- infectious_count[[1]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
      if(adaptation==''){
        pop_sizes <- c(sum(vaccinees2),sum(trial_participants2) - sum(vaccinees2)) - colSums(excluded)
        pval_binary_mle[tr] <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
        ve_est[tr]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
        vaccinated_count[[2]] <- vaccinated_count[[2]] + sum(vaccinees2)/nTrials
        enrolled_count[[2]] <- enrolled_count[[2]] + sum(trial_participants2)/nTrials
        infectious_count[[2]] <- infectious_count[[2]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
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
    power <- VE_est <- VE_sd <- c()
    power[1] <- sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
    #power[3] <- sum(pval_binary_mle3<0.05,na.rm=T)/sum(!is.na(pval_binary_mle3))
    VE_est[1] <- mean(ve_est2,na.rm=T)
    #VE_est[3] <- mean(ve_est3,na.rm=T)
    VE_sd[1] <- sd(ve_est2,na.rm=T)
    #VE_sd[3] <- sd(ve_est3,na.rm=T)
    if(adaptation==''){
      power[2] <- sum(pval_binary_mle<0.05,na.rm=T)/sum(!is.na(pval_binary_mle))
      VE_est[2] <- mean(ve_est,na.rm=T)
      VE_sd[2] <- sd(ve_est,na.rm=T)
    }
    return(list(power, VE_est, VE_sd,vaccinated_count, infectious_count, enrolled_count))
  }
  for(des in 1:nCombAdapt){
    cluster_flag <- trial_designs$cluster[des]
    direct_VE <- trial_designs$VE[des]
    adaptation <- trial_designs$adapt[des]
    trial_designs$vaccinated[des] <- trial_results[[des]][[4]][[1]]
    trial_designs$infectious[des] <- trial_results[[des]][[5]][[1]]
    trial_designs$enrolled[des] <- trial_results[[des]][[6]][[1]]
    if(adaptation==''){
      trial_designs$vaccinated[des+nComb] <- trial_results[[des]][[4]][[2]]
      trial_designs$infectious[des+nComb] <- trial_results[[des]][[5]][[2]]
      trial_designs$enrolled[des+nComb] <- trial_results[[des]][[6]][[2]]
    }
    trial_designs$power[des] <- trial_results[[des]][[1]][1]
    trial_designs$VE_est[des] <- trial_results[[des]][[2]][1]
    trial_designs$VE_sd[des] <- trial_results[[des]][[3]][1]
    if(adaptation==''){
      trial_designs$power[des+nComb] <- trial_results[[des]][[1]][2]
      trial_designs$VE_est[des+nComb] <- trial_results[[des]][[2]][2]
      trial_designs$VE_sd[des+nComb] <- trial_results[[des]][[3]][2]
      #trial_designs$ttepower[des+nComb] <- trial_results[[des]][[1]][3]
      #trial_designs$tteVE_est[des+nComb] <- trial_results[[des]][[2]][3]
      #trial_designs$tteVE_sd[des+nComb] <- trial_results[[des]][[3]][3]
    }
    #trial_designs$ttepower[des] <- trial_results[[des]][[1]][3]
    #trial_designs$tteVE_est[des] <- trial_results[[des]][[2]][3]
    #trial_designs$tteVE_sd[des] <- trial_results[[des]][[3]][3]
  }
  print(subset(trial_designs,VE==0))
  print(subset(trial_designs,VE>0))
  
  result_table <- subset(trial_designs,VE>0)[,c(2:13)[-c(10:12)]]
  #result_table_tte <- subset(trial_designs,VE>0)[,c(2:13)[-c(7:9)]]
  result_table$t1e <- subset(trial_designs,VE==0)$power
  #result_table_tte$t1e <- subset(trial_designs,VE==0)$ttepower
  result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
  #result_table_tte$VE <- paste0(round(result_table_tte$tteVE_est,2),' (',round(result_table_tte$tteVE_sd,2),')')
  result_table <- result_table[,!colnames(result_table)%in%c('VE_est','VE_sd')]
  #result_table_tte <- result_table_tte[,!colnames(result_table_tte)%in%c('tteVE_est','tteVE_sd')]
  #colnames(result_table_tte)[colnames(result_table_tte)=='ttepower'] <- 'power'
  result_table$endpoint <- 'binary'
  #result_table_tte$endpoint <- 'TTE'
  #result_table <- rbind(result_table,result_table_tte)
  result_table$adapt <- as.character(result_table$adapt)
  result_table$adapt[result_table$adapt==''] <- 'None'
  result_table$cluster[result_table$cluster==0] <- 'Individual'
  result_table$cluster[result_table$cluster==1] <- 'Cluster'
  result_table <- result_table[,c(1:3,10,4:9)]
  result_table <- subset(result_table,!(endpoint=='TTE'&weight=='binary'))
  colnames(result_table) <- c('Randomisation','Adaptation','Weighting','Endpoint','Sample size','Infectious','Vaccinated','Power','Type 1 error','VE estimate')
  print(xtable(result_table), include.rownames = FALSE)
  saveRDS(result_table,paste0('storage/rrsilo',rnd,'.Rds'))
}

