source('set_up_script.R')

## ring vaccination trial ##################################################
nClusters <- 75
nTrials <- 1000
vaccine_efficacies <- c(0,0.6)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
trial_designs <- rbind(trial_designs,trial_designs[trial_designs$adapt=='',])
trial_designs$weight[(nCombAdapt+1):(nComb*(length(adaptations)+1))] <- 'binary'
trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
registerDoParallel(cores=16)
func <- get_efficacious_probabilities
eval_day <- 20
latest_infector_time <- eval_day - 0

  trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
      func <- get_efficacious_probabilities
    cluster_flag <- trial_designs$cluster[des]
    direct_VE <- trial_designs$VE[des]
    adaptation <- trial_designs$adapt[des]
    vaccinated_count <- infectious_count <- enrolled_count <- list()
    for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- enrolled_count[[i]] <- 0
    pval_binary_mle3 <- ve_est3 <- pval_binary_mle2 <- ve_est2 <- pval_binary_mle <- ve_est <- c()
    for(tr in 1:nTrials){
      vaccinees2 <- trial_participants2 <- c()
      infectious_by_vaccine <- excluded <- matrix(0,nrow=nClusters,ncol=2)
      results_list <- list()
      allocation_ratio <- 0.5
      netwk_list <- list()
      for(iter in 1:nClusters){
        ## select random person to start
        first_infected <- sample(g_name,1)
        netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,cluster_flag=cluster_flag,allocation_ratio=allocation_ratio,direct_VE=direct_VE)
        netwk_list[[iter]] <- netwk
        results_list[[iter]] <- netwk[[1]]
        results <- results_list[[iter]]
        infectious_by_vaccine[iter,] <- c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+9),sum(!results$vaccinated&results$inTrial&results$DayInfectious>results$RecruitmentDay+9))
        excluded[iter,] <- c(sum(results$vaccinated&results$DayInfectious<results$RecruitmentDay+10),sum(!results$vaccinated&results$inTrial&results$DayInfectious<results$RecruitmentDay+10))
        
        vaccinees2[iter] <- netwk[[4]]
        trial_participants2[iter] <- netwk[[5]]
        
        ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
        if(adaptation!=''&&iter %% eval_day == 0 && sum(vaccinees2)>0){
          allocation_ratio <- response_adapt(results_list,vaccinees2,trial_participants2,adaptation,func=func)
        }
      }
      
      
      eval_list <- func(results_list,vaccinees2,trial_participants2)
      
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
    }
    print(c(des,adaptation))
    power <- VE_est <- VE_sd <- c()
    power[1] <- sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
    VE_est[1] <- mean(ve_est2,na.rm=T)
    VE_sd[1] <- sd(ve_est2,na.rm=T)
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
    }
  }
  subset(trial_designs,VE==0)
  subset(trial_designs,VE>0)
  
  result_table <- subset(trial_designs,VE>0)[,c(2:10)]
  result_table$t1e <- subset(trial_designs,VE==0)$power
  result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
  result_table <- result_table[,!colnames(result_table)%in%c('VE_est','VE_sd')]
  result_table$adapt <- as.character(result_table$adapt)
  result_table$adapt[result_table$adapt==''] <- 'None'
  result_table$cluster[result_table$cluster==0] <- 'Individual'
  result_table$cluster[result_table$cluster==1] <- 'Cluster'
  colnames(result_table) <- c('Randomisation','Adaptation','Weighting','Sample size','Infectious','Vaccinated','Power','Type 1 error','VE estimate')
  print(xtable(result_table), include.rownames = FALSE)

