source('set_up_script.R')
registerDoParallel(cores=12)

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
vaccine_efficacies <- c(0,0.7)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'binary'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
trial_designs <- rbind(trial_designs,trial_designs[trial_designs$adapt=='',])
trial_designs$weight[(nCombAdapt+1):(nComb*(length(adaptations)+1))] <- 'binary'
trial_designs$powertst <- trial_designs$VE_esttst <- trial_designs$VE_sdtst <- trial_designs$VE_estht <- trial_designs$VE_sdht <- 
  trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
func <- function(results_list,vaccinees,trial_participants,infectious_by_vaccine,excluded,max_time=10000,contact_network=2){
  pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
  pval_binary_mle <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
  ve_estimate  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
  weight_sums <- colSums(infectious_by_vaccine)
  return(list(ve_estimate[1],pop_sizes,weight_sums))
}

latest_infector_time <- eval_day - 0
trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
  set.seed(des)
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  vaccinated_count <- infectious_count <- enrolled_count <- list()
  for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- enrolled_count[[i]] <- 0
  pval_binary_mle3 <- ve_est3 <- pval_binary_mle2 <- ve_est2 <- pval_binary_mle <- ve_est <- ve_estht <- c()
  rr_list <- list()
  exports <- c()
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
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
        pop_sizes2 <- probs[[2]]
        fails <- probs[[3]]
        allocation_ratio <- response_adapt(fails,pop_sizes2,days=iter,adaptation)
        people_per_ratio <- rbind(people_per_ratio,c(sum(trial_participants),iter,allocation_ratio))
        #if(allocation_ratio==0) break
      }
    }
    rr_list[[tr]] <- people_per_ratio
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees)/nTrials
    enrolled_count[[1]] <- enrolled_count[[1]] + sum(trial_participants)/nTrials
    infectious_count[[1]] <- infectious_count[[1]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
    ## correcting for trend 
    if(adaptation!='')
      pval_binary_mle3[tr] <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                                    tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio,observed=observed)
    
    ## exports
    exports[tr] <- sum(sapply(results_list,function(x)sum(!x$inCluster)-1))
    
    eval_list <- func(results_list,vaccinees,trial_participants,infectious_by_vaccine,excluded)
    pval_binary_mle2[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]
    if(adaptation==''){
      pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
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
  power[3] <- sum(pval_binary_mle3<0.05,na.rm=T)/sum(!is.na(pval_binary_mle3))
  VE_est[1] <- mean(ve_est2,na.rm=T)
  VE_est[3] <- mean(ve_est3,na.rm=T)
  VE_sd[1] <- sd(ve_est2,na.rm=T)
  VE_sd[3] <- sd(ve_est3,na.rm=T)
  if(adaptation==''){
    power[2] <- sum(pval_binary_mle<0.05,na.rm=T)/sum(!is.na(pval_binary_mle))
    VE_est[2] <- mean(ve_est,na.rm=T)
    VE_sd[2] <- sd(ve_est,na.rm=T)
  }
  return(list(power, VE_est, VE_sd,vaccinated_count, infectious_count, enrolled_count))
}
saveRDS(trial_results,'storage/bin_trial_results.Rds')
trial_results <- readRDS('storage/bin_trial_results.Rds')
for(des in 1:nCombAdapt){
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  trial_designs$vaccinated[des] <- trial_results[[des]][[4]][[1]]
  trial_designs$infectious[des] <- trial_results[[des]][[5]][[1]]
  trial_designs$enrolled[des] <- trial_results[[des]][[6]][[1]]
  trial_designs$power[des] <- trial_results[[des]][[1]][1]
  trial_designs$VE_est[des] <- trial_results[[des]][[2]][1]
  trial_designs$VE_sd[des] <- trial_results[[des]][[3]][1]
  trial_designs$ttepower[des] <- trial_results[[des]][[1]][3]
  trial_designs$tteVE_est[des] <- trial_results[[des]][[2]][3]
  trial_designs$tteVE_sd[des] <- trial_results[[des]][[3]][3]
}
subset(trial_designs,VE==0)
subset(trial_designs,VE>0)

result_table <- subset(trial_designs,VE>0)[,c(2:13)[-c(10:12)]]
result_table_tte <- subset(trial_designs,VE>0)[,c(2:13)[-c(7:9)]]
result_table$t1e <- subset(trial_designs,VE==0)$power
result_table_tte$t1e <- subset(trial_designs,VE==0)$ttepower
result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
result_table_tte$VE <- paste0(round(result_table_tte$tteVE_est,2),' (',round(result_table_tte$tteVE_sd,2),')')
result_table <- result_table[,!colnames(result_table)%in%c('VE_est','VE_sd')]
result_table_tte <- result_table_tte[,!colnames(result_table_tte)%in%c('tteVE_est','tteVE_sd')]
colnames(result_table_tte)[colnames(result_table_tte)=='ttepower'] <- 'power'
result_table$endpoint <- 'binary'
result_table_tte$endpoint <- 'TTE'
result_table_tte$weight <- 'continuous'
result_table <- rbind(result_table,result_table_tte)
result_table$adapt <- as.character(result_table$adapt)
result_table$adapt[result_table$adapt==''] <- 'None'
result_table$cluster[result_table$cluster==0] <- 'Individual'
#result_table$cluster[result_table$cluster==1] <- 'Cluster'
result_table <- result_table[,c(1:3,10,4:9)]
#result_table <- subset(result_table,!(endpoint=='TTE'&weight=='binary'))
colnames(result_table) <- c('Randomisation','Adaptation','Weighting','Endpoint','Sample size','Infectious','Vaccinated','Power','Type 1 error','VE estimate')
print(xtable(result_table), include.rownames = FALSE)

saveRDS(result_table,'storage/binsilo.Rds')