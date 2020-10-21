source('set_up_script.R')
registerDoParallel(cores=8)

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
vaccine_efficacies <- c(0,0.7)
adaptations <- ''#c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,weight=c('none','binary','cont','continuous'),stringsAsFactors = F)
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
ref_recruit_day <- 30
get_efficacious_probabilities_orig <- get_efficacious_probabilities
get_infectee_weights_orig <- get_infectee_weights


latest_infector_time <- eval_day - 0
trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
  set.seed(des)
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  weight <- trial_designs$weight[des]
  if(weight=='none'){
    get_efficacious_probabilities <- get_efficacious_probabilities_none
    get_infectee_weights <- get_infectee_weights_none
  }else if(weight=='binary'){
    get_efficacious_probabilities <- get_efficacious_probabilities_bin
    get_infectee_weights <- get_infectee_weights_bin
    target_weight <- 26
  }else if(weight=='cont'){
    get_efficacious_probabilities <- get_efficacious_probabilities_cont
    get_infectee_weights <- get_infectee_weights_orig
  }else if(weight=='continuous'){
    get_efficacious_probabilities <- get_efficacious_probabilities_orig
    get_infectee_weights <- get_infectee_weights_orig
    target_weight <- 25
  }
  vaccinated_count <- rr_list <- list()
  for(i in 1:2) vaccinated_count[[i]] <- 0
  pval_binary_mle3 <- ve_est3 <- zval_binary_mle2 <- ve_est2 <- zval_binary_mle <- ve_est <- ve_estht <- c()
  exports <- enrolled_count <- infectious_count <- c()
  for(tr in 1:nTrials){
    randomisation_ratios <- c()
    people_per_ratio <- c()
    vaccinees <- trial_participants <- c()
    results_list <- list()
    allocation_ratio <- 0.5
    netwk_list <- list()
    weight_break <- 0
    iter <- 0
    while(weight_break<target_weight){
      iter <- iter + 1
      set.seed(iter*nTrials+tr)
      #for(iter in 1:nClusters){
      ## select random person to start
      randomisation_ratios[iter] <- allocation_ratio
      first_infected <- sample(g_name[eligible_first_person],1)
      netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,end_time=eval_day,allocation_ratio=allocation_ratio,direct_VE=direct_VE,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper)
      netwk[[1]]$seroconverted <- netwk[[9]][match(netwk[[1]]$InfectedNode,netwk[[6]])]
      netwk_list[[iter]] <- netwk
      results_list[[iter]] <- netwk[[1]]
      results_list[[iter]]$obs <- c(runif(nrow(results_list[[iter]]))<observed)
      results <- results_list[[iter]]
      vaccinees[iter] <- netwk[[4]]
      trial_participants[iter] <- netwk[[5]]
      ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
      if(sum(vaccinees)>0){
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
        weight_break <- sum(probs[[3]])
      }
    }
    rr_list[[tr]] <- people_per_ratio
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees)/nTrials
    enrolled_count[tr] <- sum(trial_participants)
    infectious_count[tr] <- sum(sapply(results_list,function(x)sum(x$obs&x$inTrial&!is.na(x$DayInfectious)&x$RecruitmentDay<x$DayInfectious)))
    ## correcting for trend 
    if(adaptation!='')
      pval_binary_mle3[tr] <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                                    tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio,observed=observed,eval_ve=F)
    
    ## exports
    exports[tr] <- sum(sapply(results_list,function(x)sum(!x$inCluster)-1))/length(results_list)*100
    
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
    zval_binary_mle2[tr]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]

    eval_list <- get_efficacious_probabilities_tte(results_list,vaccinees,trial_participants)
    zval_binary_mle[tr]  <- eval_list[[2]]
    ve_est[tr]  <- eval_list[[1]]
  }
  power <- VE_est <- VE_sd <- c()
  power[1] <- sum(zval_binary_mle2>qnorm(0.95),na.rm=T)/sum(!is.na(zval_binary_mle2))
  power[2] <- sum(zval_binary_mle>qnorm(0.95),na.rm=T)/sum(!is.na(zval_binary_mle))
  power[3] <- sum(zval_binary_mle2<pval_binary_mle3,na.rm=T)/sum(!is.na(pval_binary_mle3)&!is.na(zval_binary_mle2))
  print(c(des,adaptation,power))
  VE_est[1] <- mean(ve_est2,na.rm=T)
  VE_est[2] <- mean(ve_est,na.rm=T)
  VE_est[3] <- mean(ve_est2[!is.na(zval_binary_mle2)&zval_binary_mle2>qnorm(0.95)],na.rm=T)
  VE_sd[1] <- sd(ve_est2,na.rm=T)
  VE_sd[2] <- sd(ve_est,na.rm=T)
  VE_sd[3] <- sd(ve_est2[!is.na(zval_binary_mle2)&zval_binary_mle2>qnorm(0.95)],na.rm=T)
  if(adaptation==''){
    VE_est[2] <- mean(ve_est,na.rm=T)
    VE_sd[2] <- sd(ve_est,na.rm=T)
  }
  enrolled <- list(mean(enrolled_count),sd(enrolled_count))
  print(c(des,power))
  #saveRDS(infectious_count,paste0('storage/inf',des,'.Rds'))
  return(list(power, VE_est, VE_sd,vaccinated_count, mean(infectious_count), enrolled,rr_list,mean(exports)))
}
saveRDS(trial_results,'storage/bin_trial_results.Rds')
trial_results <- readRDS('storage/bin_trial_results.Rds')
trial_designs$prange <- trial_designs$mee <- trial_designs$powertst <- trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$power2 <- trial_designs$VE_est2 <- trial_designs$VE_sd2 <- 
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
  trial_designs$power2[des] <- trial_results[[des]][[1]][2]
  trial_designs$prange[des] <- trial_results[[des]][[1]][2]
  trial_designs$VE_est[des] <- trial_results[[des]][[2]][1]
  trial_designs$VE_sd[des] <- trial_results[[des]][[3]][1]
  trial_designs$VE_est2[des] <- trial_results[[des]][[2]][2]
  trial_designs$VE_sd2[des] <- trial_results[[des]][[3]][2]
  trial_designs$powertst[des] <- trial_results[[des]][[1]][3]
  trial_designs$mee[des] <- trial_results[[des]][[8]]
}
subset(trial_designs,VE==0)
subset(trial_designs,VE>0)
result_table <- subset(trial_designs,VE>0)[,-c(1:2)]
result_table$t1e <- subset(trial_designs,VE==0)$power
#result_table$t1etst <- subset(trial_designs,VE==0)$powertst
result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
#result_table$power <- paste0(round(result_table$power,2),' (',round(result_table$prange,2),')')
result_table$enrolled <- paste0(result_table$enrolled,' (',result_table$enrolledsd,')')
result_table$nullenrolled <- paste0(subset(trial_designs,VE==0)$enrolled,' (',subset(trial_designs,VE==0)$enrolledsd,')')
result_table$adapt <- as.character(result_table$adapt)
result_table$adapt[result_table$adapt==''] <- 'None'
result_table$nmee <- subset(trial_designs,VE==0)$mee - subset(trial_designs,VE>0)$mee
result_table <- result_table[,!colnames(result_table)%in%c('powertst','adapt','VE_est','VE_sd','enrolledsd','mee','prange')]
colnames(result_table) <- c('Weighting','Sample size','Symptomatic','Vaccinated','Power',
                            'Type 1 error','VE estimate','Null enrolled','Prevented export infections')
print(xtable(result_table,digits=c(0,0,0,0,0,2,2,0,0,2)), include.rownames = FALSE)

saveRDS(result_table,'storage/binsilo.Rds')


