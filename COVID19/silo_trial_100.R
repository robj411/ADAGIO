source('set_up_script.R')
registerDoParallel(cores=12)

## ring vaccination trial ##################################################
nClusters <- floor(target_weight*67/32)+10
nTrials <- 1000
vaccine_efficacies <- c(0,0.7)
estimates <- c('on','under','over')
cluster_flags <- 0
adaptations <- ''
endings <- c('case','cluster')
trial_designs <- expand.grid(VE=vaccine_efficacies,estimate=estimates,adapt='',ending=endings,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
ref_recruit_day <- 30
latest_infector_time <- eval_day - 0

trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
  set.seed(des)
  cluster_flag <- 0
  direct_VE <- trial_designs$VE[des]
  estimate <- trial_designs$estimate[des]
  ending <- trial_designs$ending[des]
  if(estimate=='under') beta_base <<- beta_base*4/5
  if(estimate=='over') beta_base <<- beta_base*5/4
  adaptation <- ''
  vaccinated_count <- rr_list <- list()
  for(i in 1:2) vaccinated_count[[i]] <- 0
  pval_binary_mle3 <- ve_est3 <- zval_binary_mle2 <- ve_est2 <- pval_binary_mle <- ve_est <- ve_estht <- c()
  exports <- enrolled_count <- infectious_count <- c()
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
    while((ending=='case'&weight_break<target_weight)|(ending=='cluster'&iter<nClusters)){
      iter <- iter + 1
      set.seed(iter*nTrials+tr)
      #for(iter in 1:nClusters){
      set.seed(iter*nTrials+tr)
      ## select random person to start
      randomisation_ratios[iter] <- allocation_ratio
      first_infected <- sample(g_name[eligible_first_person],1)
      netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,end_time=eval_day,allocation_ratio=allocation_ratio,direct_VE=direct_VE,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper)
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
      }else if(iter >= eval_day && sum(vaccinees)>0 && ending=='case'){
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
        weight_break <- sum(probs[[3]])
      }
    }
    rr_list[[tr]] <- people_per_ratio
    ## regular test
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1,observed=observed)
    zval_binary_mle2[tr]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]
    ## correct VE test
    #eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,randomisation_ratios=randomisation_ratios,#rbht_norm=0,
    #                                           rbht_norm=ifelse(adaptation=='',1,2),
    #                                           people_per_ratio=people_per_ratio,adaptation=adaptation,contact_network=-1,observed=observed)#adaptation=adapt if rbht_norm=2
    ve_estht[tr]  <- eval_list[[1]]
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees)/nTrials
    enrolled_count[tr] <- sum(trial_participants)
    infectious_count[tr] <- sum(sapply(results_list,function(x)sum(x$inTrial&!is.na(x$DayInfectious)&runif(nrow(x))<observed)))
    if(adaptation==''){
      pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
      pval_binary_mle[tr] <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      ve_est[tr]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      vaccinated_count[[2]] <- vaccinated_count[[2]] + sum(vaccinees)/nTrials
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
  power[1] <- sum(zval_binary_mle2>qnorm(0.95),na.rm=T)/sum(!is.na(zval_binary_mle2))
  VE_est[1] <- mean(ve_est2,na.rm=T)
  VE_sd[1] <- sd(ve_est2,na.rm=T)
  power[3] <- sum(zval_binary_mle2<pval_binary_mle3,na.rm=T)/sum(!is.na(pval_binary_mle3)&!is.na(zval_binary_mle2))
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
  saveRDS(list(zval_binary_mle2,pval_binary_mle3),paste0('storage/silo100p',des,'.Rds'))
  enrolled <- list(mean(enrolled_count),sd(enrolled_count))
  print(des)
  print(enrolled)
  return(list(power, VE_est, VE_sd,vaccinated_count, list(mean(infectious_count),sd(infectious_count)), enrolled,rr_list,mean(exports)))
}

saveRDS(trial_results,'storage/silo_100.Rds')
trial_designs$mee <- trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- 
  trial_designs$vaccinated <- trial_designs$infectioussd <- trial_designs$infectious <- trial_designs$daysd <- trial_designs$day <- trial_designs$enrolledsd <- trial_designs$enrolled <- 0
for(des in 1:nCombAdapt){
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  trial_designs$vaccinated[des] <- round(trial_results[[des]][[4]][[1]])
  trial_designs$infectious[des] <- round(trial_results[[des]][[5]][[1]])
  trial_designs$infectioussd[des] <- round(trial_results[[des]][[5]][[2]])
  trial_designs$enrolled[des] <- round(trial_results[[des]][[6]][[1]])
  trial_designs$enrolledsd[des] <- round(trial_results[[des]][[6]][[2]])
  trial_designs$day[des] <- round(trial_results[[des]][[6]][[1]]/enrolled_per_contact)
  trial_designs$daysd[des] <- round(trial_results[[des]][[6]][[2]]/enrolled_per_contact)
  trial_designs$power[des] <- trial_results[[des]][[1]][1]
  trial_designs$VE_est[des] <- trial_results[[des]][[2]][1]
  trial_designs$VE_sd[des] <- trial_results[[des]][[3]][1]
  trial_designs$mee[des] <- trial_results[[des]][[8]]
}
subset(trial_designs,VE==0)
subset(trial_designs,VE>0)
result_table <- subset(trial_designs,VE>0)[,-c(1,3)]
result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
result_table$enrolled <- paste0(result_table$enrolled,' (',result_table$enrolledsd,')')
result_table$day <- paste0(result_table$day,' (',result_table$daysd,')')
result_table$infectious <- paste0(result_table$infectious,' (',result_table$infectioussd,')')
result_table$nmee <- subset(trial_designs,VE==0)$mee - subset(trial_designs,VE>0)$mee
result_table$enrolled2 <- subset(trial_designs,VE==0)$enrolled
result_table$enrolled2 <- paste0(round(result_table$enrolled2,2),' (',round(subset(trial_designs,VE==0)$enrolledsd,2),')')
result_table$t1e <- subset(trial_designs,VE==0)$power
result_table <- result_table[,!colnames(result_table)%in%c('daysd','infectioussd','weight','VE_est','VE_sd','enrolledsd','mee','prange')]
result_table$estimate[result_table$estimate=='on'] <- 'True'
result_table$estimate[result_table$estimate=='over'] <- 'Higher'
result_table$estimate[result_table$estimate=='under'] <- 'Lower'
result_table$ending[result_table$ending=='case'] <- 'Cases'
result_table$ending[result_table$ending=='cluster'] <- 'Sample size'
colnames(result_table) <- c('True beta','Ending','Sample size','Duration','Symptomatic','Vaccinated','Power',
                            'VE estimate','NMEE','Sample size (futile)','Type 1 error')
print(xtable(result_table,digits=c(0,0,0,0,0,0,0,2,0,2,0,2)), include.rownames = FALSE)

#saveRDS(trial_results,'storage/silo_trial_results.Rds')


