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
get_efficacious_probabilities_bin <- function(results_list,vaccinees,trial_participants,max_time=10000,contact_network=2,
                                              tested=F,randomisation_ratios=NULL,rbht_norm=0,people_per_ratio=NULL,adaptation='TST',observed=1,age_counts=NULL){
  infectious_by_vaccine <- excluded <- c()
  for(iter in 1:length(results_list)){
    results <- results_list[[iter]]
    infectious_by_vaccine <- rbind(infectious_by_vaccine,
                                   c(sum((results$vaccinated&results$DayInfectious>results$RecruitmentDay+8)*c(runif(nrow(results))<observed)),
                                     sum((!results$vaccinated&results$inTrial&results$DayInfectious>results$RecruitmentDay+8)*c(runif(nrow(results))<observed))))
    excluded <- rbind(excluded,c(sum(results$vaccinated&results$DayInfectious<results$RecruitmentDay+9),
                                 sum(!results$vaccinated&results$inTrial&results$DayInfectious<results$RecruitmentDay+9)))
  }
  weight_sums <- colSums(infectious_by_vaccine,na.rm=T)
  pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
  
  pval_binary_mle <- calculate_pval(weight_sums,pop_sizes)
  ve_estimate  <- calculate_ve(weight_sums,pop_sizes)
  
  return(list(ve_estimate[1],pop_sizes,weight_sums))
}
get_infectee_weights_bin <- function(results,ve_point_est,contact_network=2,tested=F,correct_for_ve=T){
  # the day the cluster is recruited
  recruit_day <- results$RecruitmentDay
  # the day individuals became infectious
  days_infectious <- results$DayInfectious
  # the durations for which they were infectious
  weight_hh_rem <- matrix(0,ncol=2,nrow=1)
  infectee_names <- c()
  # those who were infected by someone else
  infectee_index <- !is.na(recruit_day) & days_infectious>recruit_day
  if(sum(infectee_index)>0){
    weight_hh_rem <- matrix(0,ncol=2,nrow=sum(infectee_index))
    infectees <- days_infectious[infectee_index]
    infectee_names <- results$InfectedNode[infectee_index]
    infectee_trial <- results$inTrial[infectee_index]
    infectee_vaccinated <- results$vaccinated[infectee_index]
    for(j in 1:length(infectees)){
      if(infectee_trial[j]&(days_infectious[infectee_index]>recruit_day[infectee_index]+8)[j]){
        # add to weight for vaccinated or unvaccinated
        if(infectee_vaccinated[j]){
          weight_hh_rem[j,1] <- 1
        }else{
          weight_hh_rem[j,2] <- 1
        }
      }
    }
  }
  return(list(weight_hh_rem,infectee_names))
}
get_efficacious_probabilities_none <- function(results_list,vaccinees,trial_participants,max_time=10000,contact_network=2,
                                              tested=F,randomisation_ratios=NULL,rbht_norm=0,people_per_ratio=NULL,adaptation='TST',observed=1,age_counts=NULL){
  infectious_by_vaccine <- excluded <- c()
  for(iter in 1:length(results_list)){
    results <- results_list[[iter]]
    infectious_by_vaccine <- rbind(infectious_by_vaccine,
                                   c(sum((results$vaccinated)*c(runif(nrow(results))<observed)),
                                     sum((!results$vaccinated&results$inTrial)*c(runif(nrow(results))<observed))))
  }
  weight_sums <- colSums(infectious_by_vaccine,na.rm=T)
  pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees))
  
  pval_binary_mle <- calculate_pval(weight_sums,pop_sizes)
  ve_estimate  <- calculate_ve(weight_sums,pop_sizes)
  
  return(list(ve_estimate[1],pop_sizes,weight_sums))
}
get_infectee_weights_none <- function(results,ve_point_est,contact_network=2,tested=F,correct_for_ve=T){
  # the day the cluster is recruited
  recruit_day <- results$RecruitmentDay
  # the day individuals became infectious
  days_infectious <- results$DayInfectious
  # the durations for which they were infectious
  weight_hh_rem <- matrix(0,ncol=2,nrow=1)
  infectee_names <- c()
  # those who were infected by someone else
  infectee_index <- !is.na(recruit_day) & days_infectious>recruit_day
  if(sum(infectee_index)>0){
    weight_hh_rem <- matrix(0,ncol=2,nrow=sum(infectee_index))
    infectees <- days_infectious[infectee_index]
    infectee_names <- results$InfectedNode[infectee_index]
    infectee_trial <- results$inTrial[infectee_index]
    infectee_vaccinated <- results$vaccinated[infectee_index]
    for(j in 1:length(infectees)){
      if(infectee_trial[j]){
        # add to weight for vaccinated or unvaccinated
        if(infectee_vaccinated[j]){
          weight_hh_rem[j,1] <- 1
        }else{
          weight_hh_rem[j,2] <- 1
        }
      }
    }
  }
  return(list(weight_hh_rem,infectee_names))
}
get_efficacious_probabilities_cont <- function(results_list,vaccinees,trial_participants,max_time=10000,contact_network=2,
                                          tested=F,randomisation_ratios=NULL,rbht_norm=0,people_per_ratio=NULL,adaptation='TST',observed=1,age_counts=NULL){
  controls <- trial_participants - vaccinees
  if(is.null(randomisation_ratios)) randomisation_ratios <- rep(0.5,length(trial_participants))
  
  uninf_vacc <- vaccinees - sapply(results_list,function(x)sum(x$vaccinated))
  uninf_cont <- trial_participants - vaccinees - sapply(results_list,function(x)sum(x$inTrial&!x$vaccinated))
  
  uninf <- data.frame(vaccinated=c(rep(T,sum(uninf_vacc)),rep(F,sum(uninf_cont))),
                      allocRatio=c(rep(randomisation_ratios,uninf_vacc),rep(randomisation_ratios,uninf_cont)),
                      weight=1,infected=F)
  
  ve_estimate <- c(0,1)
  break_count <- 0
  not_nas <- lapply(1:length(results_list),function(x){
    results <- results_list[[x]]
    !is.na(results$RecruitmentDay)&results$RecruitmentDay<results$DayInfectious # subset(y,RecruitmentDay<DayInfectious)
  })
  if(contact_network==-1){
    results_tab_list <- list()
    for(x in 1:length(results_list)){
      results <- results_list[[x]]
      results$startDay <- x
      results$allocRatio <- randomisation_ratios[x]
      results_tab_list[[x]] <- results
    }
    result_tab <- do.call(rbind,results_tab_list)
    result_tab <- result_tab[unlist(not_nas),]
    if(nrow(result_tab)>0) {
      result_tab$infected <- T
      #result_tab$observed <- runif(nrow(result_tab))<observed
    }
  }
  
  if(nrow(result_tab)>0){
    weights <- get_infectee_weights(results=result_tab,ve_point_est=ve_estimate[1],contact_network,tested,correct_for_ve=F)
    result_tab$weight <- rowSums(weights[[1]])*c(runif(nrow(result_tab))<observed)
  }
  
  #result_tab$weight <- rowSums(get_infectee_weights(result_tab,ve_estimate[1],contact_network,tested)[[1]])
  ve_estimate[2] <- ve_estimate[1]
  
  if(nrow(result_tab)==0) return(list(0,c(sum(vaccinees),sum(trial_participants)-sum(vaccinees)),c(0,0)))
  
  all_results <- rbind(result_tab[,match(colnames(uninf),colnames(result_tab))],uninf)
  weights <- get_weights_from_all_results(all_results)
  fails <- weights[[1]]
  pop_sizes2 <- weights[[2]]
  if(fails[2]>0&&!any(pop_sizes2==0))
    ve_estimate[1] <- calculate_ve(fails,sizes=pop_sizes2)
  break_count <- break_count + 1
  
  return(list(ve_estimate[1],pop_sizes2,fails))
}

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
  }else if(weight=='cont'){
    get_efficacious_probabilities <- get_efficacious_probabilities_cont
    get_infectee_weights <- get_infectee_weights_orig
  }else if(weight=='continuous'){
    get_efficacious_probabilities <- get_efficacious_probabilities_orig
    get_infectee_weights <- get_infectee_weights_orig
  }
  vaccinated_count <- infectious_count <- rr_list <- list()
  for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- 0
  pval_binary_mle3 <- ve_est3 <- zval_binary_mle2 <- ve_est2 <- pval_binary_mle <- ve_est <- ve_estht <- c()
  exports <- enrolled_count <- c()
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
      results <- results_list[[iter]]
      
      vaccinees[iter] <- netwk[[4]]
      trial_participants[iter] <- netwk[[5]]
      
      
      ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
      if(iter >= eval_day && sum(vaccinees)>0){
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
        weight_break <- sum(probs[[3]])
      }
    }
    rr_list[[tr]] <- people_per_ratio
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees)/nTrials
    enrolled_count[tr] <- sum(trial_participants)
    infectious_count[[1]] <- infectious_count[[1]] + (observed*sum(sapply(results_list,function(x)sum(x$inTrial&!is.na(x$DayInfectious)))))/nTrials
    ## correcting for trend 
    if(adaptation!='')
      pval_binary_mle3[tr] <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                                    tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio,observed=observed,eval_ve=F)
    
    ## exports
    exports[tr] <- sum(sapply(results_list,function(x)sum(!x$inCluster)-1))/length(results_list)*100
    
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
    zval_binary_mle2[tr]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]

    #eval_list <- get_efficacious_probabilities_orig(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
    #pval_binary_mle[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    #ve_est[tr]  <- eval_list[[1]]
  }
  power <- VE_est <- VE_sd <- c()
  power[1] <- sum(zval_binary_mle2>qnorm(0.95),na.rm=T)/sum(!is.na(zval_binary_mle2))
  power[3] <- sum(zval_binary_mle2<pval_binary_mle3,na.rm=T)/sum(!is.na(pval_binary_mle3)&!is.na(zval_binary_mle2))
  print(c(des,adaptation,power))
  VE_est[1] <- mean(ve_est2,na.rm=T)
  VE_est[3] <- mean(ve_est3,na.rm=T)
  VE_sd[1] <- sd(ve_est2,na.rm=T)
  VE_sd[3] <- sd(ve_est3,na.rm=T)
  if(adaptation==''){
    VE_est[2] <- mean(ve_est,na.rm=T)
    VE_sd[2] <- sd(ve_est,na.rm=T)
  }
  power[2] <- infotheo::entropy(discretize(zval_binary_mle2,disc='equalwidth')) #quantile(zval_binary_mle2,0.95) - quantile(zval_binary_mle2,0.05)
  enrolled <- list(mean(enrolled_count),sd(enrolled_count))
  print(c(des,power))
  print(enrolled)
  if(des%%2==0&tr==1){
    true_excluded <- colSums(sapply(c(T,F),function(y)sapply(results_list,function(x)sum(x$vaccinated==y&(x$seroconverted>=x$DayInfected|x$RecruitmentDay>=x$DayInfected),na.rm=T))))
    true_excluded
    w <- list()
    w[[1]] <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
    vaxes <- sum(vaccinees) - sapply(w,function(x)x[[2]][1])
    cont <- sum(trial_participants) - sum(vaccinees) - sapply(w,function(x)x[[2]][2])
    ss <- sum(trial_participants) - sum(sapply(w,function(x)x[[2]]))
    ss[2] <- sum(true_excluded)
    vaxes[2] <- true_excluded[1]
    cont[2] <- true_excluded[2]
    data.frame(ss,vaxes,cont)
    print(xtable(data.frame(c(weight,'True'),(ss),(vaxes),(cont)),digits=0), include.rownames = FALSE)
  }
  return(list(power, VE_est, VE_sd,vaccinated_count, infectious_count, enrolled,rr_list,mean(exports)))
}
saveRDS(trial_results,'storage/bin_trial_results.Rds')
trial_results <- readRDS('storage/bin_trial_results.Rds')
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
                            'Type 1 error','VE estimate','Null enrolled','NMEE')
print(xtable(result_table,digits=c(0,0,0,0,0,2,2,0,0,2)), include.rownames = FALSE)

saveRDS(result_table,'storage/binsilo.Rds')


