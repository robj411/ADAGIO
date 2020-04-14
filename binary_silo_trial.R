source('set_up_script.R')

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
vaccine_efficacies <- c(0,0.8)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'binary'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- length(vaccine_efficacies)*length(adaptations)
trial_designs$ttepower <- trial_designs$tteVE_est <- trial_designs$tteVE_sd <- trial_designs$power <- 
  trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
registerDoParallel(cores=10)
func <- function(results_list,vaccinees,trial_participants,infectious_by_vaccine,excluded,max_time=10000,contact_network=2){
  pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
  pval_binary_mle <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
  ve_estimate  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
  weight_sums <- colSums(infectious_by_vaccine)
  return(list(ve_estimate[1],pop_sizes,weight_sums))
}

eval_day <- 31
latest_infector_time <- eval_day - 0
rnd = 1
#for(rnd in 1:2){
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
        vaccinees2[iter] <- netwk[[4]]
        trial_participants2[iter] <- netwk[[5]]
        
        ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
        if(adaptation!=''&&iter %% eval_day == 0 && sum(vaccinees)>0){
          weights <- func(results_list,vaccinees,trial_participants,infectious_by_vaccine,excluded,max_time=10000,contact_network=2)
          allocation_ratio <- response_adapt(weights[[3]],weights[[2]],days=iter,adaptation=adaptation)
        }
      }
      
      ph_results <- iterate_ph_model(netwk_list,cluster_flag=cluster_flag,pre_randomisation=F)
      
      pval_binary_mle3[tr]  <- ph_results[1]
      ve_est3[tr]  <- ph_results[2]
      
      
      #if(rnd==1){
        eval_list <- func(results_list,vaccinees2,trial_participants2,infectious_by_vaccine,excluded)
      #}else{
      #  eval_list <- func(results_list,vaccinees,trial_participants)
      #}
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
  saveRDS(trial_results,'trial_results.Rds')
  trial_results <- readRDS('trial_results.Rds')
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
  result_table_tte$weight <- 'continous'
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
#}
