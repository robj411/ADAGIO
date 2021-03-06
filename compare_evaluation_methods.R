source('set_up_script.R')

## ring vaccination trial ##################################################
nClusters <- 70
nTrials <- 1000
vaccine_efficacies <- c(0,0.7)
ref_recruit_day <- 30
registerDoParallel(cores=2)
#func <- get_efficacious_probabilities
eval_day <- 31
latest_infector_time <- eval_day - 0
power <- VE_est <- VE_sd <- list()

res <- foreach(rnd = 1:2)%dopar%{
  cluster_flag <- 0
  direct_VE <- vaccine_efficacies[rnd]
  adaptation <- ''
  vaccinated_count <- infectious_count <- enrolled_count <- list()
  for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- enrolled_count[[i]] <- 0
  zval_binary_mle <- ve_est <- matrix(nrow=nTrials,ncol=9)
  for(tr in 1:nTrials){
    vaccinees <- trial_participants <- recruit_times <- c()
    vaccinees2 <- trial_participants2 <- c()
    infectious_by_vaccine <- excluded <- matrix(0,nrow=nClusters,ncol=2)
    results_list <- list()
    allocation_ratio <- 0.5
    netwk_list <- list()
    for(iter in 1:nClusters){
      #set.seed(iter*nTrials+tr)
      ## select random person to start
      first_infected <- sample(g_name,1)
      netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,allocation_ratio=allocation_ratio,direct_VE=direct_VE)
      netwk_list[[iter]] <- netwk
      results_list[[iter]] <- netwk[[1]]
      results <- results_list[[iter]]
      vax <- results$vaccinated
      too_early <- results$DayInfectious<results$RecruitmentDay+10
      infectious_by_vaccine[iter,] <- c(sum(vax&!too_early),sum(!vax&results$inTrial&!too_early))
      excluded[iter,] <- c(sum(vax&too_early),sum(!vax&results$inTrial&too_early))
      recruit_times[iter] <- max(netwk[[3]])
      
      ##!! weighting non-events
      rec_day <- recruit_times[iter]
      infectious_index <- results$DayInfectious<latest_infector_time+rec_day&(results$DayRemoved>rec_day|is.na(results$DayRemoved))
      infectious_names <- results$InfectedNode[infectious_index]
      infectious_ends <- pmin(results$DayRemoved[infectious_index],latest_infector_time+rec_day)
      infectious_ends[is.na(infectious_ends)] <- latest_infector_time+rec_day
      infectious_starts <- pmax(results$DayInfectious[infectious_index],rec_day)
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
      vaccinees2[iter] <- netwk[[4]]
      trial_participants2[iter] <- netwk[[5]]
      
    }
    #load(file=paste0('../AD_old/storage/results',rnd,'-',tr,'.Rdata'))
    
    # method 1: raw
    pop_sizes <- c(sum(vaccinees2),sum(trial_participants2) - sum(vaccinees2))
    zval_binary_mle[tr,1] <- calculate_zval(colSums(infectious_by_vaccine,na.rm=T)+colSums(excluded),pop_sizes)
    ve_est[tr,1]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
    # method 2: binary
    pop_sizes <- c(sum(vaccinees2),sum(trial_participants2) - sum(vaccinees2)) - colSums(excluded)
    zval_binary_mle[tr,2] <- calculate_zval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
    ve_est[tr,2]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
    # method 3: continuous
    eval_list <- get_efficacious_probabilities(results_list,vaccinees=vaccinees2,trial_participants=trial_participants2,contact_network=-1)
    zval_binary_mle[tr,3]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est[tr,3]  <- eval_list[[1]]
    # method 4: continuous + network
    eval_list <- get_efficacious_probabilities(results_list,vaccinees=vaccinees2,trial_participants=trial_participants2,contact_network=0)
    zval_binary_mle[tr,4]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est[tr,4]  <- eval_list[[1]]
    # method 5: household
    eval_list <- get_efficacious_probabilities(results_list,vaccinees2,trial_participants2,contact_network=1)
    zval_binary_mle[tr,5]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est[tr,5]  <- eval_list[[1]]
    # method 6: contact network
    eval_list <- get_efficacious_probabilities(results_list,vaccinees2,trial_participants2)
    zval_binary_mle[tr,6]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est[tr,6]  <- eval_list[[1]]
    # method 7: weight non events
    eval_list <- get_efficacious_probabilities2(netwk_list)
    zval_binary_mle[tr,7]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est[tr,7]  <- eval_list[[1]]
    # method 8: tte
    ph_results <- iterate_ph_model(netwk_list,cluster_flag=cluster_flag,pre_randomisation=F)
    zval_binary_mle[tr,8] <- ph_results[1]
    ve_est[tr,8] <- ph_results[2]
    # method 9: time-varying tte
    ph_results <- iterate_ph_model(netwk_list,cluster_flag=cluster_flag,pre_randomisation=T)
    zval_binary_mle[tr,9] <- ph_results[1]
    ve_est[tr,9] <- ph_results[2]
    
  }
  power <- apply(zval_binary_mle,2,function(x)sum(x>qnorm(1-0.05),na.rm=T)/sum(!is.na(x)))
  VE_est <- apply(ve_est,2,mean,na.rm=T)
  VE_sd <- apply(ve_est,2,sd,na.rm=T)
  return(list(power,VE_est,VE_sd))
}
for(rnd in 1:2){
  power[[rnd]] <- res[[rnd]][[1]]
  VE_est[[rnd]] <- res[[rnd]][[2]]
  VE_sd[[rnd]] <- res[[rnd]][[3]]
}

results <- data.frame(Method=c('Raw','Binary','Continuous','Continuous+cases','Household','Contact network','Non events','TTE','Time-varying TTE'),
                      Power=round(power[[2]],digits=2),T1E=round(power[[1]],digits=2))
results$ve <- paste0(round(VE_est[[2]],digits=2),' (',round(VE_sd[[2]],digits=2),')')
colnames(results) <- c('Method','Power','Type 1 error','Vaccine efficacy')
print(xtable(results))
