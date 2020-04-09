source('set_up_script.R')

## ring vaccination trial ##################################################

nTrials <- 1000
vaccine_efficacies <- c(0.7)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
trial_designs$powertst <- trial_designs$VE_esttst <- trial_designs$VE_sdtst <- trial_designs$VE_estht <- trial_designs$VE_sdht <- 
  trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
registerDoParallel(cores=32)
eval_day <- 31
latest_infector_time <- eval_day - 0
nClusters <- 50


power <- vax <- ss <- infected <- rep(0,nCombAdapt)
while(any(power<0.8)){
  nClusters <- nClusters + 10
  toupdate <- which(power<0.8)
  for(des in toupdate){
    set.seed(des)
    cluster_flag <<- trial_designs$cluster[des]
    direct_VE <<- trial_designs$VE[des]
    adaptation <<- trial_designs$adapt[des]
    res <- foreach(tr = 1:nTrials,.combine=rbind) %dopar% {
      vaccinated_count <- infectious_count <- enrolled_count <- 0
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
      vaccinated_count <- vaccinated_count + sum(vaccinees)/nTrials
      enrolled_count <- enrolled_count + sum(trial_participants)/nTrials
      infectious_count <- infectious_count + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
      ## regular test
      threshold <- 0.05
      eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1,observed=observed)
      pval  <- calculate_pval(eval_list[[3]],eval_list[[2]])
      ## correcting for trend 
      if(adaptation!=''){
        threshold <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                           tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio,observed=observed)
      }
      return(c(pval,threshold,vaccinated_count,enrolled_count,infectious_count))
    }
    power[des] <- sum(res[,1]<res[,2],na.rm=T)/sum(!is.na(res[,1])&!is.na(res[,2]))
    vax[des] <- mean(res[,3],na.rm=T)
    ss[des] <- mean(res[,4],na.rm=T)
    infected[des] <- mean(res[,5],na.rm=T)
  }
  print(c(nClusters,power))
  saveRDS(list(power,vax,ss,infected),paste0('storage/cl',nClusters,'.Rds'))
}
