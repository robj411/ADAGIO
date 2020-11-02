source('set_up_script.R')
registerDoParallel(cores=10)

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
vaccine_efficacies <- c(0,0.7)
adaptations <- c('')
cluster_flags <- 0
ratios <- c(0.9,0.5,0.1) # ratio of seen to unseen
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,ratio=ratios,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
ref_recruit_day <- 30
func <- get_efficacious_probabilities
latest_infector_time <- eval_day - 0
base_nonrandom_scalar <- nonrandom_scalar
base_random_scalar <- random_scalar
random_edges <- length(E(random_g))
nonrandom_edges <- length(E(new_g))
total_edges <- nonrandom_edges*base_nonrandom_scalar + random_edges*base_random_scalar

## use binary weight
get_efficacious_probabilities <- get_efficacious_probabilities_bin
get_infectee_weights <- get_infectee_weights_bin



trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
  func <- get_efficacious_probabilities
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  ratio <- trial_designs$ratio[des]
  random_scalar <<- total_edges*(1-ratio)/(random_edges)
  nonrandom_scalar <<- (total_edges-random_scalar*random_edges)/nonrandom_edges
  #print(length(E(new_g))*nonrandom_scalar + length(E(random_g))*random_scalar)
  vaccinated_count <- infectious_count <- list()
  for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- 0
  exports <- enrolled_count <- c()
  zval_binary_mle3 <- ve_est3 <- zval_binary_mle2 <- ve_est2 <- zval_binary_mle <- ve_est <- c()
  for(tr in 1:nTrials){
    vaccinees2 <- trial_participants2 <- c()
    results_list <- list()
    allocation_ratio <- 0.5
    netwk_list <- list()
    weight_break <- 0
    iter <- 0
    while(weight_break<target_weight){
      iter <- iter + 1
      #for(iter in 1:nClusters){
      ## select random person to start
      first_infected <- sample(g_name,1)
      netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,end_time=20,allocation_ratio=allocation_ratio,direct_VE=direct_VE,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper)
      netwk_list[[iter]] <- netwk
      results_list[[iter]] <- netwk[[1]]
      results <- results_list[[iter]]
      
      vaccinees2[iter] <- netwk[[4]]
      trial_participants2[iter] <- netwk[[5]]
      
      probs <- get_efficacious_probabilities(results_list,vaccinees2,trial_participants2,max_time=length(results_list),contact_network=-1,observed=observed)
      weight_break <- sum(probs[[3]])
    }
    
    eval_list <- func(results_list,vaccinees=vaccinees2,trial_participants=trial_participants2,tested=F,contact_network=-1)
    zval_binary_mle2[tr]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees2)/nTrials
    enrolled_count[tr] <- sum(trial_participants2)
    infectious_count[[1]] <- infectious_count[[1]] + (observed*sum(sapply(results_list,function(x)sum(x$inTrial&!is.na(x$DayInfectious)))))/nTrials
    eval_list <- func(results_list,vaccinees2,trial_participants2,tested=T,contact_network=-1)
    zval_binary_mle3[tr]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est3[tr]  <- eval_list[[1]]
  }
  print(c(des,adaptation))
  power <- VE_est <- VE_sd <- c()
  power[1] <- sum(zval_binary_mle2>qnorm(1-0.05),na.rm=T)/sum(!is.na(zval_binary_mle2))
  ve_est2 <- ve_est2[is.finite(ve_est2)]
  VE_est[1] <- mean(ve_est2,na.rm=T)
  VE_sd[1] <- sd(ve_est2,na.rm=T)
  power[3] <- sum(zval_binary_mle3<0.05,na.rm=T)/sum(!is.na(zval_binary_mle3))
  VE_est[3] <- mean(ve_est3,na.rm=T)
  VE_sd[3] <- sd(ve_est3,na.rm=T)
  if(adaptation==''){
    power[2] <- sum(zval_binary_mle<0.05,na.rm=T)/sum(!is.na(zval_binary_mle))
    VE_est[2] <- mean(ve_est,na.rm=T)
    VE_sd[2] <- sd(ve_est,na.rm=T)
  }
  enrolled <- list(mean(enrolled_count),sd(enrolled_count))
  return(list(power, VE_est, VE_sd,vaccinated_count, infectious_count,enrolled))
}

trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolledsd <- trial_designs$enrolled <- 0

for(des in 1:nCombAdapt){
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  trial_designs$vaccinated[des] <- round(trial_results[[des]][[4]][[1]])
  trial_designs$infectious[des] <- round(trial_results[[des]][[5]][[1]])
  trial_designs$enrolled[des] <- round(trial_results[[des]][[6]][[1]])
  trial_designs$enrolledsd[des] <- round(trial_results[[des]][[6]][[2]])
  trial_designs$power[des] <- trial_results[[des]][[1]][1]
  trial_designs$VE_est[des] <- trial_results[[des]][[2]][1]
  trial_designs$VE_sd[des] <- trial_results[[des]][[3]][1]
}
subset(trial_designs,VE==0)
subset(trial_designs,VE>0)

result_table <- subset(trial_designs,VE>0)[,-1]
result_table$enrolled <- paste0(result_table$enrolled,' (',result_table$enrolledsd,')')
result_table$t1e <- subset(trial_designs,VE==0)$power
result_table$infectious2 <- subset(trial_designs,VE==0)$infectious
result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
#result_table$adapt <- as.character(result_table$adapt)
#result_table$adapt[result_table$adapt==''] <- 'None'
#result_table$cluster[result_table$cluster==0] <- 'Individual'
#result_table$cluster[result_table$cluster==1] <- 'Cluster'
result_table <- result_table[,!colnames(result_table)%in%c('cluster','weight','powertst','adapt','VE_est','VE_sd','enrolledsd','mee','prange')]

colnames(result_table) <- c('Ratio','Sample size','Symptomatic','Vaccinated','Power','Type 1 error','Symptomatic (VE=0)','VE estimate')
print(xtable(result_table,digits=c(0,1,0,0,0,2,2,0,0)), include.rownames = FALSE)

