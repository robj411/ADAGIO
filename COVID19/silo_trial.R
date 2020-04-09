source('set_up_script.R')

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
vaccine_efficacies <- c(0,0.7)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
trial_designs <- rbind(trial_designs,trial_designs[trial_designs$adapt=='',])
trial_designs$weight[(nCombAdapt+1):(nComb*(length(adaptations)+1))] <- 'binary'
trial_designs$powertst <- trial_designs$VE_esttst <- trial_designs$VE_sdtst <- trial_designs$VE_estht <- trial_designs$VE_sdht <- 
  trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
registerDoParallel(cores=12)
eval_day <- 31
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
    if(tr<6) rr_list[[tr]] <- people_per_ratio
    ## regular test
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1,observed=observed)
    pval_binary_mle2[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]
    ## correct VE test
    #eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,randomisation_ratios=randomisation_ratios,#rbht_norm=0,
    #                                           rbht_norm=ifelse(adaptation=='',1,2),
    #                                           people_per_ratio=people_per_ratio,adaptation=adaptation,contact_network=-1,observed=observed)#adaptation=adapt if rbht_norm=2
    ve_estht[tr]  <- eval_list[[1]]
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees)/nTrials
    enrolled_count[[1]] <- enrolled_count[[1]] + sum(trial_participants)/nTrials
    infectious_count[[1]] <- infectious_count[[1]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
    if(adaptation==''){
      pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
      pval_binary_mle[tr] <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      ve_est[tr]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      vaccinated_count[[2]] <- vaccinated_count[[2]] + sum(vaccinees)/nTrials
      enrolled_count[[2]] <- enrolled_count[[2]] + sum(trial_participants)/nTrials
      infectious_count[[2]] <- infectious_count[[2]] + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
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
    exports[tr] <- sum(sapply(results_list,function(x)sum(!x$inCluster)-1))
  }
  power <- VE_est <- VE_sd <- c()
  power[1] <- sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
  VE_est[1] <- mean(ve_est2,na.rm=T)
  VE_sd[1] <- sd(ve_est2,na.rm=T)
  power[3] <- sum(pval_binary_mle2<pval_binary_mle3,na.rm=T)/sum(!is.na(pval_binary_mle3)&!is.na(pval_binary_mle2))
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
  return(list(power, VE_est, VE_sd,vaccinated_count, infectious_count, enrolled_count,rr_list,mean(exports)))
}
trial_designs$mee <- 0
for(des in 1:nCombAdapt){
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  trial_designs$vaccinated[des] <- trial_results[[des]][[4]][[1]]
  trial_designs$infectious[des] <- trial_results[[des]][[5]][[1]]
  trial_designs$enrolled[des] <- trial_results[[des]][[6]][[1]]
  #if(adaptation==''){
  #  trial_designs$vaccinated[des+nComb] <- trial_results[[des]][[4]][[2]]
  #  trial_designs$infectious[des+nComb] <- trial_results[[des]][[5]][[2]]
  #  trial_designs$enrolled[des+nComb] <- trial_results[[des]][[6]][[2]]
  #}
  trial_designs$power[des] <- trial_results[[des]][[1]][1]
  trial_designs$VE_est[des] <- trial_results[[des]][[2]][1]
  trial_designs$VE_sd[des] <- trial_results[[des]][[3]][1]
  trial_designs$powertst[des] <- trial_results[[des]][[1]][3]
  trial_designs$VE_esttst[des] <- trial_results[[des]][[2]][3]
  trial_designs$VE_sdtst[des] <- trial_results[[des]][[3]][3]
  trial_designs$VE_estht[des] <- trial_results[[des]][[2]][4]
  trial_designs$VE_sdht[des] <- trial_results[[des]][[3]][4]
  trial_designs$mee[des] <- trial_results[[des]][[8]]
  #if(adaptation==''){
  #  trial_designs$power[des+nComb] <- trial_results[[des]][[1]][2]
  #  trial_designs$VE_est[des+nComb] <- trial_results[[des]][[2]][2]
  #  trial_designs$VE_sd[des+nComb] <- trial_results[[des]][[3]][2]
  #  trial_designs$powertst[des+nComb] <- trial_results[[des]][[1]][3]
  #  trial_designs$VE_esttst[des+nComb] <- trial_results[[des]][[2]][3]
  #  trial_designs$VE_sdtst[des+nComb] <- trial_results[[des]][[3]][3]
  #  trial_designs$VE_estht[des+nComb] <- trial_results[[des]][[2]][4]
  #  trial_designs$VE_sdht[des+nComb] <- trial_results[[des]][[3]][4]
  #  trial_designs$mee[des+nComb] <- trial_results[[des]][[8]]
  #}
}
subset(trial_designs,VE==0)
subset(trial_designs,VE>0)
saveRDS(trial_designs,'storage/silo_trials.Rds')
result_table <- subset(trial_designs,VE>0)[,c(3:15)]
result_table$t1e <- subset(trial_designs,VE==0)$power
result_table$t1etst <- subset(trial_designs,VE==0)$powertst
result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
result_table <- result_table[,!colnames(result_table)%in%c('VE_est','VE_sd')]
result_table$htVE <- paste0(round(result_table$VE_estht,2),' (',round(result_table$VE_sdht,2),')')
result_table <- result_table[,!colnames(result_table)%in%c('VE_estht','VE_sdht')]
#result_table$tstVE <- paste0(round(result_table$VE_esttst,2),' (',round(result_table$VE_sdtst,2),')')
result_table <- result_table[,!colnames(result_table)%in%c('VE_esttst','VE_sdtst')]
result_table$adapt <- as.character(result_table$adapt)
result_table$adapt[result_table$adapt==''] <- 'None'
result_table$nmee <- subset(trial_designs,VE==0)$mee - subset(trial_designs,VE>0)$mee
#result_table$cluster[result_table$cluster==0] <- 'Individual'
#result_table$cluster[result_table$cluster==1] <- 'Cluster'
colnames(result_table) <- c('Adaptation','Weighting','Sample size','Infectious','Vaccinated','Power','Power (corrected)',
                            'Type 1 error','Type 1 error (corrected)','VE estimate','VE estimate (TH)','NMEE')
print(xtable(result_table), include.rownames = FALSE)


change_days <- trial_results[[1]][[7]][[1]][,2]
adaptation_days <- c(1:change_days[1],c(sapply(2:length(change_days),function(x)(1+change_days[x-1]):change_days[x])))
get_allocation_vector <- function(x){
  sapply(1:5,function(y){
    vals <- trial_results[[x]][[7]][[y]][,3]
    alloc <- c(rep(0.5,change_days[1]),c(sapply(2:length(change_days),function(x)rep(vals[x-1],(change_days[x]-change_days[x-1])))))
    alloc
  })
}
cols <- rainbow(4)
{pdf('allocation_probability.pdf',width=10,height=5);
  #x11(width=10,height=5);
  par(mfrow=c(1,2),mar=c(5,5,2,2))
  matplot(adaptation_days,get_allocation_vector(1),typ='l',col=adjustcolor(cols[ceiling(1/2)],0.5),frame=F,lty=1,lwd=2,xlab='Day',ylab='Allocation probability (VE=0)',cex.axis=1.5,cex.lab=1.5,ylim=0:1)
  for(j in seq(3,7,by=2)){
    matplot(adaptation_days,get_allocation_vector(j),typ='l',col=adjustcolor(cols[ceiling(j/2)],0.5),lty=1,lwd=2,add=T)
  }    
  legend(x=-0,y=1.05,legend=c('Ney','Ros','TST','TS'),col=cols,lwd=2,bty='n')
  matplot(adaptation_days,get_allocation_vector(2),typ='l',col=adjustcolor(cols[ceiling(2/2)],0.5),frame=F,lty=1,lwd=2,xlab='Day',ylab='Allocation probability (VE=0.7)',cex.axis=1.5,cex.lab=1.5,ylim=0:1)
  for(j in seq(4,8,by=2)){
    matplot(adaptation_days,get_allocation_vector(j),typ='l',col=adjustcolor(cols[ceiling(j/2)],0.5),lty=1,lwd=2,add=T)
  }    
  dev.off()
}
