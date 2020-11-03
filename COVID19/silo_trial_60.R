source('set_up_script.R')
registerDoParallel(cores=5)

## ring vaccination trial ##################################################
nClusters <- 60
nTrials <- 10000
vaccine_efficacies <- c(0.7)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
#trial_designs <- rbind(trial_designs,trial_designs[trial_designs$adapt=='',])
#trial_designs$weight[(nCombAdapt+1):(nComb*(length(adaptations)+1))] <- 'binary'
ref_recruit_day <- 30
latest_infector_time <- eval_day - 0

trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
  set.seed(des)
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  vaccinated_count <- infectious_count <- rr_list <- list()
  for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- 0
  pval_binary_mle3 <- zval_binary_mle3 <- ve_est3 <- pval_binary_mle2 <- zval_binary_mle2 <- ve_est2 <- pval_binary_mle <- ve_est <- ve_estht <- c()
  exports <- enrolled_count <- offline_allocation_ratios <- c()
  for(tr in 1:nTrials){
    randomisation_ratios <- c()
    people_per_ratio <- c()
    vaccinees <- trial_participants <- c()
    infectious_by_vaccine <- excluded <- c()
    results_list <- list()
    allocation_ratio <- offline_allocation_ratio <- 0.5
    netwk_list <- list()
    weight_break <- 0
    iter <- 0
    while(iter<nClusters&offline_allocation_ratio<0.99){
      iter <- iter + 1
      set.seed(iter*nTrials+tr)
    #for(iter in 1:nClusters){
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
        allocation_ratios <- response_adapt(fails,pop_sizes2,days=iter,adaptation)
        allocation_ratio <- allocation_ratios[1]
        ## uncomment to stop TS methods early
        offline_allocation_ratio <- allocation_ratios[2]
        people_per_ratio <- rbind(people_per_ratio,c(sum(trial_participants),iter,allocation_ratio))
        #if(allocation_ratio==0) break
        weight_break <- sum(probs[[3]])
      }else if(iter >= eval_day && sum(vaccinees)>0){
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
        weight_break <- sum(probs[[3]])
      }
    }
    rr_list[[tr]] <- people_per_ratio
    ## regular test
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1,observed=observed)
    zval_binary_mle2[tr]  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    pval_binary_mle2[tr]  <- dnorm(zval_binary_mle2[tr])
    ve_est2[tr]  <- eval_list[[1]]
    offline_allocation_ratios[tr] <- offline_allocation_ratio
    ## correct VE test
    #eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,randomisation_ratios=randomisation_ratios,#rbht_norm=0,
    #                                           rbht_norm=ifelse(adaptation=='',1,2),
    #                                           people_per_ratio=people_per_ratio,adaptation=adaptation,contact_network=-1,observed=observed)#adaptation=adapt if rbht_norm=2
    ve_estht[tr]  <- eval_list[[1]]
    vaccinated_count[[1]] <- vaccinated_count[[1]] + sum(vaccinees)/nTrials
    enrolled_count[tr] <- sum(trial_participants)
    infectious_count[[1]] <- infectious_count[[1]] + (observed*sum(sapply(results_list,function(x)sum(x$inTrial&!is.na(x$DayInfectious)))))/nTrials
    if(adaptation==''){
      pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
      pval_binary_mle[tr] <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      ve_est[tr]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      vaccinated_count[[2]] <- vaccinated_count[[2]] + sum(vaccinees)/nTrials
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
      zval_binary_mle3[tr] <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                                 tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio,observed=observed)
    
    ## exports
    exports[tr] <- sum(sapply(results_list,function(x)sum(!x$inCluster)-1))/length(results_list)*100
  }
  power <- VE_est <- VE_sd <- c()
  power[1] <- sum(zval_binary_mle2>qnorm(0.95)|offline_allocation_ratios>0.99,na.rm=T)/sum(!is.na(zval_binary_mle2))
  VE_est[1] <- mean(ve_est2,na.rm=T)
  VE_sd[1] <- sd(ve_est2,na.rm=T)
  power[3] <- sum(zval_binary_mle2>zval_binary_mle3|offline_allocation_ratios>0.99,na.rm=T)/sum(!is.na(zval_binary_mle3)&!is.na(zval_binary_mle2))
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
  power[2] <- infotheo::entropy(discretize(pval_binary_mle2,disc='equalwidth')) #quantile(pval_binary_mle2,0.95) - quantile(pval_binary_mle2,0.05)
  saveRDS(list(zval_binary_mle2,zval_binary_mle3),paste0('storage/silop',des,'.Rds'))
  enrolled <- list(mean(enrolled_count),sd(enrolled_count))
  print(c(des,power))
  return(list(power, VE_est, VE_sd,vaccinated_count, infectious_count, enrolled,rr_list,mean(exports)))
}
saveRDS(trial_results,'storage/silo_60.Rds')
trial_designs$mee <- trial_designs$powertst <- trial_designs$power <- trial_designs$vaccinated <- trial_designs$infectious <- 
  trial_designs$daysd <- trial_designs$day <- trial_designs$enrolledsd <- trial_designs$enrolled <- 0
for(des in 1:nCombAdapt){
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  trial_designs$vaccinated[des] <- round(trial_results[[des]][[4]][[1]])
  trial_designs$infectious[des] <- round(trial_results[[des]][[5]][[1]])
  trial_designs$enrolled[des] <- round(trial_results[[des]][[6]][[1]])
  trial_designs$enrolledsd[des] <- round(trial_results[[des]][[6]][[2]])
  trial_designs$day[des] <- round(trial_results[[des]][[6]][[1]]/enrolled_per_contact)+eval_day
  trial_designs$daysd[des] <- round(trial_results[[des]][[6]][[2]]/enrolled_per_contact)
  trial_designs$power[des] <- trial_results[[des]][[1]][1]
  trial_designs$powertst[des] <- trial_results[[des]][[1]][3]
  if(adaptation=='')
    trial_designs$powertst[des] <- trial_results[[des]][[1]][1]
  trial_designs$mee[des] <- trial_results[[des]][[8]]
}
subset(trial_designs,VE>0)
result_table <- subset(trial_designs,VE>0)[,-c(1:2)]
result_table$enrolled <- paste0(result_table$enrolled,' (',result_table$enrolledsd,')')
result_table$vax2 <- round((100-result_table$day)*enrolled_per_contact*trial_designs$powertst)
result_table$mee <- round((result_table$day-eval_day)*result_table$mee/100)
result_table$vax3 <- result_table$vax2 + result_table$vaccinated
result_table$vax4 <- round((100-result_table$day)*320*trial_designs$powertst) + result_table$vaccinated
result_table$day <- paste0(result_table$day,' (',result_table$daysd,')')
result_table$adapt <- as.character(result_table$adapt)
result_table$adapt[result_table$adapt==''] <- 'FR'
result_table$adapt[result_table$adapt=='Ros'] <- 'Ros. et al.'
result_table$adapt[result_table$adapt=='Ney'] <- 'Ney.'
#result_table$nmee <- subset(trial_designs,VE==0)$mee - subset(trial_designs,VE>0)$mee
result_table <- result_table[,!colnames(result_table)%in%c('vax2','daysd','weight','VE_est','VE_sd','enrolledsd','prange')]
colnames(result_table) <- c('Adaptation','Sample size','Duration','Number of confirmed cases','Vaccinated in trial','Power','Power (adjusted)','Export infections in trial','assuming 32 per day','assuming 320 per day')
print(xtable(result_table,digits=c(0,0,0,0,0,0,2,2,0,0,0)), include.rownames = FALSE)


# first index is the trial; 5:8 is TS designs
# second index is 7, extracting the rr_list
# third index is y over the number of samples
# fourth index ,3 extracts the ratio column
## propensities to stop for each 
stop_indices <- sapply(5:8,function(x)
      sapply(1:length(trial_results[[x]][[7]]),function(y){
        if(any(trial_results[[x]][[7]][[y]][,3]>0.99))
          which(trial_results[[x]][[7]][[y]][,3]>0.99)[1]
        else 0
    })
)
print(t(sapply(1:max(stop_indices),function(x)
  apply(stop_indices,2,function(y)
    sum(y<=x&y>0)/length(y)
    )
  )
))

## overall propensity to stop
print(sapply(5:8,function(x)
  sum(sapply(1:length(trial_results[[x]][[7]]),function(y){
    any(trial_results[[x]][[7]][[y]][,3]>0.99)
  }))
  /length(trial_results[[x]][[7]])
))


## sample size for each 
print(sapply(5:8,function(x){
  sss <- sapply(1:length(trial_results[[x]][[7]]),function(y){
    rr <- trial_results[[x]][[7]][[y]]
    rs <- rr[,3]
    ss <- rr[,1]
    if(any(rs>0.99)){
      ind <- which(rs>0.99)[1]
      ss[ind]
    }else{
      trial_results[[x]][[6]][[1]]
    }
  })
  c(mean(sss),sd(sss))
}
))


n_change_days <- sapply(trial_results[[1]][[7]],nrow)
change_days <- trial_results[[1]][[7]][[which.max(n_change_days)]][,2]
adaptation_days <- c(1:change_days[1],
                     c(sapply(2:length(change_days),
                              function(x)(1+change_days[x-1]):change_days[x])))
get_allocation_vector <- function(x){
  sapply(1:5,function(y){
    vals <- trial_results[[x]][[7]][[y]][,3]
    alloc <- c(rep(0.5,change_days[1]),c(sapply(2:length(change_days),function(x)rep(vals[x-1],(change_days[x]-change_days[x-1])))))
    alloc
  })
}
cols <- rainbow(4)
{pdf('figures/allocation_probability.pdf',width=10,height=5);
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


