source('set_up_script.R')

## ring vaccination trial ##################################################

nTrials <- 10
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
eval_days <- c(31,46,61,76,91)

nClusters <- 160

cls <- seq(60,160,by=10)
for(eval_day in eval_day){

#eval_day <- 31
latest_infector_time <- eval_day - 0

res_list <- list()
for(des in 1:5){
  set.seed(des)
  cluster_flag <<- trial_designs$cluster[des]
  direct_VE <<- trial_designs$VE[des]
  adaptation <<- trial_designs$adapt[des]
  res <- foreach(tr = 1:nTrials) %dopar% {
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
    return(netwk_list)
  }
  res_list[[des]] <- res
}



for(i in 1:length(cls)){
  cl <- cls[i]
  power <- vax <- ss <- rep(0,nCombAdapt)
  for(des in 1:5){
    res <- res_list[[des]]
    enrolled_count <- vaccinated_count <- 0
    result_mat <- matrix(0,nrow=nTrials,ncol=4)
    for(tr in 1:nTrials){
      netwk_list <- res[[tr]][1:cl]
      vaccinees <- sapply(netwk_list,function(netwk)netwk[[4]])
      trial_participants <- sapply(netwk_list,function(netwk)netwk[[5]])
      vaccinated_count <- vaccinated_count + sum(vaccinees)/nTrials
      enrolled_count <- enrolled_count + sum(trial_participants)/nTrials
      results_list <- lapply(netwk_list,function(x)x[[1]])
      ## regular test
      threshold <- 0.05
      eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1,observed=observed)
      pval  <- calculate_pval(eval_list[[3]],eval_list[[2]])
      ## correcting for trend 
      if(adaptation!=''){
        threshold <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                           tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio,observed=observed)
      }
      result_mat[tr,1] <- pval
      result_mat[tr,2] <- threshold
      result_mat[tr,3] <- vaccinated_count
      result_mat[tr,4] <- enrolled_count
    }
    power[des] <- sum(result_mat[,1]<result_mat[,2],na.rm=T)/sum(!is.na(result_mat[,1])&!is.na(result_mat[,2]))
    vax[des] <- sum(result_mat[,3],na.rm=T)
    ss[des] <- sum(result_mat[,4],na.rm=T)
  }
  print(c(cl,power))
  saveRDS(list(power,vax,ss),paste0('storage/cl',eval_day,cl,'.Rds'))
}

powers <- vax <- ss <- matrix(0,nrow=5,ncol=length(cls))
for(i in 1:length(cls)){
  cl <- cls[i]
  lst <- readRDS(paste0('storage/cl',eval_day,cl,'.Rds'))
  powers[,i] <- lst[[1]]
  vax[,i] <- lst[[2]]
  ss[,i] <- lst[[3]]
}
pdf(paste0('figures/sspower',eval_day,'.pdf'))

#x11(); 

par(mar=c(5,5,2,2))
ind <- which(!duplicated(powers[5,]))
plot(ss[5,ind],powers[5,ind],typ='l',lwd=2,ylim=range(powers),xlim=c(0,max(ss)),frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Sample size',ylab='Power')
cols <- rainbow(4)
for(i in 1:4) {
  ind <- which(!duplicated(powers[i,]))
  lines(ss[i,ind],powers[i,ind],col=cols[i],lwd=2)
}
legend(bty='n',x=0,y=max(powers),cex=1.25,col=c('black',cols),lwd=2,lty=1,legend=c('iRCT','Ros','Ney','TST','TS'))
adapt_days <- floor(cl/eval_day)
for(ad in 1:adapt_days)
  abline(v=max(ss)/cl*ad*eval_day,col='grey',lwd=2,lty=2)


dev.off()

}