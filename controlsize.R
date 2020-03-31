source('set_up_script.R')

## ring vaccination trial ##################################################
nClusters <- 150
nTrials <- 1000
vaccine_efficacies <- c(0,0.8)
ref_recruit_day <- 30
registerDoParallel(cores=2)
#func <- get_efficacious_probabilities
eval_day <- 31
latest_infector_time <- eval_day - 0
power <- VE_est <- VE_sd <- list()

cluster_flag <- 0
direct_VE <- vaccine_efficacies[2]
adaptation <- ''
vaccinated_count <- infectious_count <- enrolled_count <- list()
for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- enrolled_count[[i]] <- 0
pval_binary_mle <- ve_est <- matrix(nrow=nTrials,ncol=9)
netwk_list <- list()
for(tr in 1:nTrials){
  vaccinees <- trial_participants <- recruit_times <- c()
  vaccinees2 <- trial_participants2 <- c()
  infectious_by_vaccine <- excluded <- matrix(0,nrow=nClusters,ncol=2)
  results_list <- list()
  allocation_ratio <- 0.5
  for(iter in 1:nClusters){
    ## select random person to start
    first_infected <- sample(g_name,1)
    inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
    #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
    hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
    inf_time <- min(inf_period,hosp_time)
    netwk <- simulate_contact_network(first_infected,inf_time,cluster_flag=cluster_flag,allocation_ratio=allocation_ratio,direct_VE=direct_VE)
    results_list[[iter]] <- netwk[[1]]
    results <- results_list[[iter]]
    recruit_times[iter] <- max(netwk[[3]])
    
    
    vaccinees2[iter] <- netwk[[4]]
    trial_participants2[iter] <- netwk[[5]]
    
  }
  netwk_list[[tr]] <- list(results_list,vaccinees2,trial_participants2)
  
}
 

saveRDS(netwk_list,'storage/control.Rds')
netwk_list <- readRDS('storage/control.Rds')
cont <- c()
pval <- c()
sizes <- seq(10,150,by=10)
for(i in 1:length(sizes)){   
  clus <- sizes[i]
  pval_binary_mle <- controls <- 0
  for(tr in 1:nTrials){
    netwk <- netwk_list[[tr]]
    results <- netwk[[1]][1:clus]
    vaccinees2 <- netwk[[2]][1:clus]
    trial_participants2 <- netwk[[3]][1:clus]
    eval_list <- get_efficacious_probabilities(results,vaccinees2,trial_participants2,contact_network=-1)
    pval_binary_mle[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    controls[tr] <- sum(trial_participants2-vaccinees2)
  }
  cont[i] <- mean(controls)
  pval[i] <- sum(pval_binary_mle<0.05,na.rm=T)/sum(!is.na(pval_binary_mle))
  print(c(i,mean(controls),sum(pval_binary_mle<0.05,na.rm=T)/sum(!is.na(pval_binary_mle))))
}

  

saveRDS(list(cont,pval),'storage/controlresults.Rds')
res_list <- readRDS('storage/controlresults.Rds')

cont <- res_list[[1]]
pval <- res_list[[2]]


result_table <- readRDS('storage/silo1.Rds')

controls <- result_table$`Sample size` - result_table$Vaccinated
power <- result_table$Power
labels <- result_table$Adaptation

pdf('figures/control.pdf'); par(mar=c(5,5,2,2))
cols <- rainbow(4)
plot(cont,pval,typ='l',lwd=2,cex.axis=1.5,cex.lab=1.5,frame=F,xlab='Controls',ylab='Power',xlim=range(c(cont,controls)),ylim=range(c(pval,power)))
text(x=controls,y=power,labels=labels,col=cols,cex=1.5)
for(i in 1:length(cols)){
  abline(v=controls[i],col=adjustcolor(cols[i],0.3),lwd=3)
  abline(h=power[i],col=adjustcolor(cols[i],0.3),lwd=3)
}
dev.off()