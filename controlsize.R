source('set_up_script.R')
get_infectee_weights_original <- get_infectee_weights
get_infectee_weights_binary <- function(results,ve_point_est,contact_network=2,tested=F){
  resrec <- results$RecruitmentDay
  nonna <- results[!is.na(resrec) & results$DayInfectious>resrec,]
  nonnavac <- nonna$vaccinated==T
  nonnainf <- nonna$DayInfectious
  nonnarec <- nonna$RecruitmentDay+10
  incl <- nonnainf >= nonnarec
  weight_hh_rem <- cbind(nonnavac & incl,
                         !nonnavac & incl)
  infectee_names <- nonna$InfectedNode
  return(list(weight_hh_rem,infectee_names))
}

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
# for(i in 1:2) vaccinated_count[[i]] <- infectious_count[[i]] <- enrolled_count[[i]] <- 0
# pval_binary_mle <- ve_est <- matrix(nrow=nTrials,ncol=9)
# netwk_list <- list()
# for(tr in 1:nTrials){
#   vaccinees <- trial_participants <- recruit_times <- c()
#   vaccinees2 <- trial_participants2 <- c()
#   infectious_by_vaccine <- excluded <- matrix(0,nrow=nClusters,ncol=2)
#   results_list <- list()
#   allocation_ratio <- 0.5
#   for(iter in 1:nClusters){
#     ## select random person to start
#     first_infected <- sample(g_name,1)
#     inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
#     #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
#     hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
#     netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,allocation_ratio=allocation_ratio,direct_VE=direct_VE)
#     results_list[[iter]] <- netwk[[1]]
#     results <- results_list[[iter]]
#     recruit_times[iter] <- max(netwk[[3]])
#     
#     
#     vaccinees2[iter] <- netwk[[4]]
#     trial_participants2[iter] <- netwk[[5]]
#     
#   }
#   netwk_list[[tr]] <- list(results_list,vaccinees2,trial_participants2)
#   
# }
#  
# 
# saveRDS(netwk_list,'storage/control.Rds')
netwk_list <- readRDS('storage/control.Rds')
cont <- contb <- c()
pval <- pvalb <- c()
sizes <- seq(20,150,by=10)
for(i in 1:length(sizes)){   
  clus <- sizes[i]
  pval_binary_mle <- pval_binary_mle2 <- controls <- 0
  #infectious_by_vaccine <- excluded <- matrix(0,nrow=nTrials,ncol=2)
  for(tr in 1:nTrials){
    netwk <- netwk_list[[tr]]
    results <- netwk[[1]][1:clus]
    vaccinees2 <- netwk[[2]][1:clus]
    trial_participants2 <- netwk[[3]][1:clus]
    ## 2 binary
    tab <- do.call(rbind,results)
    infectious_by_vaccine <- c(sum(tab$vaccinated&tab$DayInfectious>tab$RecruitmentDay+9,na.rm=T),
                               sum(!tab$vaccinated&tab$inTrial&tab$DayInfectious>tab$RecruitmentDay+9,na.rm=T))
    excluded <- c(sum(tab$vaccinated&tab$DayInfectious<tab$RecruitmentDay+10,na.rm=T),
                  sum(!tab$vaccinated&tab$inTrial&tab$DayInfectious<tab$RecruitmentDay+10,na.rm=T))
    pop_sizes <- c(sum(vaccinees2),sum(trial_participants2) - sum(vaccinees2)) - excluded
    pval_binary_mle2[tr]  <- calculate_pval(infectious_by_vaccine,pop_sizes)
    ## 3 continuous
    #eval_list <- get_efficacious_probabilities(results,vaccinees2,trial_participants2,contact_network=-1)
    #pval_binary_mle[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    #controls[tr] <- sum(trial_participants2-vaccinees2)
  }
  #cont[i] <- mean(controls)
  #pval[i] <- sum(pval_binary_mle<0.05,na.rm=T)/sum(!is.na(pval_binary_mle))
  pvalb[i] <- sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
}

  

#saveRDS(list(cont,pval),'storage/controlresults.Rds')
saveRDS(pvalb,'storage/controlbinresults.Rds')
res_list <- readRDS('storage/controlresults.Rds')
pvalb <- readRDS('storage/controlbinresults.Rds')
print(pvalb)
cont <- res_list[[1]][-1]
pval <- res_list[[2]][-1]


result_table <- readRDS('storage/silo1.Rds')
tokeep <- subset(result_table,Adaptation%in%c('Ney','Ros','TST','TS')&Randomisation=='Individual')
result_b <- readRDS('storage/binsilo.Rds')
tokeepb <- subset(result_b,Adaptation%in%c('Ney','Ros','TST','TS')&Randomisation=='Individual')

controls <- tokeepb$`Sample size` - tokeepb$Vaccinated
power <- tokeepb$Power
labels <- tokeepb$Adaptation

pdf('figures/controlbin.pdf'); par(mar=c(5,5,2,2))
cols <- rainbow(4)
plot(cont,pvalb,typ='l',lwd=2,cex.axis=1.5,cex.lab=1.5,frame=F,xlab='Controls',ylab='Power',xlim=range(c(cont,controls)),ylim=range(c(pval,pvalb,power)))
lines(cont,pval,typ='l',lwd=2,col='grey')
text(x=controls,y=power,labels=labels,col=cols,cex=1.5)
for(i in 1:length(cols)){
  abline(v=controls[i],col=adjustcolor(cols[i],0.3),lwd=3)
  abline(h=power[i],col=adjustcolor(cols[i],0.3),lwd=3)
}
dev.off()