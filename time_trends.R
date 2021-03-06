source('set_up_script.R')
## cores
registerDoParallel(cores=32)
## create storage/timetrend*

## can we infer a trend? ##################################################
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


direct_VE <- 0.0
reps <- 1000
nIter <- 100
adaptation <- 'TST'
pval_binary_mle2 <- pval_binary_mle21 <- ve_est2 <- ve_est21 <- pval_threshold <- c()
eval_day <<- 31
latest_infector_time <- eval_day - 0
func <- get_efficacious_probabilities
rates <- -seq(5e-7,5e-6,by=1e-6)
t1e <- t1e1 <- c()
nClusters <- nIter
pval_binary_mle <- pval_binary_mle1 <- matrix(0,nrow=reps,ncol=length(rates))
#t1elist <- list()
t1elist <- foreach(i = rep(1:length(rates),2),j=rep(1:2,each=length(rates))) %do% { #for(i in 1:length(rates)){
  #t1elist[[5*(j-1)+i]] <- matrix(0,nrow=reps,ncol=3)
  direct_VE <- c(0,0.7)[j]
  per_time_step <- rates[i]
  base_rate <- - 130 * rates[i]
  if(file.exists(paste0('storage/timetrend',i,j,'.Rds'))){
    all_reps <- readRDS(paste0('storage/timetrend',i,j,'.Rds'))
    pval_binary_mle2 <- all_reps[,1]
    pval_binary_mle21 <- all_reps[,2]
    pval_threshold <- all_reps[,3]
    return(c(sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2)), 
             sum(pval_binary_mle21<0.05,na.rm=T)/sum(!is.na(pval_binary_mle21)),
             sum(pval_binary_mle2<pval_threshold,na.rm=T)/sum(!is.na(pval_binary_mle2))))
  }
  all_reps <- foreach(rp = 1:reps,.combine=rbind) %dopar% {
    #profvis({
    allocation_ratio <- 0.5
    results_list <- netwk_list <- list()
    #vaccinees <- trial_participants <- 
    people_per_ratio <- c()
    vaccinees2 <- trial_participants2 <- randomisation_ratios <- c()
    infectious_by_vaccine <- excluded <- matrix(0,nrow=nIter,ncol=2)
    for(iter in 1:nIter){
      randomisation_ratios[iter] <- allocation_ratio
      ## select random person to start
      first_infected <- sample(g_name,1)
      netwk <- simulate_contact_network(first_infected,end_time=eval_day,start_day=iter,from_source=per_time_step,
                                        cluster_flag=0,allocation_ratio=allocation_ratio,direct_VE=direct_VE,base_rate=base_rate)
      netwk_list[[iter]] <- netwk
      results_list[[iter]] <- netwk[[1]]
      cluster_size[iter] <- netwk[[2]]
      recruit_times[iter] <- max(netwk[[3]])
      results <- results_list[[iter]]
      vax <- results$vaccinated
      too_early <- results$DayInfectious<results$RecruitmentDay+10
      infectious_by_vaccine[iter,] <- c(sum(vax&!too_early),sum(!vax&results$inTrial&!too_early))
      excluded[iter,] <- c(sum(vax&too_early),sum(!vax&results$inTrial&too_early))
      
      vaccinees2[iter] <- netwk[[4]]
      trial_participants2[iter] <- netwk[[5]]
      if(adaptation!=''&&iter %% eval_day == 0){
        get_infectee_weights <- get_infectee_weights_original ##binary##
        #probs <- get_efficacious_probabilities2(netwk_list,max_time=length(results_list))#func(results_list,vaccinees2,trial_participants2,max_time=length(results_list),contact_network=-1)
        probs <- get_efficacious_probabilities(results_list,vaccinees2,trial_participants2,max_time=length(results_list),contact_network=-1)
        pop_sizes2 <- probs[[2]]
        fails <- probs[[3]]
        allocation_ratio <- response_adapt(fails,pop_sizes2,days=iter,adaptation=adaptation)
        people_per_ratio <- rbind(people_per_ratio,c(sum(trial_participants2),iter,allocation_ratio))
        #0.9^(iter/nIter)/(0.9^(iter/nIter)+0.1^(iter/nIter))#
      }
    }
    #print(allocation_ratio)
    #})
    # method 3: continuous
    #eval_list <- get_efficacious_probabilities(results_list,vaccinees2,trial_participants2,contact_network=-1)
    #pval_binary_mle2[rp]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    #ve_est2[rp]  <- eval_list[[1]]
    # method 2: binary
    pop_sizes <- c(sum(vaccinees2),sum(trial_participants2) - sum(vaccinees2)) - colSums(excluded)
    pval_binary_mle2[rp]  <- calculate_pval(fails=colSums(infectious_by_vaccine,na.rm=T),sizes=pop_sizes)
    ve_est2[rp] <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
    get_infectee_weights <- get_infectee_weights_binary
    pval_threshold[rp] <- trend_robust_function(results_list,vaccinees=vaccinees2,trial_participants=trial_participants2,contact_network=-1,
                                                 tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio)
    # method 7: weight non events
    get_infectee_weights <- get_infectee_weights_original
    eval_list <- get_efficacious_probabilities2(netwk_list)
    pval_binary_mle21[rp]  <- calculate_pval(fails=eval_list[[3]],sizes=eval_list[[2]])
    ve_est21[rp]  <- eval_list[[1]]
    #print(c(pval_binary_mle2,ve_est2,allocation_ratio))
    return(c(pval_binary_mle2[rp],pval_binary_mle21[rp],pval_threshold[rp]))
  }
  saveRDS(all_reps,paste0('storage/timetrend',i,j,'.Rds'))
  pval_binary_mle2 <- all_reps[,1]
  pval_binary_mle21 <- all_reps[,2]
  pval_threshold <- all_reps[,3]
  
  #t1e[i] <- sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
  #pval_binary_mle[,i] <- pval_binary_mle2
  #t1e1[i] <- sum(pval_binary_mle21<0.05,na.rm=T)/sum(!is.na(pval_binary_mle21))
  #pval_binary_mle1[,i] <- pval_binary_mle21
  #print(c(i,t1e[i],mean(ve_est2),sd(ve_est2)))
  #return(pval_threshold)
  return(c(sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2)), 
           sum(pval_binary_mle21<0.05,na.rm=T)/sum(!is.na(pval_binary_mle21)),
           sum(pval_binary_mle2<pval_threshold,na.rm=T)/sum(!is.na(pval_binary_mle2))))
  #hist(rpois(1000,mean(counts-1))+1)
  #hist(counts)
}
print(t1elist)
saveRDS(t1elist,'storage/t1es.Rds')
t1elist <- readRDS('storage/t1es.Rds')
t1e <- sapply(t1elist,function(x)x[1])[1:length(rates)]
t1e1 <- sapply(t1elist,function(x)x[2])[1:length(rates)]
t1e3 <- sapply(t1elist,function(x)x[3])[1:length(rates)]
power <- sapply(t1elist,function(x)x[1])[1:length(rates)+length(rates)]
power2 <- sapply(t1elist,function(x)x[2])[1:length(rates)+length(rates)]
power3 <- sapply(t1elist,function(x)x[3])[1:length(rates)+length(rates)]
cols <- c('darkorange2','navyblue','hotpink','grey','turquoise')
#pdf('trendt1e2.pdf',height=5,width=10); par(mar=c(5,5,2,2),mfrow=c(1,2))
#matplot(sapply(rates,function(x)1:130*x - 130*x),typ='l',col=cols,lwd=3,lty=1,xlab='Day',ylab='Background rate',cex.lab=1.5,cex.axis=1.5,frame=F)
#plot(-rates,t1e,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Type 1 error',ylim=c(0.03,0.15),xaxt='n')
#axis(1,-rates,-rates,cex.axis=1.5)
#dev.off()
#pdf('trendt1e2adj.pdf',height=5,width=10); par(mar=c(5,5,2,2),mfrow=c(1,2))
#matplot(sapply(rates,function(x)1:130*x - 130*x),typ='l',col=cols,lwd=3,lty=1,xlab='Day',ylab='Background rate',cex.lab=1.5,cex.axis=1.5,frame=F)
#plot(-rates,t1e3,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Type 1 error',ylim=c(0.03,0.15),xaxt='n')
#axis(1,-rates,-rates,cex.axis=1.5)
#dev.off()
#pdf('trendt1e6.pdf',height=5,width=10); par(mar=c(5,5,2,2),mfrow=c(1,2))
#matplot(sapply(rates,function(x)1:130*x - 130*x),typ='l',col=cols,lwd=3,lty=1,xlab='Day',ylab='Background rate',cex.lab=1.5,cex.axis=1.5,frame=F)
#plot(-rates,t1e1,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Type 1 error',ylim=c(0.03,0.15),xaxt='n')
#axis(1,-rates,-rates,cex.axis=1.5)
#dev.off()

#pdf('figures/trendt1e.pdf',height=5,width=10); par(mar=c(5,5,2,2),mfrow=c(1,2))
#matplot(sapply(rates,function(x)1:130*x - 130*x),typ='l',col=cols,lwd=3,lty=1,xlab='Day',ylab='Background rate',cex.lab=1.5,cex.axis=1.5,frame=F)
#plot(-rates,t1e,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Type 1 error',ylim=c(0.0,0.15),xaxt='n')
#axis(1,-rates,-rates,cex.axis=1.5)
#points(-rates,t1e1,typ='p',cex=2,pch=17,col=cols)
#points(-rates,t1e3,typ='p',cex=2,pch=15,col=cols)
#legend(x=-rates[1],0.16,bty='n',legend=c('Method 2','Method 6','Method 2 corrected'),col=cols[1],pch=c(19,17,15),cex=1.5)
#dev.off()


pdf('figures/trend.pdf',height=5,width=15); 
#x11(height=5,width=15); 
par(mar=c(5,5,2,2),mfrow=c(1,3))
matplot(sapply(rates,function(x)1:130*x - 130*x),typ='l',col=cols,lwd=3,lty=1,xlab='Day',ylab='Background rate',cex.lab=1.5,cex.axis=1.5,frame=F)
plot(-rates,t1e,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Type 1 error',ylim=c(0.0,0.15),xaxt='n')
axis(1,-rates,-rates,cex.axis=1.5)
points(-rates,t1e1,typ='p',cex=2,pch=17,col=cols)
points(-rates,t1e3,typ='p',cex=2,pch=15,col=cols)
legend(x=-rates[1],0.16,bty='n',legend=c('Method 2','Method 7','Method 2 corrected'),col=cols[1],pch=c(19,17,15),cex=1.5)
plot(-rates,power,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Power',ylim=c(0.29,1),xaxt='n')
points(-rates,power2,typ='p',cex=2,pch=17,col=cols)
points(-rates,power3,typ='p',cex=2,pch=15,col=cols)
axis(1,-rates,-rates,cex.axis=1.5)
#legend(x=-rates[1],0.75,bty='n',legend=c('Method 2','Method 6','Method 2 corrected'),col=cols[1],pch=c(19,17,15))
dev.off()
