source('set_up_script.R')
registerDoParallel(cores=32)
## saves to 'storage/timetrend',i,j,'.Rds'
## can we infer a trend? ##################################################

direct_VE <- 0.0
reps <- 1000
nIter <- 100
adaptation <- 'TST'
zval_binary_mle2 <- ve_est2 <- zval_threshold <- c()
latest_infector_time <- eval_day - 0
func <- get_efficacious_probabilities
rates <- -seq(5e-7,5e-6,by=1e-6)
t1e <- t1e1 <- c()
nClusters <- nIter
pval_binary_mle <- pval_binary_mle1 <- matrix(0,nrow=reps,ncol=length(rates))
t1elist <- foreach(i = rep(1:length(rates),2),j=rep(1:2,each=length(rates))) %do% { #for(i in 1:length(rates)){
  direct_VE <- c(0,0.7)[j]
  per_time_step <- rates[i]
  base_rate <- - 130 * rates[i]
  if(file.exists(paste0('storage/timetrend',i,j,'.Rds'))){
    all_reps <- readRDS(paste0('storage/timetrend',i,j,'.Rds'))
    zval_binary_mle2 <- all_reps[,1]
    zval_threshold <- all_reps[,2]
    return(c(sum(dnorm(zval_binary_mle2)<0.05,na.rm=T)/sum(!is.na(zval_binary_mle2)), 
             sum(zval_binary_mle2>zval_threshold,na.rm=T)/sum(!is.na(zval_binary_mle2))))
  }
  all_reps <- foreach(rep = 1:reps,.combine=rbind) %dopar% {
    #profvis({
    allocation_ratio <- 0.5
    results_list <- list()
    vaccinees <- trial_participants <- people_per_ratio <- randomisation_ratios <- c()
    
    weight_break <- 0
    iter <- 0
    while(weight_break<target_weight){
      iter <- iter + 1
      #for(iter in 1:nIter){
      randomisation_ratios[iter] <- allocation_ratio
      ## select random person to start
      first_infected <- sample(g_name[eligible_first_person],1)
      netwk <- simulate_contact_network(first_infected,end_time=eval_day,start_day=iter,from_source=per_time_step,individual_recruitment_times=T,
                                        cluster_flag=0,allocation_ratio=allocation_ratio,direct_VE=direct_VE,base_rate=base_rate,spread_wrapper=covid_spread_wrapper)
      
      results_list[[iter]] <- netwk[[1]]
      cluster_size[iter] <- netwk[[2]]
      recruit_times[iter] <- ifelse(length(netwk[[3]])>0,max(netwk[[3]]),0)
      results <- results_list[[iter]]
      vax <- results$vaccinated
      too_early <- results$DayInfectious<results$RecruitmentDay+10
      
      ##!! weighting non-events
      rec_day <- recruit_times[iter]
      infectious_index <- results$DayInfectious<latest_infector_time+rec_day&(results$DayRemoved>rec_day|is.na(results$DayRemoved))
      infectious_names <- results$InfectedNode[infectious_index]
      infectious_ends <- pmin(results$DayRemoved[infectious_index],latest_infector_time+rec_day)
      infectious_ends[is.na(infectious_ends)] <- latest_infector_time+rec_day
      infectious_starts <- pmax(results$DayInfectious[infectious_index],rec_day)
      vaccinees[iter] <- netwk[[4]]
      trial_participants[iter] <- netwk[[5]]
      if(adaptation!=''&&iter %% eval_day == 0){
        probs <- func(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1)
        pop_sizes2 <- probs[[2]]
        fails <- probs[[3]]
        allocation_ratio <- response_adapt(fails,pop_sizes2,days=iter,adaptation=adaptation)
        people_per_ratio <- rbind(people_per_ratio,c(sum(trial_participants),iter,allocation_ratio))
        #0.9^(iter/nIter)/(0.9^(iter/nIter)+0.1^(iter/nIter))#
        weight_break <- sum(probs[[3]])
      }else if(iter >= eval_day && sum(vaccinees)>0){
        probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
        weight_break <- sum(probs[[3]])
      }
    }
    #print(allocation_ratio)
    #})
    # 
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1)
    zval_binary_mle2  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    ve_est2  <- eval_list[[1]]
    zval_threshold <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                                 tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio)
    #print(c(zval_binary_mle2,ve_est2,allocation_ratio))
    return(c(zval_binary_mle2,zval_threshold,people_per_ratio[,3],people_per_ratio[,1]))
  }
  saveRDS(all_reps,paste0('storage/timetrend',i,j,'.Rds'))
  zval_binary_mle2 <- all_reps[,1]
  zval_threshold <- all_reps[,2]
  
  #t1e[i] <- sum(zval_binary_mle2<0.05,na.rm=T)/sum(!is.na(zval_binary_mle2))
  #pval_binary_mle[,i] <- zval_binary_mle2
  #t1e1[i] <- sum(zval_binary_mle21<0.05,na.rm=T)/sum(!is.na(zval_binary_mle21))
  #pval_binary_mle1[,i] <- zval_binary_mle21
  #print(c(i,t1e[i],mean(ve_est2),sd(ve_est2)))
  #return(zval_threshold)
  return(c(sum(dnorm(zval_binary_mle2)<0.05,na.rm=T)/sum(!is.na(zval_binary_mle2)), 
           sum(zval_binary_mle2>zval_threshold,na.rm=T)/sum(!is.na(zval_binary_mle2))))
  #hist(rpois(1000,mean(counts-1))+1)
  #hist(counts)
}
saveRDS(t1elist,'storage/t1es.Rds')
t1elist <- readRDS('storage/t1es.Rds')
t1e <- sapply(t1elist,function(x)x[1])[1:length(rates)]
t1e3 <- sapply(t1elist,function(x)x[2])[1:length(rates)]
power <- sapply(t1elist,function(x)x[1])[1:length(rates)+length(rates)]
power3 <- sapply(t1elist,function(x)x[2])[1:length(rates)+length(rates)]
cols <- c('darkorange2','navyblue','hotpink','grey','turquoise')


pdf('figures/trendC19.pdf',height=5,width=15); 
#x11(height=5,width=15); 
par(mar=c(5,5,2,2),mfrow=c(1,3))
matplot(sapply(rates,function(x)1:130*x - 130*x),typ='l',col=cols,lwd=3,lty=1,xlab='Day',ylab='Background rate',cex.lab=1.5,cex.axis=1.5,frame=F)
plot(-rates,t1e,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Type 1 error',ylim=c(0.0,max(t1e,t1e3)),xaxt='n')
axis(1,-rates,-rates,cex.axis=1.5)
points(-rates,t1e3,typ='p',cex=2,pch=15,col=cols)
legend(x=-rates[1],max(t1e,t1e3),bty='n',legend=c('Without correction','With correction'),col=cols[1],pch=c(19,15),cex=1.5)
plot(-rates,power,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Power',ylim=c(0.5,1),xaxt='n')
points(-rates,power3,typ='p',cex=2,pch=15,col=cols)
axis(1,-rates,-rates,cex.axis=1.5)
#legend(x=-rates[1],0.75,bty='n',legend=c('Method 2','Method 6','Method 2 corrected'),col=cols[1],pch=c(19,17,15))
dev.off()



## plot early stopping for allocation probabilities

powers <- lapply(1:2,function(j)
  sapply(1:length(rates),function(i){
    all_reps <- readRDS(paste0('storage/timetrend',i,j,'.Rds'))
    zval_binary_mle2 <- all_reps[,1]
    zval_threshold <- all_reps[,2]
    
    c(
      sum(apply(all_reps,1,function(x)x[3]>0.99))/nrow(all_reps),
      sum(apply(all_reps,1,function(x)x[4]>0.99))/nrow(all_reps),
      sum(apply(all_reps,1,function(x)x[1]>x[2]))/nrow(all_reps),
      sum(apply(all_reps,1,function(x)x[1]>x[2]|any(x[3:4]>0.99)))/nrow(all_reps)
      )
  })
)

cols <- c('darkorange2','navyblue','hotpink','grey')

legend_text <- c(paste0('Day ',eval_day),paste0('Day ',eval_day*2),'p value','All')
pdf('figures/trendearlystopping.pdf',width=10,height=5)
#x11();
par(mar=c(5,5,2,2),mfrow=c(1,2))
for(j in 2:1){
  matplot(t(powers[[j]]),typ='l',col=cols,lwd=3,lty=1,xaxt='n',ylim=c(0,c(0.1,1)[j]),ylab=c('Type 1 error','Power')[j],xlab='Background rate',cex.lab=1.5,cex.axis=1.5,frame=F)
  axis(1,at=1:5,labels=-rates,cex.axis=1.5)
}
legend(x=1,y=.1,col=rev(cols),lwd=3,bty='n',legend=rev(legend_text),cex=1.25)
dev.off()



ss <- lapply(1:2,function(j)
  sapply(1:length(rates),function(i){
    all_reps <- readRDS(paste0('storage/timetrend',i,j,'.Rds'))
    zval_binary_mle2 <- all_reps[,1]
    zval_threshold <- all_reps[,2]
    apply(all_reps,1,function(x){
      if(x[3]>0.99) x[5]
      else if(x[4]>0.99) x[6]
      else x[6]/(2*eval_day)*100
    })
  })
)

print(lapply(1:2,function(j)
  sapply(1:length(rates),function(i){
    all_reps <- readRDS(paste0('storage/timetrend',i,j,'.Rds'))
    zval_binary_mle2 <- all_reps[,1]
    zval_threshold <- all_reps[,2]
    sum(apply(all_reps,1,function(x){
      if(x[3]>0.99) x[5]
      else if(x[4]>0.99) x[6]
      else x[6]/(2*eval_day)*100
    }))/nrow(all_reps)
  })
))
