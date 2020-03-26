source('set_up_script.R')
registerDoParallel(cores=10)
## can we infer a trend? ##################################################

direct_VE <- 0.0
reps <- 1000
nIter <- 100
adaptation <- 'TST'
pval_binary_mle2 <- pval_binary_mle21 <- ve_est2 <- ve_est21 <- pval_threshold <- c()
eval_day <- 31
latest_infector_time <- eval_day - 0
func <- get_efficacious_probabilities
rates <- -seq(5e-7,5e-6,by=1e-6)
t1e <- t1e1 <- c()
nClusters <- nIter
pval_binary_mle <- pval_binary_mle1 <- matrix(0,nrow=reps,ncol=length(rates))
t1elist <- foreach(i = rep(1:length(rates),2),j=rep(1:2,each=length(rates))) %dopar% { #for(i in 1:length(rates)){
  direct_VE <- c(0,0.8)[j]
  per_time_step <- rates[i]
  base_rate <- - 130 * rates[i]
  for(rep in 1:reps){
    #profvis({
    allocation_ratio <- 0.5
    results_list <- list()
    vaccinees <- trial_participants <- people_per_ratio <- c()
    vaccinees2 <- trial_participants2 <- randomisation_ratios <- c()
    infectious_by_vaccine <- excluded <- matrix(0,nrow=nIter,ncol=2)
    for(iter in 1:nIter){
      randomisation_ratios[iter] <- allocation_ratio
      ## select random person to start
      first_infected <- sample(g_name,1)
      inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
      # Add them to e_nodes and remove from s_nodes and v_nodes
      #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
      hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
      inf_time <- min(inf_period,hosp_time)
      netwk <- simulate_contact_network(first_infected,inf_time,end_time=eval_day,start_day=iter,from_source=per_time_step,
                                        cluster_flag=0,allocation_ratio=allocation_ratio,direct_VE=direct_VE,base_rate=base_rate)
      
      results_list[[iter]] <- netwk[[1]]
      cluster_size[iter] <- netwk[[2]]
      recruit_times[iter] <- netwk[[3]]
      results <- results_list[[iter]]
      vax <- results$vaccinated
      too_early <- results$DayInfectious<results$RecruitmentDay+10
      infectious_by_vaccine[iter,] <- c(sum(vax&!too_early),sum(!vax&results$inTrial&!too_early))
      excluded[iter,] <- c(sum(vax&too_early),sum(!vax&results$inTrial&too_early))
      
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
      if(adaptation!=''&&iter %% eval_day == 0){
        probs <- func(results_list,vaccinees2,trial_participants2,max_time=length(results_list))
        pop_sizes2 <- probs[[2]]
        fails <- probs[[3]]
        allocation_ratio <- response_adapt(fails,pop_sizes2,days=iter,adaptation=adaptation)
        people_per_ratio <- rbind(people_per_ratio,c(sum(trial_participants2),iter,allocation_ratio))
        #0.9^(iter/nIter)/(0.9^(iter/nIter)+0.1^(iter/nIter))#
      }
    }
    #print(allocation_ratio)
    #})
    # method 2: binary
    pop_sizes <- c(sum(vaccinees2),sum(trial_participants2) - sum(vaccinees2)) - colSums(excluded)
    pval_binary_mle2[rep]  <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
    ve_est2[rep] <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
    pval_threshold[rep] <- trend_robust_function(results_list,vaccinees2,trial_participants2,tested=F,randomisation_ratios=randomisation_ratios,
                                            adaptation=adaptation)
    # method 6: weight non events
    eval_list <- get_efficacious_probabilities2(results_list,vaccinees,trial_participants)
    pval_binary_mle21[rep]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    ve_est21[rep]  <- eval_list[[1]]
    #print(c(pval_binary_mle2,ve_est2,allocation_ratio))
  }
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
legend(x=-rates[1],0.16,bty='n',legend=c('Method 2','Method 6','Method 2 corrected'),col=cols[1],pch=c(19,17,15),cex=1.5)
plot(-rates,power,typ='p',cex=2,pch=19,col=cols,frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Background rate',ylab='Power',ylim=c(0.3,0.8),xaxt='n')
points(-rates,power2,typ='p',cex=2,pch=17,col=cols)
points(-rates,power3,typ='p',cex=2,pch=15,col=cols)
axis(1,-rates,-rates,cex.axis=1.5)
#legend(x=-rates[1],0.75,bty='n',legend=c('Method 2','Method 6','Method 2 corrected'),col=cols[1],pch=c(19,17,15))
dev.off()
