sample_size <- 1000000
exp_time <- pmin(rgamma(sample_size,shape=infperiod_shape,rate=infperiod_rate) , 
                 rtruncnorm(sample_size,a=0,mean=hosp_mean_index,sd=hosp_sd_index))
rec_time <- ceiling(rtruncnorm(sample_size,a=0,mean=recruit_mean,sd=recruit_sd))
time_diff <- exp_time - rec_time
prob_overhang <- sum(time_diff > 0)/sample_size 
overhang <- time_diff[time_diff > 0]
weight <- c()
max_day <- 40
for(inf_day in 1:max_day){
  probability_to_exclude <- sapply(overhang,function(x)pgamma(max(inf_day-x,0),shape=incperiod_shape,rate=incperiod_rate))
  probability_to_zero <- pgamma(inf_day,shape=incperiod_shape,rate=incperiod_rate)
  prob <- (probability_to_zero - probability_to_exclude)/(1-probability_to_exclude)
  #plot(overhang,prob*prob_overhang)
  weight[inf_day] <- mean(prob*prob_overhang)
}
plot(1:max_day,weight)

days <- 80
lags <- 80
recruitment_time <- 30 # ceiling(rtruncnorm(sample_size,a=0,mean=recruit_mean,sd=recruit_sd))
probability_after_day_0 <- probability_by_lag <- matrix(0,nrow=lags,ncol=days)
inc_period <- rgamma(sample_size,shape=incperiod_shape,rate=incperiod_rate)
for(j in 1:days){
  exp_time <- pmin(rgamma(sample_size,shape=infperiod_shape,rate=infperiod_rate),
                   rtruncnorm(sample_size,a=0,mean=hosp_mean,sd=hosp_sd))
  if(j<recruitment_time) exp_time <- pmin(exp_time,recruitment_time)
  #overhang <- exp_time[exp_time+j>recruitment_time] # samples still infectious after recruitment
  for(i in 1:lags) {
    # probability to be infected given a difference in time of infectiousness of i 
    # and time of first person infected relative to recruitment of j-10
    denom <- inc_period<i&inc_period>(i-exp_time)
    probability_by_lag[i,j] <- sum(denom)/sample_size
    # probability infected after day 0 (10) given a difference in time of infectiousness of i 
    # and time of first person infected relative to recruitment of j-10
    if((i+j)>recruitment_time&length(overhang)>0){
      #probability_to_exclude <- sapply(overhang,function(x)pgamma(max(i-x,0),shape=incperiod_shape,rate=incperiod_rate))
      #probability_to_zero <- pgamma(i+j-recruitment_time,shape=incperiod_shape,rate=incperiod_rate)
      #prob <- (probability_to_zero - probability_to_exclude)/(1-probability_to_exclude)
      infectious_after_zero <- exp_time+j>recruitment_time
      infected_after_infectious <- inc_period<i
      infected_while_infectious <- inc_period>(i-exp_time)
      infected_after_zero <- j+i-inc_period > recruitment_time
      probability_after_day_0[i,j] <- sum(infected_after_zero&infected_after_infectious&infected_while_infectious)/
        sum(infected_after_infectious&infected_while_infectious)
      #mean(prob)#sum(exp_time+j>recruitment_time&(j+i-inc_period)>10&inc_period>(i-exp_time))/sample_size/probability_by_lag[i,j]
    }
  }
}
probability_after_day_0[is.na(probability_after_day_0)] <- 0
saveRDS(probability_by_lag,paste0('probability_by_lag_',lags,days,'.Rds'))
saveRDS(probability_after_day_0,paste0('probability_after_day_0_',lags,days,'.Rds'))

{pdf('person_prob.pdf',height=10,width=8); 
  par(mar=c(8,10,3.5,10))
  xlabs <- 1:days-recruitment_time
  ylabs <- 1:lags
  heatplot(mat=probability_by_lag,xlabs=xlabs,ylabs=ylabs,cols="Blues",ncols=9,nbreaks=12,
           title='Probability Person 2 infected by Person 1',
           text1='Day Person 1 becomes infectious',text2='Number of days Person 2 becomes infectious after Person 1')
  dev.off()}
{pdf('day_prob.pdf',height=10,width=8); 
  par(mar=c(8,10,3.5,10))
  xlabs <- 1:days-recruitment_time
  ylabs <- 1:lags
  heatplot(mat=probability_after_day_0,xlabs=xlabs,ylabs=ylabs,cols="Reds",ncols=9,nbreaks=12,
           title='Probability Person 2 infected after\n recruitment (day 0) given infected by Person 1',
           text1='Day Person 1 becomes infectious',text2='Number of days Person 2 becomes infectious after Person 1')
  dev.off()}
{pdf('person_day_prob.pdf',height=10,width=8); 
  par(mar=c(8,10,3.5,10))
  xlabs <- 1:days-recruitment_time
  ylabs <- 1:lags
  heatplot(mat=probability_after_day_0*probability_by_lag,xlabs=xlabs,ylabs=ylabs,cols="Purples",ncols=9,nbreaks=12,
           title='Probability Person 2 infected by Person 1\n after recruitment (day 0)',
           text1='Day Person 1 becomes infectious',text2='Number of days Person 2 becomes infectious after Person 1')
  dev.off()}
heatplot <- function(mat,xlabs,ylabs,cols="Blues",ncols=9,nbreaks=12,title='Heatplot',text1='',text2=''){
  get.pal=colorRampPalette(brewer.pal(ncols,cols))
  redCol=rev(get.pal(nbreaks))
  bkT <- seq(max(mat)+1e-10, 0,length=nbreaks+1)
  cex.lab <- 1.5
  maxval <- round(bkT[1],digits=1)
  col.labels<- c(0,maxval/2,maxval)
  cellcolors <- vector()
  for(ii in 1:length(unlist(mat)))
    cellcolors[ii] <- redCol[tail(which(unlist(mat[ii])<bkT),n=1)]
  color2D.matplot(mat,cellcolors=cellcolors,main="",xlab="",ylab="",cex.lab=2,axes=F,border='white')
  fullaxis(side=2,las=2,at=0:(length(ylabs)-1)+1/2,labels=rev(ylabs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.2)
  fullaxis(side=1,las=2,at=0:(length(xlabs)-1)+0.5,labels=xlabs,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.2)
  mtext(3,text=title,line=1,cex=1.3)
  mtext(1,text=text1,line=3,cex=1.3)
  mtext(2,text=text2,line=3,cex=1.3)
  color.legend(length(xlabs)+0.5,0,length(xlabs)+0.8,length(ylabs),col.labels,rev(redCol),gradient="y",cex=1.2,align="rb")
  for(i in seq(0,length(ylabs),by=1)) abline(h=i)
  for(i in seq(0,length(xlabs),by=1)) abline(v=i)
}

probability_after_day_0_given_removal <- probability_by_lag_given_removal <- list()
inc_period <- rgamma(sample_size,shape=incperiod_shape,rate=incperiod_rate)
for(exp_time in 1:20){
  probability_after_day_0_given_removal[[exp_time]] <- probability_by_lag_given_removal[[exp_time]] <- matrix(0,nrow=lags,ncol=days)
  for(j in 1:days){
    for(i in 1:lags) {
      # probability to be infected given a difference in time of infectiousness of i 
      # and time of first person infected relative to recruitment of j-10
      denom <- inc_period<i&inc_period>(i-exp_time)
      probability_by_lag_given_removal[[exp_time]][i,j] <- sum(denom)/sample_size
      # probability infected after day 0 (10) given a difference in time of infectiousness of i 
      # and time of first person infected relative to recruitment of j-10
      if((i+j)>recruitment_time){
        #probability_to_exclude <- sapply(overhang,function(x)pgamma(max(i-x,0),shape=incperiod_shape,rate=incperiod_rate))
        #probability_to_zero <- pgamma(i+j-recruitment_time,shape=incperiod_shape,rate=incperiod_rate)
        #prob <- (probability_to_zero - probability_to_exclude)/(1-probability_to_exclude)
        infected_after_infectious <- inc_period<i
        infected_while_infectious <- inc_period>(i-exp_time)
        infected_after_zero <- j+i-inc_period > recruitment_time
        probability_after_day_0_given_removal[[exp_time]][i,j] <- sum(infected_after_zero&infected_after_infectious&infected_while_infectious)/
          sum(infected_after_infectious&infected_while_infectious)
      }
    }
  }
}
for(i in 1:length(probability_after_day_0_given_removal)) 
  probability_after_day_0_given_removal[[i]][is.na(probability_after_day_0_given_removal[[i]])] <- 0
saveRDS(probability_by_lag_given_removal,paste0('probability_by_lag_given_removal_',lags,days,'.Rds'))
saveRDS(probability_after_day_0_given_removal,paste0('probability_after_day_0_given_removal_',lags,days,'.Rds'))
