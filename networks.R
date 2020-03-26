source('set_up_script.R')


## start ############################################################
nIter <- 100
#profvis({
for(iter in 1:nIter){
  ## select random person to start
  first_infected <- sample(g_name,1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
  hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
  inf_time <- min(inf_period,hosp_time)
  netwk <- simulate_contact_network(first_infected,inf_time,start_day=iter,from_source=0,cluster_flag=0)
  
  results_list[[iter]] <- netwk[[1]]
  cluster_size[iter] <- netwk[[2]]
  recruit_times[iter] <- netwk[[3]][1]
  hosp_times[iter] <- inf_time
  
}
#})

## results #######################################################

number_infectious <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster)
}
)

number_infectious_after_randomisation <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfectious>results$RecruitmentDay,na.rm=T)
}
)
number_infected_after_randomisation <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfected>results$RecruitmentDay,na.rm=T)
}
)

number_infectious_before_day10 <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfectious<results$RecruitmentDay+10&results$DayInfectious>results$RecruitmentDay-2,na.rm=T)
}
)
contacts_infectious_before_day10 <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$contact&results$DayInfectious<results$RecruitmentDay+10&results$DayInfectious>results$RecruitmentDay-2,na.rm=T)
}
)
c(sum(contacts_infectious_before_day10),sum(number_infectious_before_day10))
number_infectious_after_randomisation_before_day10 <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10,na.rm=T)
}
)
number_infected_after_randomisation_before_day10 <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfected>results$RecruitmentDay&results$DayInfected<results$RecruitmentDay+10,na.rm=T)
}
)

secondary_infections <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$DayInfected>hosp_times[iter])
}
)
infectious_on_recruitment <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$DayInfectious[-1]<recruit_times[iter])
}
)
infectious_by_vaccine <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+9),sum(results$inTrial&results$DayInfectious>results$RecruitmentDay+9))
}
)

c(sum(contacts_infectious_before_day10),sum(number_infectious_before_day10))
df <- data.frame(col1=0,col2=c(0.29,0.59,0.88),col3=c(0.14,0.57,0.70))
df[1,1] <- sum(number_infectious_before_day10-number_infectious_after_randomisation_before_day10)/sum(cluster_size)*100
df[2,1] <- sum(number_infectious_after_randomisation_before_day10)/sum(cluster_size)*100
df[3,1] <- sum(number_infectious_before_day10)/sum(cluster_size)*100
rownames(df) <- c('Infectious before randomisation','Infectious between days','Infectious before day 10')
print(xtable(df,digits=2))

col1 <- 0:6
col2 <- sapply(0:6,function(x)sum(number_infectious-number_infectious_before_day10==x))/10
col3 <- c(78.7,10.6,4.3,2.1,2.1,0.0,2.1)
col4 <- c(84.4,8.5,3.4,1.7,0.9,0.0,0.9)
print(xtable(data.frame(col1,col2,col3,col4),digits=1), include.rownames = FALSE)


quantile(cluster_size,c(0.25,0.5,0.75))

variables <- list(cs=cluster_size,ht=hosp_times,rt=recruit_times)
outcomes <- list(number_infectious_before_day10,number_infectious_after_randomisation_before_day10,number_infected_after_randomisation_before_day10,secondary_infections,infectious_on_recruitment)
zeros <- sapply(outcomes,function(x)sum(x==0)/length(x))*100
means <- sapply(outcomes,function(x)mean(x))
ranges <- apply(sapply(outcomes,function(x)quantile(x[x>0],c(0.05,0.95))),2,function(x)paste0(x,collapse='--'))
outcome_labels <- c('Number infectious, day < 10','Number infectious, 0 < day < 10',
                    'Number infected, 0 < day < 10','Number of secondary infections','Number infectious on recruitment')
df <- data.frame(cbind(zeros,means,ranges),row.names=outcome_labels,stringsAsFactors = F)
for(i in 1:2) df[,i] <- as.numeric(df[,i])
xtable(df,digits=c(1,1,2,1))

##!! doesn't include the 12 excluded
delay_before_day_ten <- c(rep(0,30+3),4,3,3,3,rep(2,5),rep(1,4))
denominator <- 3096+1461
sum(delay_before_day_ten)/denominator
hist(delay_before_day_ten)
den <- c(rnbinom(1000,3.5,0.1),rep(0,2000))
points(sort(unique(den)),sapply(sort(unique(den)),function(x)sum(den==x))/50)

probability_by_lag <- readRDS(paste0('probability_by_lag_8080.Rds'))
probability_after_day_0 <- readRDS(paste0('probability_after_day_0_8080.Rds'))
probability_by_lag_given_removal <- readRDS(paste0('probability_by_lag_given_removal_8080.Rds'))
probability_after_day_0_given_removal <- readRDS(paste0('probability_after_day_0_given_removal_8080.Rds'))

ref_recruit_day <- 30
true_vs_weight <- c(0,0,0,0,0)
how_surprising <- c()
for(i in 1:length(results_list)){
  results <- results_list[[i]]
  if(nrow(results)>1){
    weights <- get_weighted_results(results,how_surprising)
    true_vs_weight <- rbind(true_vs_weight,weights[[1]])
    how_surprising <- weights[[2]]
  }
}
hist(how_surprising)
quantile(how_surprising,c(0.95,0.975,0.99))
sapply(2:5,function(x)sum(abs(true_vs_weight[,1]-true_vs_weight[,x])))
colSums(true_vs_weight)
nrow(true_vs_weight)
pdf('count_v_estimate.pdf'); par(mfrow=c(1,1),mar=c(5,5,2,2),cex=1.5,pch=20)
plot(true_vs_weight[,1],true_vs_weight[,3],col='navyblue',ylab='Estimated count',xlab='True count',frame=F)
points(true_vs_weight[,1],true_vs_weight[,2],cex=0.8,col='hotpink')
lines(c(0,max(true_vs_weight)),c(0,max(true_vs_weight)),col='grey',lty=2,lwd=2)
legend(x=0,y=22,legend=c('Binary weights','Continuous weights'),col=c('navyblue','hotpink'),pch=20,cex=0.75,bty='n')
dev.off()

plot(sapply(results_list,nrow))


## infections #############################################
## if we know removal day
removal_days <- c()
for(i in 1:length(results_list)){
  results <- results_list[[i]]
  removed <- subset(results,!is.na(DayRemoved))
  removal_days <- c(removal_days,removed$DayRemoved-removed$DayInfectious)
}


rm(results_list)

## mi #################################################################

evppi <- sapply(variables,function(x)
  sapply(outcomes,
         function(y)mutinformation(discretize(x),y)
  )
)

{pdf('evppi.pdf',height=10,width=4+length(outcomes)); 
#  {x11(height=10,width=4+length(outcomes)); 
    par(mar=c(12,24,3.5,10))
  labs <- rev(c('Cluster size','Index removal time','Recruitment time'))
  get.pal=colorRampPalette(brewer.pal(9,"Reds"))
  redCol=rev(get.pal(12))
  bkT <- seq(max(evppi)+1e-10, 0,length=13)
  cex.lab <- 1.5
  maxval <- round(bkT[1],digits=1)
  col.labels<- c(0,maxval/2,maxval)
  cellcolors <- vector()
  for(ii in 1:length(unlist(evppi)))
    cellcolors[ii] <- redCol[tail(which(unlist(evppi[ii])<bkT),n=1)]
  color2D.matplot(evppi,cellcolors=cellcolors,main="",xlab="",ylab="",cex.lab=2,axes=F,border='white')
  fullaxis(side=2,las=2,at=0:(length(outcomes)-1)+1/2,labels=rev(outcome_labels),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.2)
  fullaxis(side=1,las=2,at=(length(labs)-1):0+0.5,labels=labs,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.2)
  mtext(3,text='Mutual information',line=1,cex=1.3)
  color.legend(length(labs)+0.5,0,length(labs)+0.8,length(outcomes),col.labels,rev(redCol),gradient="y",cex=1.2,align="rb")
  for(i in seq(0,length(outcomes),by=1)) abline(h=i)
  for(i in seq(0,length(labs),by=1)) abline(v=i)
  #dev.off()
  }


## abc ##################################
data <- c()
# are there any cases before day 0
data[1] <- 16/11833
# % before day 0 = 16 per 11833 people = 0.0014 per person
data[2] <- 16/11833
# 0 in cluster days days 0 to 10
data[3:7] <- sapply(0:4,function(x)sum(delay_before_day_ten==x)/length(delay_before_day_ten))
# pc infectees are contacts
infected_clusters <- sum(c(11,6,10,21,6))
of_which_contacts <- sum(c(10,6,10,20,6))
data[8] <- of_which_contacts/infected_clusters
# pc infectees are high risk
of_which_hr <- sum(c(10,5,10,20,5))
data[9] <- of_which_hr/infected_clusters
cluster_sum <- sapply(0:4,function(x)sum(delay_before_day_ten==x))

zeros <- c()
sumll <- c()

distnce <- distnce_sample <- c()
distnce[1] <- 50
parameter_samples <- list(beta_hyper=list(vals=c(),star=c()),
                          nb_hyper=list(vals=c(),star=c()),
                          hr_hyper=list(vals=c(),star=c())
)
param_names=c('beta_base','neighbour_scalar','high_risk_scalar')
param_quantiles <- c(qlnorm,qbeta,function(x,...){1/qbeta(x,...)})
param_hypers <- list(c(log(0.01),0.7),c(2,2),c(5,5))
successes <- rep(0,length(parameter_samples))
step_size <- rep(0.05,length(parameter_samples))
for(i in 1:length(parameter_samples)){
  parameter_samples[[i]]$vals[1] <- parameter_samples[[i]]$star <- runif(1)
}
direct_VE <- 0
#profvis({
for(iter in 1:10000){
  
  
  for(param in 1:length(parameter_samples)){
    parameter_sample <- parameter_samples[[param]]
    star <- parameter_sample$star
    val <- parameter_sample$val
    assign(param_names[param], 
           param_quantiles[param][[1]](star,param_hypers[[param]][1],param_hypers[[param]][2]))
  
    
    est <- matrix(0,ncol=10,nrow=length(delay_before_day_ten))
    for(i in 1:length(delay_before_day_ten)){
      ## select random person to start
      first_infected <- sample(g_name,1)
      inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
      # Add them to e_nodes and remove from s_nodes and v_nodes
      #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
      hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
      inf_time <- min(inf_period,hosp_time)
      netwk <- simulate_contact_network(first_infected,inf_time,end_time=10,cluster_flag=cluster_flag)
      
      results <- netwk[[1]]
      cluster_sz <- netwk[[2]]
      #recruit_times[iter] <- netwk[[3]]
      #hosp_times[iter] <- inf_time
      
      # number before day 0
      est[i,1] <- sum(results$inCluster&results$DayInfectious<results$RecruitmentDay+1)
      
      # 0 in cluster days days 0 to 10
      est[i,2:6] <- sapply(0:4,function(x)sum(results$inCluster&results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10)==x)
      # cases in cluster days days 0 to 10
      est[i,7] <- sum(results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10&results$inCluster)
      # pc infectees are contacts
      est[i,8] <- sum(results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10&results$contact)
      # pc infectees are high risk
      est[i,9] <- sum(results$inCluster&results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10&results$highrisk==T)
      # cluster size
      est[i,10] <- cluster_sz
    }
    distnce_tmp <- c()
    # distance from number before day 0
    distnce_tmp[1] <- abs(data[1]-sum(est[,1])/sum(est[,10]))*100
    # distances for number of clusters with 0, 1, 2 etc cases
    distnce_tmp[2:6] <- abs(colSums(est[,2:6]) - cluster_sum)
    distnce_tmp[7] <- abs(sum(est[,8])/sum(est[,7]) - of_which_contacts/infected_clusters)*10
    distnce_tmp[8] <- abs(sum(est[,9])/sum(est[,7]) - of_which_hr/infected_clusters)*10
    distncestar <- sum(distnce_tmp[is.finite(distnce_tmp)])
    #print(c(distnce_tmp,distncestar))
    
    distnce_sample[iter] <- distncestar
    if(distncestar < 30){ 
      val[iter+1] <- star
      distnce[iter+1] <- distncestar
      successes[param] <- successes[param] + 1
    }else{
      val[iter+1] <- val[iter]
      distnce[iter+1] <- distnce[iter]
    }
    star <- val[iter+1] + rnorm(1,0,step_size[param])
    while(star<0) star <- star+ceiling(abs(star))
    while(star>1) star <- star-floor(star)
    
    if(iter%%20==0){
      acceptance_ratio <- successes[param]/20
      step_size[param] <- step_size[param]*((acceptance_ratio+0.05)/0.4)^0.25
      successes[param] <- 0
    }
    parameter_sample$star <- star
    parameter_sample$val <- val
    parameter_samples[[param]] <- parameter_sample
  }
  
  zeros[iter] <- est[2]==0
  sumll[iter] <- sum(distnce_tmp[is.finite(distnce_tmp)])
}
  #})
par(mar=c(4,4,2,2))
hist(qlnorm(runif(1000),log(0.01),0.7))
hist(qlnorm(parameter_samples[[1]]$val,log(0.01),0.7))
mean(qlnorm(parameter_samples[[1]]$val,log(0.01),0.7))
plot(parameter_samples[[1]]$val)
plot(parameter_samples[[2]]$val)
hist(qbeta(runif(1000),2,2))
hist(qbeta(parameter_samples[[2]]$val,2,2))
mean(qbeta(parameter_samples[[2]]$val,2,2))
plot(parameter_samples[[2]]$val,parameter_samples[[3]]$val)
plot(parameter_samples[[3]]$val)
hist(1/qbeta(runif(1000),5,5))
hist(1/qbeta(parameter_samples[[3]]$val,5,5))
mean(1/qbeta(parameter_samples[[3]]$val,5,5))
sum(zeros)
plot(distnce_sample)
