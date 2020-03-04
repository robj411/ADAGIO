set.seed(1)
source('set_up_script.R')

nIter <- 10000
draws <- 10000
cores <<- 16
registerDoParallel(cores=cores)

highrisk_scalar_samples <<- runif(draws,0,1)
neighbour_scalar_samples <<- runif(draws,0,1)
recall_samples <<- rbeta(draws,5,1)
variables <- list(highrisk_scalar_samples,neighbour_scalar_samples,recall_samples)

number_sampled <- 100#sample(range_informative_clusters,1)

true_contact_list <<- contact_list
true_contact_of_contact_list <<- contact_of_contact_list
true_household_list <<- household_list
true_high_risk_list <<- high_risk_list

trim_contact_networks <- function(recall){
  contact_list <<- trim_contact_network(true_contact_list,recall)
  contact_of_contact_list <<- trim_contact_network(true_contact_of_contact_list,recall)
  household_list <<- trim_contact_network(true_household_list,recall)
  high_risk_list <<- trim_contact_network(true_high_risk_list,recall)
}

trim_contact_network <- function(full_contact_list,recall){
  trimmed_contact_list <- full_contact_list
  for(i in 1:length(full_contact_list)){
    size <- length(full_contact_list[[i]])
    if(size>0){
      ssize <- rbinom(1,size,recall)
      if(ssize>0)
        trimmed_contact_list[[i]] <- sample(full_contact_list[[i]],ssize)
    }
  }
  return(trimmed_contact_list)
}

## type 1 error ############################################################
direct_VE <<- 0
netwk_list <- list()

neighbour_scalar <<- 0.39
highrisk_scalar <<- 2.17

for(iter in 1:nIter){
  ## select random person to start
  first_infected <- sample(g_name,1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
  hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
  inf_time <- min(inf_period,hosp_time)
  netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,inf_time,start_day=iter,from_source=0,cluster_flag=0)
  netwk_list[[iter]] <- netwk
  
}

trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=0)
tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
trial_summary <- c()
ntwks <- unique(tte$cluster)

#profvis({
par_results <- do.call(rbind,mclapply(1:draws,function(cl){
  
  clusters_sampled <- sample(ntwks,number_sampled,replace=F)
  
  neighbour_scalar <<- qbeta(neighbour_scalar_samples[cl],3,3)
  high_risk_scalar <<- 1/qbeta(highrisk_scalar_samples[cl],3,3)
  recall <- recall_samples[cl]
  trim_contact_networks(recall)
  
  pv <- iterate_ph_model(netwk_list[clusters_sampled])
  trial_summary <- list()
  #for(i in 1:length(clusters_sampled)) trial_summary[[i]] <- summarise_trial(netwk_list[[clusters_sampled[i]]],ve_est_temp=pv[2])
  #netwk_list <- c()
  for(cluster in 1:length(clusters_sampled)) {
    trial_summary[[cluster]] <- summarise_trial(netwk_list[[clusters_sampled[cluster]]],ve_est_temp=pv[2])
    #if(!is.null(trial_summary[[cluster]]))
    #  trial_summary[[cluster]] <- cbind(trial_summary[[cluster]],cluster)
  }
  
  tte <- do.call(bind_rows,trial_summary)
  #tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
  trial_summary <- c()
  ttemat <- as.matrix(tte)
  
  survmodel <- coxph(Surv(time,outcome)~vaccinated,weights=weight,tte)
  vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
  zval <- survmodel$coefficient/sqrt(survmodel$var)
  pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
  return(c(pval,sum(ttemat[,7]==T),nrow(ttemat),vaccEffEst[1], sum(ttemat[,2]), sum(ttemat[,6]) )) ## output weights and exposures (time)
},mc.cores=cores))
#})
#netwk_list <- c()

variables <- list(highrisk_scalar_samples,neighbour_scalar_samples,recall_samples,par_results[,5],par_results[,6],par_results[,2])

saveRDS(list(par_results,variables),'storage/mit1e.Rds')



t1e <- par_results[,1]
VEt1e <- par_results[,4]

outcomes <- list(par_results[,1],par_results[,4])
evppi <- sapply(variables,function(x)
  sapply(outcomes,
         function(y)mutinformation(discretize(x),discretize(y))
  )
)
colnames(evppi) <- c('High-risk scalar','Neighbour scalar','Contact recall','Total weight','Total exposure','Cases')
rownames(evppi) <- c('p value (VE=0)','VE estimate (VE=0)')
print(evppi)

x=readRDS('storage/mit1e.Rds')
t1e=x[[1]][,1]
VEt1e=x[[1]][,4]
variables<-x[[2]]

{pdf('figures/outcomebidensitiest1e.pdf',width=15); par(mfrow=c(2,5),mar=c(5,5,2,2))
  for(i in 1:length(outcomes)) 
    for(j in 1:length(variables))
      plot(variables[[j]],outcomes[[i]],frame=F,xlab=colnames(evppi)[j],ylab=rownames(evppi)[i],main='',cex.lab=1.5,cex.axis=1.5,pch=15)
  dev.off()}


## power ############################################################
direct_VE <<- 0.8
netwk_list <- list()

neighbour_scalar <<- 0.39
high_risk_scalar <<- 2.17

contact_list <<- true_contact_list
contact_of_contact_list <<- true_contact_of_contact_list
household_list <<- true_household_list
high_risk_list <<- true_high_risk_list

set.seed(1)
for(iter in 1:nIter){
  ## select random person to start
  first_infected <- sample(g_name,1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
  hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
  inf_time <- min(inf_period,hosp_time)
  netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,inf_time,start_day=iter,from_source=0,cluster_flag=0,direct_VE = direct_VE)
  netwk_list[[iter]] <- netwk
  
}

trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=0.8)
tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
trial_summary <- c()
ntwks <- unique(tte$cluster)


#profvis({
par_results <- do.call(rbind,mclapply(1:draws,function(cl){
  
  clusters_sampled <- sample(ntwks,number_sampled,replace=F)
  neighbour_scalar <<- qbeta(neighbour_scalar_samples[cl],3,3)
  high_risk_scalar <<- 1/qbeta(highrisk_scalar_samples[cl],3,3)
  recall <- recall_samples[cl]
  trim_contact_networks(recall)
  
  pv <- iterate_ph_model(netwk_list[clusters_sampled])
  trial_summary <- list()
  #for(i in 1:length(clusters_sampled)) trial_summary[[i]] <- summarise_trial(netwk_list[[clusters_sampled[i]]],ve_est_temp=pv[2])
  #netwk_list <- c()
  for(cluster in 1:length(clusters_sampled)) {
    trial_summary[[cluster]] <- summarise_trial(netwk_list[[clusters_sampled[cluster]]],ve_est_temp=pv[2])
    #if(!is.null(trial_summary[[cluster]]))
    #  trial_summary[[cluster]] <- cbind(trial_summary[[cluster]],cluster)
  }
  
  tte <- do.call(bind_rows,trial_summary)
  #tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
  trial_summary <- c()
  ttemat <- as.matrix(tte)
  
  survmodel <- coxph(Surv(time,outcome)~vaccinated,weights=weight,tte)
  vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
  zval <- survmodel$coefficient/sqrt(survmodel$var)
  pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
  return(c(pval,sum(ttemat[,7]==T),nrow(ttemat),vaccEffEst[1], sum(ttemat[,2]), sum(ttemat[,6])  )) ## output weights and exposures (time)
},mc.cores=cores))
#})
#netwk_list <- c()

variables <- list(highrisk_scalar_samples,neighbour_scalar_samples,recall_samples,par_results[,5],par_results[,6],par_results[,2])

saveRDS(list(par_results,variables),'storage/mipower.Rds')




power <- par_results[,1]
VE <- par_results[,4]
outcomes2 <- list(power,VE)

evppi2 <- sapply(variables,function(x)
  sapply(outcomes2,
         function(y)mutinformation(infotheo::discretize(x),infotheo::discretize(y))/infotheo::entropy(infotheo::discretize(y))
  )
)
colnames(evppi2) <- c('High-risk scalar','Neighbour scalar','Contact recall','Total weight','Total exposure','Cases')
rownames(evppi2) <- c('p value (VE=0.8)','VE estimate (VE=0.8)')
outcomes3 <- list(t1e,VEt1e,power,VE)
print(evppi2)
evppi3 <- rbind(evppi,evppi2)*100
print(evppi3)

sapply(1:4,function(i)infotheo::entropy(infotheo::discretize(outcomes3[[i]])))

#1/(c(var(old_par_results[,1]),var(old_par_results[,4]))/c(var(par_results[,1]),var(par_results[,4])))

{pdf('figures/scalarmi.pdf',height=6,width=2+length(outcomes3)); 
#    {x11(height=6,width=2+length(outcomes3)); 
  par(mar=c(15,12,3.5,10))
      outcome_labels <- rownames(evppi3)
  labs <- rev(colnames(evppi3))
  get.pal=colorRampPalette(brewer.pal(9,"Reds"))
  redCol=rev(get.pal(12))
  bkT <- seq(max(evppi3)+1e-10, 0,length=13)
  cex.lab <- 1.5
  maxval <- round(bkT[1],digits=1)
  col.labels<- c(0,maxval/2,maxval)
  cellcolors <- vector()
  for(ii in 1:length(unlist(evppi3)))
    cellcolors[ii] <- redCol[tail(which(unlist(evppi3[ii])<bkT),n=1)]
  color2D.matplot(evppi3,cellcolors=cellcolors,main="",xlab="",ylab="",cex.lab=2,axes=F,border='white')
  fullaxis(side=2,las=2,at=0:(length(outcomes3)-1)+1/2,labels=rev(outcome_labels),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.2)
  fullaxis(side=1,las=2,at=(length(labs)-1):0+0.5,labels=labs,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.2)
  mtext(3,text='Mutual information',line=1,cex=1.3)
  color.legend(length(labs)+0.5,0,length(labs)+0.8,length(outcomes3),col.labels,rev(redCol),gradient="y",cex=1.2,align="rb")
  for(i in seq(0,length(outcomes3),by=1)) abline(h=i)
  for(i in seq(0,length(labs),by=1)) abline(v=i)
  dev.off()
}


{pdf('figures/outcomedensities.pdf'); par(mfrow=c(2,2),mar=c(3,5,5,2))
for(i in 1:length(outcomes3)) plot(density(outcomes3[[i]]),frame=F,lwd=2,xlab='',ylab='Density',main=rownames(evppi3)[i],cex.lab=1.5,cex.axis=1.5)
dev.off()
}

dis_fun <- function (x, numBins, r = range(x), b = NULL){
  if(is.null(b)) b = seq(from = r[1], to = r[2], length.out = numBins + 1)
  y = table(cut(x, breaks = b, include.lowest = TRUE))
  return(list(y,b))
}


x=readRDS('storage/mipower.Rds')
power=x[[1]][,1]
VE=x[[1]][,4]
outcomes3 <- list(power,VE)

{pdf('figures/outcomebidensitiespower.pdf',width=15); par(mfrow=c(2,5),mar=c(5,5,2,2))
for(i in 1:length(outcomes3)) 
  for(j in 1:length(variables))
    plot(variables[[j]],outcomes3[[i]],frame=F,xlab=colnames(evppi3)[j],ylab=rownames(evppi3)[i],main='',cex.lab=1.5,cex.axis=1.5,pch=15)
dev.off()}

#old_power <- old_par_results[,1]
#numBins <- length(unlist(unique(infotheo::discretize(old_power))))
#bins <- dis_fun(old_power,numBins=numBins)[[2]]
#entropy(dis_fun(old_power,numBins=numBins,b=bins)[[1]])
#entropy(dis_fun(power,numBins=numBins,b=bins)[[1]])
