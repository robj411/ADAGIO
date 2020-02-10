source('set_up_script.R')


## type 1 error ############################################################
direct_VE <<- 0
netwk_list <- list()
nIter <- 10000
for(iter in 1:nIter){
  ## select random person to start
  first_infected <- sample(g_name,1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
  hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
  inf_time <- min(inf_period,hosp_time)
  netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,inf_time,start_day=iter,from_source=0,cluster_flag=0)
  netwk_list[[iter]] <- netwk
  results_list[[iter]] <- netwk[[1]]
  cluster_size[iter] <- netwk[[2]]
  recruit_times[iter] <- netwk[[3]]
  hosp_times[iter] <- inf_time
  
}

trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=ves[1])

tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))

ntwks <- unique(tte$cluster)
sapply(ntwks,function(x)sum(tte$time[tte$cluster==x]))
sum(tte$time)

pvals <- c()
xs <- c()
ys <- c()
for(cl in 1:10000){
  # sample between 20 and 500
  number_sampled <- sample(20:200,1)
  clusters_sampled <- sample(ntwks,number_sampled,replace=F)
  subtte <- subset(tte,cluster%in%clusters_sampled)
  survmodel <- coxph(Surv(time,outcome)~vaccinated,weights=weight,subtte)
  vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
  zval <- survmodel$coefficient/sqrt(survmodel$var)
  pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
  pvals[cl] <- pval
  ys[cl] <- sum(subset(subtte,outcome==T)$time)
  xs[cl] <- sum(subtte$time)
}
pvals[is.na(pvals)] <- 1
x_lower <- floor(min(xs))
y_lower <- floor(min(ys))
x_upper <- ceiling(max(xs))
y_upper <- ceiling(max(ys))
x_points <- 10
y_points <- 10
grid_pval <- matrix(NA,nrow=x_points,ncol=y_points)
x_labs <- seq(x_lower,x_upper,length.out=x_points)
y_labs <- seq(y_lower,y_upper,length.out=y_points)

binned_x <- discretize(xs,"equalwidth", x_points)
binned_y <- discretize(ys,"equalwidth", y_points)
for(i in 1:x_points){
  for(j in 1:y_points){
    grid_pvals <- pvals[binned_x==i&binned_y==j]
    if(length(grid_pvals)>0){
      grid_pval[i,j] <- sum(grid_pvals<0.05)/length(grid_pvals)
    }
  }
}

grid_pval <- grid_pval[nrow(grid_pval):1,]
get.pal=colorRampPalette(brewer.pal(9,"RdBu"))
redCol=rev(get.pal(5))
bkT <- seq(max(grid_pval,na.rm=T)+1e-10, min(grid_pval,na.rm=T)-1e-10,length=length(redCol)+1)
cex.lab <- 1.5
maxval <- round(bkT[1],digits=1)
col.labels<- c(0,maxval/2,maxval)
cellcolors <- vector()
for(ii in 1:length(unlist(grid_pval)))
  if(!is.na(grid_pval[ii]))
    cellcolors[ii] <- redCol[tail(which(unlist(grid_pval[ii])<bkT),n=1)]
pdf('phtype1error.pdf')
color2D.matplot(grid_pval,cellcolors=cellcolors,main="",xlab="Sample size",ylab="Events",cex.lab=1,axes=F,border=NA)
fullaxis(side=2,las=1,at=1:nrow(grid_pval),labels=round(y_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
fullaxis(side=1,las=2,at=1:ncol(grid_pval),labels=round(x_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
color.legend(ncol(grid_pval)+0.5,0,ncol(grid_pval)+2,nrow(grid_pval),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
dev.off()

## power ############################################################
direct_VE <<- 0.8
netwk_list <- list()
for(iter in 1:nIter){
  ## select random person to start
  first_infected <- sample(g_name,1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
  hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
  inf_time <- min(inf_period,hosp_time)
  netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,inf_time,start_day=iter,from_source=0,cluster_flag=0,direct_VE=direct_VE)
  netwk_list[[iter]] <- netwk
  results_list[[iter]] <- netwk[[1]]
  cluster_size[iter] <- netwk[[2]]
  recruit_times[iter] <- netwk[[3]]
  hosp_times[iter] <- inf_time
  
}

trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=ves[1])

tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))

ntwks <- unique(tte$cluster)
sapply(ntwks,function(x)sum(tte$time[tte$cluster==x]))
sum(tte$time)

pvals <- c()
xs <- c()
ys <- c()
for(cl in 1:10000){
  # sample between 20 and 500
  number_sampled <- sample(20:200,1)
  clusters_sampled <- sample(ntwks,number_sampled,replace=F)
  subtte <- subset(tte,cluster%in%clusters_sampled)
  survmodel <- coxph(Surv(time,outcome)~vaccinated,weights=weight,subtte)
  vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
  zval <- survmodel$coefficient/sqrt(survmodel$var)
  pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
  pvals[cl] <- pval
  ys[cl] <- sum(subset(subtte,outcome==T)$time)
  xs[cl] <- sum(subtte$time)
}
pvals[is.na(pvals)] <- 1
x_lower <- floor(min(xs))
y_lower <- floor(min(ys))
x_upper <- ceiling(max(xs))
y_upper <- ceiling(max(ys))
x_points <- 10
y_points <- 10
grid_pval <- matrix(NA,nrow=x_points,ncol=y_points)
x_labs <- seq(x_lower,x_upper,length.out=x_points)
y_labs <- seq(y_lower,y_upper,length.out=y_points)

binned_x <- discretize(xs,"equalwidth", x_points)
binned_y <- discretize(ys,"equalwidth", y_points)
for(i in 1:x_points){
  for(j in 1:y_points){
    grid_pvals <- pvals[binned_x==i&binned_y==j]
    if(length(grid_pvals)>0){
      grid_pval[i,j] <- sum(grid_pvals<0.05)/length(grid_pvals)
    }
  }
}

grid_pval <- grid_pval[nrow(grid_pval):1,]
get.pal=colorRampPalette(brewer.pal(9,"RdBu"))
redCol=rev(get.pal(5))
bkT <- seq(max(grid_pval,na.rm=T)+1e-10, min(grid_pval,na.rm=T)-1e-10,length=length(redCol)+1)
cex.lab <- 1.5
maxval <- round(bkT[1],digits=1)
col.labels<- c(0,maxval/2,maxval)
cellcolors <- vector()
for(ii in 1:length(unlist(grid_pval)))
  if(!is.na(grid_pval[ii]))
    cellcolors[ii] <- redCol[tail(which(unlist(grid_pval[ii])<bkT),n=1)]
pdf('phpower.pdf')
color2D.matplot(grid_pval,cellcolors=cellcolors,main="",xlab="Sample size",ylab="Events",cex.lab=1,axes=F,border=NA)
fullaxis(side=2,las=1,at=1:nrow(grid_pval),labels=round(y_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
fullaxis(side=1,las=2,at=1:ncol(grid_pval),labels=round(x_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
color.legend(ncol(grid_pval)+0.5,0,ncol(grid_pval)+2,nrow(grid_pval),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
dev.off()
