source('set_up_script.R')

nIter <- 10000
range_informative_clusters <- 10:60
draws <- 100000

## type 1 error ############################################################
direct_VE <<- 0
netwk_list <- list()

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
pvals <- c()
xs <- c()
ys <- c()

##!! we should be estimating VE in the loop over cl.
trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=0)
tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
trial_summary <- c()
ntwks <- unique(tte$cluster)
registerDoParallel(cores=4)

#profvis({
par_results <- do.call(rbind,lapply(1:10,function(cl){
  number_sampled <- sample(range_informative_clusters,1)
  clusters_sampled <- sample(ntwks,number_sampled,replace=F)
  
  pv <- iterate_ph_model(netwk_list[clusters_sampled])
  trial_summary <- list()
  for(i in 1:length(clusters_sampled)) trial_summary[[i]] <- summarise_trial(netwk_list[[clusters_sampled[i]]],ve_est_temp=pv[2])
  #netwk_list <- c()
  tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
  trial_summary <- c()
  ttemat <- as.matrix(tte)
  
  survmodel <- coxph(Surv(time,outcome)~vaccinated,weights=weight,tte)
  vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
  zval <- survmodel$coefficient/sqrt(survmodel$var)
  pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
  return(c(pval,sum(ttemat[ttemat[,4]==T,2]),sum(ttemat[,2]) ))
}))
#})
netwk_list <- c()
pvals <- par_results[,1]
ys <- par_results[,2]
xs <- par_results[,3]
par_results <- ttemat <- tte <- c()
pvals[is.na(pvals)] <- 1
x_lower <- floor(min(xs))
y_lower <- floor(min(ys))
x_upper <- ceiling(max(xs))
y_upper <- ceiling(max(ys))
x_points <- 10
y_points <- 10
x_labs <- seq(x_lower,x_upper,length.out=x_points+1)
y_labs <- seq(y_lower,y_upper,length.out=y_points+1)

binned_x <- discretize(xs,"equalwidth", x_points)
binned_y <- discretize(ys,"equalwidth", y_points)

grid_pval <- matrix(NA,nrow=x_points,ncol=y_points)
for(i in 1:x_points){
  for(j in 1:y_points){
    grid_pvals <- pvals[binned_x==i&binned_y==j]
    if(length(grid_pvals)>0&length(grid_pvals)>8){
      grid_pval[i,j] <- sum(grid_pvals<0.05)/length(grid_pvals)
    }
  }
}

saveRDS(list(binned_x,binned_y,pvals,x_labs,y_labs),'t1e.Rds')

grid_pval <- grid_pval[nrow(grid_pval):1,]
get.pal=colorRampPalette(brewer.pal(9,"Spectral"))
redCol=rev(get.pal(10))
bkT <- seq(max(grid_pval,na.rm=T)+1e-10, min(grid_pval,na.rm=T)-1e-10,length=length(redCol)+1)
cex.lab <- 1.5
maxval <- round(bkT[1],digits=1)
col.labels<- c(0,maxval/2,maxval)
cellcolors <- vector()
for(ii in 1:length(unlist(grid_pval)))
  if(!is.na(grid_pval[ii]))
    cellcolors[ii] <- redCol[tail(which(unlist(grid_pval[ii])<bkT),n=1)]
pdf('phtype1error.pdf'); par(mar=c(6,6,2,6))
color2D.matplot(grid_pval,cellcolors=cellcolors,main="",xlab="Total exposure",ylab="Case exposure",cex.lab=1,axes=F,border=NA)
fullaxis(side=2,las=1,at=0:nrow(grid_pval),labels=round(y_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
fullaxis(side=1,las=2,at=0:ncol(grid_pval),labels=round(x_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
color.legend(ncol(grid_pval)+0.5,0,ncol(grid_pval)+1,nrow(grid_pval),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
dev.off()


binned_x<-binned_y<-pvals<-xs<-ys<-c()

sort(sapply(ls(),function(x)object.size(get(x))))/10000000


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
}

pvals <- c()
xs <- c()
ys <- c()

##!! we should be estimating VE in the loop over cl.
trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=0.8)
tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
ntwks <- unique(tte$cluster)

#profvis({
par_results <- do.call(rbind,mclapply(1:draws,function(cl){
  number_sampled <- sample(range_informative_clusters,1)
  clusters_sampled <- sample(ntwks,number_sampled,replace=F)
  
  pv <- iterate_ph_model(netwk_list[clusters_sampled])
  trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=pv[2])
  netwk_list <- c()
  tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
  trial_summary <- c()
  ttemat <- as.matrix(tte)
  
  survmodel <- coxph(Surv(time,outcome)~vaccinated,weights=weight,tte)
  vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
  zval <- survmodel$coefficient/sqrt(survmodel$var)
  pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
  return(c(pval,sum(ttemat[ttemat[,4]==T,2]),sum(ttemat[,2]) ))
},mc.cores=32))
#})
netwk_list <- tte <- c()
pvals <- par_results[,1]
ys <- par_results[,2]
xs <- par_results[,3]
par_results <- c()
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
saveRDS(list(binned_x,binned_y,pvals,x_labs,y_labs),'power.Rds')

grid_pval <- grid_pval[nrow(grid_pval):1,]
get.pal=colorRampPalette(brewer.pal(9,"Spectral"))
redCol=rev(get.pal(10))
bkT <- seq(max(grid_pval,na.rm=T)+1e-10, min(grid_pval,na.rm=T)-1e-10,length=length(redCol)+1)
cex.lab <- 1.5
maxval <- round(bkT[1],digits=1)
col.labels<- c(0,maxval/2,maxval)
cellcolors <- vector()
for(ii in 1:length(unlist(grid_pval)))
  if(!is.na(grid_pval[ii]))
    cellcolors[ii] <- redCol[tail(which(unlist(grid_pval[ii])<bkT),n=1)]
pdf('phpower.pdf'); par(mar=c(6,6,2,6))
color2D.matplot(grid_pval,cellcolors=cellcolors,main="",xlab="Total exposure",ylab="Case exposure",cex.lab=1,axes=F,border=NA)
fullaxis(side=2,las=1,at=0:nrow(grid_pval),labels=round(y_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
fullaxis(side=1,las=2,at=0:ncol(grid_pval),labels=round(x_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
color.legend(ncol(grid_pval)+0.5,0,ncol(grid_pval)+1,nrow(grid_pval),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
dev.off()
