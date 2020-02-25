source('set_up_script.R')

nIter <<- 1000
range_informative_clusters <<- 20:100
draws <<- 500
cores <<- 2
registerDoParallel(cores=cores)

grid_plot <- function(type,metric,par_results){
  pvals <- par_results[,1]
  ycol <- which(colnames(par_results)==paste0('case_',metric))
  xcol <- which(colnames(par_results)==paste0('',metric))
  ys <- par_results[,ycol]
  xs <- par_results[,xcol]
  #par_results <- ttemat <- tte <- c()
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
      if(length(grid_pvals)>0){
        grid_pval[i,j] <- sum(grid_pvals<0.05)/length(grid_pvals)
      }
    }
  }
  print(grid_pval)
  
  saveRDS(list(binned_x,binned_y,pvals,x_labs,y_labs),paste0(type,metric,'.Rds'))
  
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
  pdf(paste0('ph',type,metric,'.pdf')); par(mar=c(6,6,2,6))
      color2D.matplot(grid_pval,cellcolors=cellcolors,main="",xlab=paste0("Total ",metric),ylab=paste0("Case ",metric),cex.lab=1,axes=F,border=NA)
      fullaxis(side=2,las=1,at=0:nrow(grid_pval),labels=round(y_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
      fullaxis(side=1,las=2,at=0:ncol(grid_pval),labels=round(x_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
      color.legend(ncol(grid_pval)+0.5,0,ncol(grid_pval)+1,nrow(grid_pval),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
      dev.off()
      
      
}

compute_grid <- function(type){
  
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
  
  trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=0)
  tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
  trial_summary <- c()
  ntwks <- unique(tte$cluster)
  
  #profvis({
  par_results <- do.call(rbind,mclapply(1:draws,function(cl){
    number_sampled <- sample(range_informative_clusters,1)
    clusters_sampled <- sample(ntwks,number_sampled,replace=F)
    
    pv <- iterate_ph_model(netwk_list[clusters_sampled])
    trial_summary <- list()
    #for(i in 1:length(clusters_sampled)) trial_summary[[i]] <- summarise_trial(netwk_list[[clusters_sampled[i]]],ve_est_temp=pv[2])
    #netwk_list <- c()
    for(cluster in 1:length(clusters_sampled)) {
      trial_summary[[cluster]] <- summarise_trial(netwk_list[[clusters_sampled[cluster]]],ve_est_temp=pv[2])
      if(!is.null(trial_summary[[cluster]]))
        cbind(trial_summary[[cluster]],cluster)
    }
    
    tte <- do.call(bind_rows,trial_summary)
    #tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
    trial_summary <- c()
    ttemat <- as.matrix(tte)
    
    survmodel <- coxph(Surv(time,outcome)~vaccinated,weights=weight,tte)
    vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
    zval <- survmodel$coefficient/sqrt(survmodel$var)
    pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
    weight_ind <- which(colnames(tte)=='weight')
    time_ind <- which(colnames(tte)=='time')
    inf_ind <- which(colnames(tte)=='outcome')
    return(c(pval,
             sum(ttemat[ttemat[,inf_ind]==T,weight_ind]),sum(ttemat[,weight_ind]),
             sum(ttemat[ttemat[,inf_ind]==T,time_ind]),sum(ttemat[,time_ind]) )) ## output weights and exposures (time)
  },mc.cores=cores))
  #})
  netwk_list <- c()
  colnames(par_results) <- c('pval','case_weight','weight','case_exposure','exposure')
  
  for(metric in c('weight','exposure')) grid_plot(type,metric,par_results)
  
  sort(sapply(ls(),function(x)object.size(get(x))))/10000000
}

## type 1 error ############################################################
direct_VE <<- 0
type <- 't1e'
compute_grid(type)

## power ############################################################
direct_VE <<- 0.8
type <- 'power'
compute_grid(type)


