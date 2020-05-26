source('set_up_script.R')

nIter <<- 100000
range_informative_clusters <<- 20:125
draws <<- 1000000
cores <<- 32
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
  
  binned_x <- discretize(xs,"equalfreq", x_points)
  binned_y <- discretize(ys,"equalfreq", y_points)
  
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
  
  saveRDS(list(binned_x,binned_y,pvals,x_labs,y_labs),paste0('storage/',type,metric,'.Rds'))
  
  grid_pval <- grid_pval[nrow(grid_pval):1,]
  get.pal=colorRampPalette(brewer.pal(9,"Spectral"))
  redCol=rev(get.pal(10))
  bkT <- seq(max(grid_pval,na.rm=T)+1e-10, 0-1e-10,length=length(redCol)+1)
  cex.lab <- 1.5
  maxval <- round(bkT[1],digits=1)
  col.labels<- c(0,maxval/2,maxval)
  cellcolors <- vector()
  for(ii in 1:length(unlist(grid_pval)))
    if(!is.na(grid_pval[ii]))
      cellcolors[ii] <- redCol[tail(which(unlist(grid_pval[ii])<bkT),n=1)]
  pdf(paste0('figures/ph',type,metric,'.pdf')); par(mar=c(6,6,2,6))
      color2D.matplot(grid_pval,cellcolors=cellcolors,main="",xlab=paste0("Total ",metric),ylab=paste0("Case ",metric),cex.lab=1,axes=F,border=NA)
      fullaxis(side=2,las=1,at=0:nrow(grid_pval),labels=round(y_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
      fullaxis(side=1,las=2,at=0:ncol(grid_pval),labels=round(x_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
      color.legend(ncol(grid_pval)+0.5,0,ncol(grid_pval)+1,nrow(grid_pval),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
      dev.off()
      
      
}

compute_grid <- function(type){
  
  results_list <- list()
  vaccinees <- trial_participants <- c()
  for(iter in 1:nIter){
    ## select random person to start
    first_infected <- sample(g_name[eligible_first_person],1)
    netwk <- simulate_contact_network(first_infected,start_day=iter,from_source=0,cluster_flag=0,direct_VE=direct_VE,individual_recruitment_times = T,spread_wrapper = covid_spread_wrapper,end_time=eval_day)
    results_list[[iter]] <- netwk[[1]]
    vaccinees[iter] <- netwk[[4]]
    trial_participants[iter] <- netwk[[5]]
    
  }
  
  #profvis({
  par_results <- do.call(rbind,mclapply(1:draws,function(cl){
    number_sampled <- sample(range_informative_clusters,1)
    clusters_sampled <- sample(1:nIter,number_sampled,replace=F)
    
    results <- results_list[clusters_sampled]
    vaccinees2 <- vaccinees[clusters_sampled]
    trial_participants2 <- trial_participants[clusters_sampled]
    eval_list <- get_efficacious_probabilities(results,vaccinees2,trial_participants2,contact_network=-1)
    zval  <- calculate_zval(eval_list[[3]],eval_list[[2]])
    
    return(c(zval,sum(eval_list[[3]]),sum(eval_list[[2]]))) ## output weights 
  },mc.cores=cores))
  #})
  results_list <- c()
  colnames(par_results) <- c('zval','case_weight','weight')
  saveRDS(par_results,paste0('storage/',type,'_results.Rds'))
  metric <- 'weight'
  grid_plot(type,metric,par_results)
  
  sort(sapply(ls(),function(x)object.size(get(x))))/10000000
}

## type 1 error ############################################################
#direct_VE <<- 0
#type <- 't1e'
#compute_grid(type)

## power ############################################################
direct_VE <<- 0.7
type <- 'power'
compute_grid(type)


## plot #########################################################
source('../plotting_functions.R')
