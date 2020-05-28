source('set_up_script.R')

nIter <<- 1000
range_informative_clusters <<- 20:100
draws <<- 500
cores <<- 32
registerDoParallel(cores=cores)

grid_plot <- function(mat){
  
  grid_zval <- mat[nrow(mat):1,]
  get.pal=colorRampPalette(brewer.pal(9,"Spectral"))
  redCol=rev(get.pal(10))
  bkT <- seq(max(grid_zval,na.rm=T)+1e-10, 0-1e-10,length=length(redCol)+1)
  cex.lab <- 1.5
  maxval <- round(bkT[1],digits=1)
  col.labels<- c(0,maxval/2,maxval)
  cellcolors <- vector()
  for(ii in 1:length(unlist(grid_zval)))
    if(!is.na(grid_zval[ii]))
      cellcolors[ii] <- redCol[tail(which(unlist(grid_zval[ii])<bkT),n=1)]
  pdf(paste0('figures/es',type,'.pdf')); par(mar=c(6,6,2,6))
  color2D.matplot(grid_zval,cellcolors=cellcolors,main="",xlab=paste0("Total weight"),ylab=paste0("Case weight"),cex.lab=1,axes=F,border=NA)
  fullaxis(side=2,las=1,at=0:nrow(grid_zval),labels=round(y_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
  fullaxis(side=1,las=2,at=0:ncol(grid_zval),labels=round(x_labs),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
  color.legend(ncol(grid_zval)+0.5,0,ncol(grid_zval)+1,nrow(grid_zval),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
  dev.off()
  
  
}


compute_grid <- function(type){
  
  results_list <- list()
  vaccinees <- trial_participants <- c()
  for(iter in 1:nIter){
    ## select random person to start
    first_infected <- sample(g_name[eligible_first_person],1)
    netwk <- simulate_contact_network(first_infected,start_day=iter,end_time=eval_day,from_source=0,cluster_flag=0,direct_VE=direct_VE,individual_recruitment_times = T,spread_wrapper = covid_spread_wrapper)
    results_list[[iter]] <- netwk[[1]]
    vaccinees[iter] <- netwk[[4]]
    trial_participants[iter] <- netwk[[5]]
    
  }
  
  #profvis({
  par_results <- do.call(rbind,mclapply(1:draws,function(cl){
    #number_sampled <- sample(range_informative_clusters,1)
    clusters_sampled <- sample(1:nIter,300,replace=F)
    unlisted <- do.call(rbind,lapply(1:length(clusters_sampled),function(x)cbind(results_list[[clusters_sampled[x]]],x)))
    up_to <- unlisted$x[which(unlisted$inTrial)[ceiling(first_threshold*2)]]
    
    cumulative_participants <- trial_participants[clusters_sampled]
    
    case_weight <- 0
    trial_participants2 <- trial_participants[clusters_sampled[1:up_to]]
    results <- results_list[clusters_sampled[1:up_to]]
    vaccinees2 <- vaccinees[clusters_sampled[1:up_to]]
    while(case_weight<first_threshold){
      eval_list <- get_efficacious_probabilities(results,vaccinees2,trial_participants2,contact_network=-1)
      case_weight <- sum(eval_list[[3]])
      if(case_weight<first_threshold){
        trial_participants2 <- trial_participants[clusters_sampled[1:(length(results)+2)]]
        vaccinees2 <- vaccinees[clusters_sampled[1:(length(results)+2)]]
        results <- results_list[clusters_sampled[1:(length(results)+2)]]
      }else{
        zval  <- calculate_zval(eval_list[[3]],eval_list[[2]])
      }
    }
    tp <- sum(trial_participants2)
    
    up_to <- unlisted$x[which(unlisted$inTrial)[ceiling(second_threshold*length(results)/first_threshold)]]
    while(case_weight<second_threshold|!exists('zval2')){
      eval_list2 <- get_efficacious_probabilities(results,vaccinees2,trial_participants2,contact_network=-1)
      case_weight <- sum(eval_list2[[3]])
      if(case_weight<second_threshold){
        trial_participants2 <- trial_participants[clusters_sampled[1:(length(results)+2)]]
        vaccinees2 <- vaccinees[clusters_sampled[1:(length(results)+2)]]
        results <- results_list[clusters_sampled[1:(length(results)+2)]]
      }else{
        zval2  <- calculate_zval(eval_list2[[3]],eval_list2[[2]])
      }
    }
    
    return(c(zval,sum(eval_list[[3]]),sum(eval_list[[2]]),
             zval2,sum(eval_list2[[3]]),sum(eval_list2[[2]]),
             tp,sum(trial_participants2),eval_list[[3]][1],eval_list2[[3]][1])) ## output weights 
  },mc.cores=cores))
  #})
  results_list <- c()
  return(par_results)
  #colnames(par_results) <- c('zval','case_weight','weight')
  #saveRDS(par_results,paste0('storage/',type,'_results.Rds'))
}


first_thresholds <- seq(13,25,by=5)
second_thresholds <- seq(31,41,by=5)

## power ############################################################
#direct_VE <<- 0
#type <- 't1e'
#res1 <- compute_grid(type)
#sum(res1[,1]<0.03)
#sum(res1[,4]<0.02)
#sum(res1[,1]<0.03|res[,1]>0.03&res[,4]<0.02)

types <- c('t1e','power')
for(ty in 1:length(types)){
  type <- types[ty]
  direct_VE <<- c(0,0.7)[ty]
  powers <- halfways <- ss <- matrix(0,nrow=length(first_thresholds),ncol=length(second_thresholds))
  total_iterations <- 100
  results <- c()
  for(ti in 1:total_iterations){
    for(i in 1:length(first_thresholds)){
      first_threshold <<- first_thresholds[i]
      for(j in 1:length(second_thresholds)){
        second_threshold <<- second_thresholds[j]
        res2 <- compute_grid(type)
        fst <- sum(res2[,1]>qnorm(1-0.03))
        snd <- sum(res2[,4]>qnorm(1-0.02))
        total <- sum(res2[,1]>qnorm(1-0.03)|res2[,1]<qnorm(1-0.03)&res2[,4]>qnorm(1-0.02))
        halfways[i,j] <- halfways[i,j] + fst/draws/total_iterations
        powers[i,j] <- powers[i,j] + total/draws/total_iterations
        sample_size <- res2[,6]
        sample_size[res2[,1]>qnorm(1-0.03)] <- res2[res2[,1]>qnorm(1-0.03),3]
        ss[i,j] <- ss[i,j] + sum(sample_size)/draws/total_iterations
        results <- rbind(results,res2)
      }
    }
    print(c(ti))
  }
  print(halfways)
  print(powers)
  print(ss)
  saveRDS(results,paste0('storage/es',type,'results.Rds'))
  saveRDS(list(halfways,powers,ss),paste0('storage/run_es_',type,'.Rds'))
  results <- readRDS(paste0('storage/es',type,'results.Rds'))
  resultsdf <- as.data.frame(results)
  colnames(resultsdf) <- c('earlyzval','earlyweight','V3','latezval','lateweight','V6','earlyss','latess','earlyexpweight','lateexpweight')
  bounds <- c(14,16,18,20)
  bounds2 <- c(37,39,41)
  earlycaseweightboundary <- 8
  ## power
  print(sapply(2:length(bounds),function(x){
    subtab <- subset(resultsdf,earlyweight<bounds[x]&earlyweight>bounds[x-1])
    sapply(2:length(bounds2),function(y){
      subtab2 <- subset(subtab,lateweight<bounds2[y]&lateweight>bounds2[y-1])
      sum(subtab2$earlyzval>qnorm(1-0.03)|
            (subtab2$earlyzval<qnorm(1-0.03)&subtab2$latezval>qnorm(1-0.02)&subtab2$earlyexpweight<earlycaseweightboundary))/nrow(subtab2)
    })
  }))
  ## early power
  print(sapply(2:length(bounds),function(x){
    subtab <- subset(resultsdf,earlyweight<bounds[x]&earlyweight>bounds[x-1])
    sapply(2:length(bounds2),function(y){
      subtab2 <- subset(subtab,lateweight<bounds2[y]&lateweight>bounds2[y-1])
      sum(subtab2$earlyzval>qnorm(1-0.03))/nrow(subtab2)
    })
  }))
  ## expected sample size
  print(sapply(2:length(bounds),function(x){
    subtab <- subset(resultsdf,earlyweight<bounds[x]&earlyweight>bounds[x-1])
    sapply(2:length(bounds2),function(y){
      subtab2 <- subset(subtab,lateweight<bounds2[y]&lateweight>bounds2[y-1])
      sample_size <- subtab2$latess
      sample_size[subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary] <- 
        subtab2$earlyss[subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary]
      mean(sample_size)
    })
  }))
  
  ## sd sample size
  print(sapply(2:length(bounds),function(x){
    subtab <- subset(resultsdf,earlyweight<bounds[x]&earlyweight>bounds[x-1])
    sapply(2:length(bounds2),function(y){
      subtab2 <- subset(subtab,lateweight<bounds2[y]&lateweight>bounds2[y-1])
      sample_size <- subtab2$latess
      sample_size[subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary] <- 
        subtab2$earlyss[subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary]
      sd(sample_size)
    })
  }))
  
  ## early and late sample sizes for y=4
  print(sapply(2:length(bounds),function(x){
    subtab <- subset(resultsdf,earlyweight<bounds[x]&earlyweight>bounds[x-1])
    y <- 2
    subtab2 <- subset(subtab,lateweight<bounds2[y]&lateweight>bounds2[y-1])
    early_index <- subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary
    earlies <- subtab2$earlyss[early_index]
    lates <- subtab2$latess[!early_index]
    c(mean(earlies),mean(lates))
  }))
  
  ## early and late sample size sds for y=4
  print(sapply(2:length(bounds),function(x){
    subtab <- subset(resultsdf,earlyweight<bounds[x]&earlyweight>bounds[x-1])
    y <- 2
    subtab2 <- subset(subtab,lateweight<bounds2[y]&lateweight>bounds2[y-1])
    early_index <- subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary
    earlies <- subtab2$earlyss[early_index]
    lates <- subtab2$latess[!early_index]
    c(sd(earlies),sd(lates))
  }))
  
  ## stopping for futility
  print(sapply(2:length(bounds),function(x){
    subtab <- subset(resultsdf,earlyweight<bounds[x]&earlyweight>bounds[x-1])
    sapply(2:length(bounds2),function(y){
      subtab2 <- subset(subtab,lateweight<bounds2[y]&lateweight>bounds2[y-1])
      sum(subtab2$earlyexpweight>earlycaseweightboundary)/nrow(subtab2)
    })
  }))
  ## stopping erroneously for futility
  print(sapply(2:length(bounds),function(x){
    subtab <- subset(resultsdf,earlyweight<bounds[x]&earlyweight>bounds[x-1])
    sapply(2:length(bounds2),function(y){
      subtab2 <- subset(subtab,lateweight<bounds2[y]&lateweight>bounds2[y-1])
      sum(subtab2$earlyexpweight>earlycaseweightboundary&subtab2$latezval>qnorm(1-0.02))/nrow(subtab2)
      #sum(subtab2$earlyexpweight>10)/nrow(subtab2)
    })
  }))
  
  
  ## without interim
  bounds2 <- c(30:40)
  ## power
  print(
    sapply(2:length(bounds2),function(y){
      subtab2 <- subset(resultsdf,lateweight<bounds2[y]&lateweight>bounds2[y-1])
      sum(subtab2$latezval>qnorm(1-0.05))/nrow(subtab2)
  }))
  y <- 7
  subtab2 <- subset(resultsdf,lateweight<bounds2[y]&lateweight>bounds2[y-1])
  print(mean(subtab2$latess))
  print(sd(subtab2$latess))
  
}
