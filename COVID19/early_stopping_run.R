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
  
  #profvis({
  par_results <- do.call(rbind,mclapply(1:draws,function(cl){
    #set.seed(cl*total_iterations + seed)
    #number_sampled <- sample(range_informative_clusters,1)
    clusters_sampled <- sample(1:nIter,300,replace=F)
    unlisted <- do.call(rbind,lapply(1:length(clusters_sampled),function(x)cbind(results_list[[clusters_sampled[x]]],x)))
    result_tab <- subset(unlisted,inTrial&DayInfectious>RecruitmentDay)
    n_infected <- sapply(1:length(clusters_sampled),function(x)sum(result_tab$x==x))
    n_v_infected <- sapply(1:length(clusters_sampled),function(x)sum(result_tab$x==x&result_tab$vaccinated))
    raw_weights <- get_infectee_weights(results=result_tab,ve_point_est=0,contact_network=-1,tested=F,correct_for_ve=F)
    result_tab$raw_weight <- rowSums(raw_weights[[1]])
    x <- result_tab$raw_weight
    y <- 1 - x
    up_to <- result_tab$x[ceiling(first_threshold)]
    reordered_participants <- trial_participants[clusters_sampled] - n_infected
    reordered_vaccinees <- vaccinees[clusters_sampled] - n_v_infected
    case_weight <- 0
    trial_participants2 <- reordered_participants[1:up_to]
    vaccinees2 <- reordered_vaccinees[1:up_to]
    n_infected2 <- n_infected[1:up_to]
    vax_flag <- result_tab$vaccinated[result_tab$x<=up_to]
    x_up_to <- x[result_tab$x<=up_to]
    y_up_to <- 1 - x_up_to
    ve_estimate <- c(0,1)
    break_count <- 0
    while(case_weight<first_threshold){
      while(abs(ve_estimate[1]-ve_estimate[2])>0.005&&break_count<5){
        ve_estimate[2] <- ve_estimate[1]
        new_weights <- x_up_to[vax_flag]
        new_weights <- (1-ve_estimate[1])*new_weights/(y_up_to[vax_flag]+(1-ve_estimate[1])*new_weights)
        fails <- c(sum(new_weights*rbinom(length(new_weights),1,observed)),sum(x_up_to[!vax_flag]*rbinom(sum(!vax_flag),1,observed)))
        pop_sizes2 <- c(sum(vaccinees2),sum(trial_participants2)-sum(vaccinees2)) + fails
        if(fails[2]>0&&!any(pop_sizes2==0))
          ve_estimate[1] <- calculate_ve(fails,sizes=pop_sizes2)
        break_count <- break_count + 1
      }
      case_weight <- sum(fails)
      if(case_weight<first_threshold){
        up_to <- up_to + 1
        trial_participants2 <- reordered_participants[1:up_to]
        vaccinees2 <- reordered_vaccinees[1:up_to]
        x_up_to <- x[result_tab$x<=up_to]
        y_up_to <- 1 - x_up_to
        vax_flag <- result_tab$vaccinated[result_tab$x<=up_to]
        break_count <- 0
        ve_estimate[2] <- ve_estimate[1]- 0.006
      }else{
        zval  <- calculate_zval(fails,sizes=pop_sizes2)
      }
    }
    early_ss <- sum(pop_sizes2)
    early_case <- sum(fails)
    early_fails <- fails
    tp <- sum(trial_participants2)
    #eval_list <- get_efficacious_probabilities(unlisted,vaccinees2,trial_participants2,contact_network = -1,observed = observed)
    
    up_to <- ceiling(second_threshold*up_to/first_threshold*0.7)
    trial_participants2 <- reordered_participants[1:up_to]
    vaccinees2 <- reordered_vaccinees[1:up_to]
    n_infected2 <- n_infected[1:up_to]
    vax_flag <- result_tab$vaccinated[result_tab$x<=up_to]
    x_up_to <- x[result_tab$x<=up_to]
    y_up_to <- 1 - x_up_to
    ve_estimate <- c(0,1)
    break_count <- 0
    while(case_weight<second_threshold|!exists('zval2')){
      while(abs(ve_estimate[1]-ve_estimate[2])>0.005&&break_count<5){
        ve_estimate[2] <- ve_estimate[1]
        new_weights <- x_up_to[vax_flag]
        new_weights <- (1-ve_estimate[1])*new_weights/(y_up_to[vax_flag]+(1-ve_estimate[1])*new_weights)
        fails <- c(sum(new_weights*rbinom(length(new_weights),1,observed)),sum(x_up_to[!vax_flag]*rbinom(sum(!vax_flag),1,observed)))
        pop_sizes2 <- c(sum(vaccinees2),sum(trial_participants2)-sum(vaccinees2)) + fails
        if(fails[2]>0&&!any(pop_sizes2==0))
          ve_estimate[1] <- calculate_ve(fails,sizes=pop_sizes2)
        break_count <- break_count + 1
      }
      case_weight <- sum(fails)
      if(case_weight<second_threshold){
        up_to <- up_to + 1
        trial_participants2 <- reordered_participants[1:up_to]
        vaccinees2 <- reordered_vaccinees[1:up_to]
        x_up_to <- x[result_tab$x<=up_to]
        y_up_to <- 1 - x_up_to
        vax_flag <- result_tab$vaccinated[result_tab$x<=up_to]
        break_count <- 0
        ve_estimate[2] <- ve_estimate[1] - 0.006
      }else{
        zval2  <- calculate_zval(fails,sizes=pop_sizes2)
      }
    }
    
    return(c(zval,early_case,early_ss,
             zval2,sum(fails),sum(pop_sizes2),
             tp,sum(trial_participants2),early_fails[1],fails[1])) ## output weights 
  },mc.cores=cores))
  #})
  results_list <- c()
  return(par_results)
  #colnames(par_results) <- c('zval','case_weight','weight')
  #saveRDS(par_results,paste0('storage/',type,'_results.Rds'))
}

bounds <- c(12,14,16,24)
bounds2 <- c(24,30,32)

first_thresholds <- bounds
second_thresholds <- bounds2


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
  total_iterations <<- 100
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
  trial_participants <<- trial_participants
  vaccinees <<- vaccinees
  results_list <<- results_list
  
  powers <- halfways <- ss <- matrix(0,nrow=length(first_thresholds),ncol=length(second_thresholds))
  results <- c()
  for(ti in 1:total_iterations){
    print(c(ti))
    for(i in 1:length(first_thresholds)){
      first_threshold <<- first_thresholds[i]
      for(j in 1:length(second_thresholds)){
        #print(j)
        seed <<- ti
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
  }
  print(halfways)
  print(powers)
  print(ss)
  saveRDS(results,paste0('storage/es',type,'results.Rds'))
  saveRDS(list(halfways,powers,ss),paste0('storage/run_es_',type,'.Rds'))
  results <- readRDS(paste0('storage/es',type,'results.Rds'))
  resultsdf <- as.data.frame(results)
  colnames(resultsdf) <- c('earlyzval','earlyweight','V3','latezval','lateweight','V6','earlyss','latess','earlyexpweight','lateexpweight')
  earlycaseweightboundary <- 8
  res_by_threshold <- list()
  for(i in 1:length(first_thresholds)){
    res_by_threshold[[i]] <- list()
    for(j in 1:length(second_thresholds)){
      res_by_threshold[[i]][[j]] <- matrix(nrow=0,ncol=ncol(resultsdf))
    }
  }
  nrows <- draws
  start_index <- 1
  for(ti in 1:total_iterations){
    for(i in 1:length(first_thresholds)){
      for(j in 1:length(second_thresholds)){
        res_by_threshold[[i]][[j]] <- rbind(res_by_threshold[[i]][[j]],resultsdf[start_index:(start_index+nrows-1),])
        start_index <- start_index + nrows
      }
    }
  }
  
  
  ## power
  final_power <- matrix(sapply(1:length(bounds),function(x){
    sapply(1:length(bounds2),function(y){
      subtab2 <- res_by_threshold[[x]][[y]]
      sum(subtab2$earlyzval>qnorm(1-0.03)|
            (subtab2$earlyzval<qnorm(1-0.03)&subtab2$latezval>qnorm(1-0.02)&
               subtab2$earlyexpweight<earlycaseweightboundary))/nrow(subtab2)
    })
  }),nrow=length(bounds2))
  print(final_power)
  ## early power
  early_power <- matrix(sapply(1:length(bounds),function(x){
    sapply(1:length(bounds2),function(y){
      subtab2 <- res_by_threshold[[x]][[y]]
      sum(subtab2$earlyzval>qnorm(1-0.03))/nrow(subtab2)
    })
  }),nrow=length(bounds2))
  print(early_power)
  ## expected sample size
  ess <- matrix(sapply(1:length(bounds),function(x){
    sapply(1:length(bounds2),function(y){
      subtab2 <- res_by_threshold[[x]][[y]]
      sample_size <- subtab2$latess
      sample_size[subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary] <- 
        subtab2$earlyss[subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary]
      mean(sample_size)
    })
  }),nrow=length(bounds2))
  print(ess)
  
  ## sd sample size
  sdss <- matrix(sapply(1:length(bounds),function(x){
    sapply(1:length(bounds2),function(y){
      subtab2 <- res_by_threshold[[x]][[y]]
      sample_size <- subtab2$latess
      sample_size[subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary] <- 
        subtab2$earlyss[subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary]
      sd(sample_size)
    })
  }),nrow=length(bounds2))
  print(sdss)
  
  ## early and late sample sizes for y=4
  y <- 3
  elss <- as.matrix(sapply(1:length(bounds),function(x){
    subtab2 <- res_by_threshold[[x]][[y]]
    early_index <- subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary
    earlies <- subtab2$earlyss[early_index]
    lates <- subtab2$latess[!early_index]
    c(mean(earlies),mean(lates))
  }))
  print(elss)
  
  ## early and late sample size sds for y=4
  sdelss <- as.matrix(sapply(1:length(bounds),function(x){
    subtab2 <- res_by_threshold[[x]][[y]]
    early_index <- subtab2$earlyzval>qnorm(1-0.03)|subtab2$earlyexpweight>earlycaseweightboundary
    earlies <- subtab2$earlyss[early_index]
    lates <- subtab2$latess[!early_index]
    c(sd(earlies),sd(lates))
  }))
  print(sdelss)
  
  ## stopping for futility
  futility <- matrix(sapply(1:length(bounds),function(x){
    sapply(1:length(bounds2),function(y){
      subtab2 <- res_by_threshold[[x]][[y]]
      sum(subtab2$earlyexpweight>earlycaseweightboundary)/nrow(subtab2)
    })
  }),nrow=length(bounds2))
  print(futility)
  ## stopping erroneously for futility
  futility_error <- matrix(sapply(1:length(bounds),function(x){
    sapply(1:length(bounds2),function(y){
      subtab2 <- res_by_threshold[[x]][[y]]
      sum(subtab2$earlyexpweight>earlycaseweightboundary&subtab2$latezval>qnorm(1-0.02))/nrow(subtab2)
      #sum(subtab2$earlyexpweight>10)/nrow(subtab2)
    })
  }),nrow=length(bounds2))
  print(futility_error)
  
  
  print_tab <- data.frame(bounds,early_power[3,],futility[3,],final_power[3,])
  print_tab$ess <- paste0(round(ess[3,]),' (',round(sdss[3,]),')')
  print_tab$lss <- paste0(round(elss[1,]),' (',round(sdelss[1,]),')')
  print_tab$uss <- paste0(round(elss[2,]),' (',round(sdelss[2,]),')')
  print(xtable(print_tab,digits=c(0,0,2,2,2,0,0,0)), include.rownames = FALSE)
  
  
  
  ## without interim
  ## power
  print(
    sapply(1:length(bounds2),function(y){
      subtab2 <- res_by_threshold[[y]][[1]]
      sum(subtab2$latezval>qnorm(1-0.05))/nrow(subtab2)
  }))
  y <- 1
  subtab2 <- res_by_threshold[[1]][[1]]
  print(mean(subtab2$latess))
  print(sd(subtab2$latess))
  
}
