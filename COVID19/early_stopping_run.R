source('set_up_script.R')

nIter <<- 1000
range_informative_clusters <<- 20:100
draws <<- 1000
cores <<- 32
registerDoParallel(cores=cores)

compute_grid <- function(type){
  
  results_list <- list()
  vaccinees <- trial_participants <- c()
  for(iter in 1:nIter){
    ## select random person to start
    first_infected <- sample(g_name[eligible_first_person],1)
    netwk <- simulate_contact_network(first_infected,start_day=iter,from_source=0,cluster_flag=0,direct_VE=direct_VE,individual_recruitment_times = T,spread_wrapper = covid_spread_wrapper)
    results_list[[iter]] <- netwk[[1]]
    vaccinees[iter] <- netwk[[4]]
    trial_participants[iter] <- netwk[[5]]
    
  }
  
  #profvis({
  par_results <- do.call(rbind,mclapply(1:draws,function(cl){
    #number_sampled <- sample(range_informative_clusters,1)
    clusters_sampled <- sample(1:nIter,150,replace=F)
    
    cumulative_participants <- trial_participants[clusters_sampled]
    
    up_to <- which(cumsum(cumulative_participants)>first_threshold)[1]
    case_weight <- 0
    trial_participants2 <- trial_participants[clusters_sampled[1:up_to]]
    results <- results_list[clusters_sampled[1:up_to]]
    vaccinees2 <- vaccinees[clusters_sampled[1:up_to]]
    while(case_weight<15){
      eval_list <- get_efficacious_probabilities(results,vaccinees2,trial_participants2,contact_network=-1)
      case_weight <- sum(eval_list[[3]])
      pval  <- calculate_pval(eval_list[[3]],eval_list[[2]])
      trial_participants2 <- trial_participants[clusters_sampled[1:(length(results)+2)]]
      vaccinees2 <- vaccinees[clusters_sampled[1:(length(results)+2)]]
      results <- results_list[clusters_sampled[1:(length(results)+2)]]
    }
    
    up_to <- which(cumsum(cumulative_participants)>second_threshold)[1]
    if(is.na(up_to)) up_to <- length(cumulative_participants)
    trial_participants2 <- trial_participants[clusters_sampled[1:up_to]]
    results <- results_list[clusters_sampled[1:up_to]]
    vaccinees2 <- vaccinees[clusters_sampled[1:up_to]]
    eval_list2 <- get_efficacious_probabilities(results,vaccinees2,trial_participants2,contact_network=-1)
    pval2  <- calculate_pval(eval_list2[[3]],eval_list2[[2]])
    
    return(c(pval,sum(eval_list[[3]]),sum(eval_list[[2]]),pval2,sum(eval_list2[[3]]),sum(eval_list2[[2]]))) ## output weights 
  },mc.cores=cores))
  #})
  results_list <- c()
  return(par_results)
  #colnames(par_results) <- c('pval','case_weight','weight')
  #saveRDS(par_results,paste0('storage/',type,'_results.Rds'))
}
results <- c()


first_thresholds <- seq(800,1200,by=100)
second_thresholds <- seq(1800,2200,by=100)

## power ############################################################
direct_VE <<- 0.7
type <- 'power'
powers <- halfways <- matrix(0,nrow=length(first_thresholds),ncol=length(second_thresholds))
total_interations <- 1000
for(ti in 1:total_interations){
  for(i in 1:length(first_thresholds)){
    first_threshold <<- first_thresholds[i]
    for(j in 1:length(second_thresholds)){
      second_threshold <<- second_thresholds[j]
      res2 <- compute_grid(type)
      fst <- sum(res2[,1]<0.03)
      snd <- sum(res2[,4]<0.02)
      total <- sum(res2[,1]<0.03|res2[,1]>0.03&res2[,4]<0.02)
      halfways[i,j] <- halfways[i,j] + fst/draws/total_interations
      powers[i,j] <- powers[i,j] + total/draws/total_interations
      results <- rbind(results,c(first_threshold,second_threshold,fst,total))
    }
  }
  print(c(ti))
}
print(halfways)
print(powers)

saveRDS(list(halfways,powers),'storage/run_es.Rds')

## type 1 error ############################################################
direct_VE <<- 0
type <- 't1e'
res1 <- compute_grid(type)
sum(res1[,1]<0.03)
sum(res1[,4]<0.02)
sum(res1[,1]<0.03|res[,1]>0.03&res[,4]<0.02)


