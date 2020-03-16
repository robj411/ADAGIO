source('set_up_script.R')


## start ############################################################
nIter <- 1000
e_order <- list()
#profvis({
for(iter in 1:nIter){
  ## select random person to start
  first_infected <- sample(g_name[eligible_first_person],1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,start_day=iter,from_source=0,cluster_flag=0)
  
  results_list[[iter]] <- netwk[[1]]
  cluster_size[iter] <- netwk[[2]]
  #recruit_times[iter] <- netwk[[3]]
  e_order[[iter]] <- netwk[[8]][!duplicated(netwk[[8]])]
}
#})

## results #######################################################

# match the generation interval (time between infection events in an infector-infectee pair) 5.20 (3.78, 6.78) 1.72 (0.91, 3.93)
# and serial interval (time between symptom onsets in an infector-infectee pair) 5.21 (-3.35, 13.94) 4.32 (4.06, 5.58)
##!! nb only picking up the first interval, not the average interval
gi <- unlist(sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  if(nrow(results)>1){
    second <- which(results$InfectedNode==e_order[[iter]][2])
    first <- which(results$InfectedNode==e_order[[iter]][1])
    results$DayInfected[second] - results$DayInfected[first]
  }
}
))
mean(gi)
sd(gi)

si <- unlist(sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  if(nrow(results)>1){
    second <- which(results$InfectedNode==e_order[[iter]][2])
    first <- which(results$InfectedNode==e_order[[iter]][1])
    results$DayInfectious[second] - results$DayInfectious[first]
  }
}
))
mean(si)
sd(si)


