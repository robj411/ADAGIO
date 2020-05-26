source('set_up_script.R')


## start ############################################################
nIter <- 1000
e_order <- list()
#profvis({
for(iter in 1:nIter){
  ## select random person to start
  first_infected <- sample(g_name[eligible_first_person],1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  netwk <- simulate_contact_network(first_infected,start_day=iter,from_source=0,cluster_flag=0,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper)
  
  results_list[[iter]] <- netwk[[1]]
  cluster_size[iter] <- netwk[[2]]
  #recruit_times[iter] <- netwk[[3]]
  e_order[[iter]] <- netwk[[8]][!duplicated(netwk[[8]])]
}
#})

###############################################################
## count e infections
covid_spread_wrapper <- function(i_nodes_info,s_nodes,v_nodes,e_nodes_info,direct_VE){
  # to contacts
  # e infects house and work and anyone - only enodes infected one day ago or more, and only enodes with one day left incubating
  ##!! a subset of i_nodes are nonsymptomatic and therefore continue to infect contacts. these should be a fixed list, not sampled randomly every time.
  current_infectious <- c(e_nodes_info[e_nodes_info[,2]>=e_nodes_info[,3],1])
  start_s <- sum(s_nodes)
  if(length(current_infectious)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_list,beta_scalar=nonrandom_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=random_list,beta_scalar=random_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
  }
  end_s <- sum(s_nodes)
  e_count <<- e_count + (start_s - end_s)
  
  ##!! a subset of i_nodes are nonsymptomatic and therefore continue to infect contacts. these should be a fixed list, not sampled randomly every time.
  current_infectious <- c(i_nodes_info[c(runif(nrow(i_nodes_info))<observed),1])
  if(length(current_infectious)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_list,beta_scalar=nonrandom_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=random_list,beta_scalar=random_scalar)
    s_nodes[e_nodes_info[,1]] <- 0
  }
  
  # i infects house
  current_infectious <- c(i_nodes_info[,1])
  start_e <- nrow(e_nodes_info)
  if(length(current_infectious)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=household_list,beta_scalar=1)
  }
  end_e <- nrow(e_nodes_info)
  i_count <<- i_count + (end_e - start_e)
  return(e_nodes_info)
}

e_count <<- 0
e_h_count <<- 0
i_count <<- 0
for(iter in 1:nIter){
  first_infected <- sample(g_name[eligible_first_person],1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  netwk <- simulate_contact_network(first_infected,start_day=iter,from_source=0,cluster_flag=0,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper)
  
  results_list[[iter]] <- netwk[[1]]
  cluster_size[iter] <- netwk[[2]]
}

print(c(e_count,i_count)/(e_count+i_count))
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



#############################################################
## weight comparison

get_efficacious_probabilities_orig <- get_efficacious_probabilities
get_efficacious_probabilities <- function(results_list,vaccinees,trial_participants,max_time=10000,contact_network=2,
                                          tested=F,randomisation_ratios=NULL,rbht_norm=0,people_per_ratio=NULL,adaptation='TST',observed=1,age_counts=NULL){
  infectious_by_vaccine <- excluded <- c()
  for(iter in 1:length(results_list)){
    results <- results_list[[iter]]
    infectious_by_vaccine <- rbind(infectious_by_vaccine,
                                   c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+6)*observed,
                                     sum(!results$vaccinated&results$inTrial&results$DayInfectious>results$RecruitmentDay+6)*observed))
    excluded <- rbind(excluded,c(sum(results$vaccinated&results$DayInfectious<results$RecruitmentDay+7)*observed,
                                 sum(!results$vaccinated&results$inTrial&results$DayInfectious<results$RecruitmentDay+7)*observed))
  }
  weight_sums <- colSums(infectious_by_vaccine,na.rm=T)
  pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
  pval_binary_mle <- calculate_pval(weight_sums,pop_sizes)
  ve_estimate  <- calculate_ve(weight_sums,pop_sizes)
  
  return(list(ve_estimate[1],pop_sizes,weight_sums))
}

# generate results
vaccinees <- trial_participants <- c()
nIter <- 100
direct_VE <- 0.7
results_list <- list()
for(iter in 1:nIter){
  first_infected <- sample(g_name[eligible_first_person],1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  netwk <- simulate_contact_network(first_infected,direct_VE=direct_VE,start_day=iter,from_source=0,cluster_flag=0,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper)
  results_list[[iter]] <- netwk[[1]]
  cluster_size[iter] <- netwk[[2]]
  vaccinees[iter] <- netwk[[4]]
  trial_participants[iter] <- netwk[[5]]
}
# get continuous estimate for vaccine efficacy
eval_list <- get_efficacious_probabilities_orig(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
conve <- eval_list[[1]]
# summarise results
true_vs_weight <- true_vs_weight_con <- c(0,0,0)
for(i in 1:length(results_list)){
  results <- results_list[[i]]
  if(nrow(results)>1){
    trueweight <- sum(results$vaccinated&results$DayInfected>results$RecruitmentDay+1,na.rm=T)
    binweight <- sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+7,na.rm=T)
    conweight <- sum(get_infectee_weights(results=results,ve_point_est=conve,contact_network=-1,tested=F)[[1]][,1])
    true_vs_weight <- rbind(true_vs_weight,c(trueweight,binweight,conweight))
    trueweight <- sum(!results$vaccinated&results$DayInfected>results$RecruitmentDay+1,na.rm=T)
    binweight <- sum(!results$vaccinated&results$DayInfectious>results$RecruitmentDay+7,na.rm=T)
    conweight <- sum(get_infectee_weights(results=results,ve_point_est=conve,contact_network=-1,tested=F)[[1]][,2])
    true_vs_weight_con <- rbind(true_vs_weight_con,c(trueweight,binweight,conweight))
  }
}

sapply(2:3,function(x)sum(abs(true_vs_weight[,1]-true_vs_weight[,x])))
sapply(2:3,function(x)sum(abs(true_vs_weight_con[,1]-true_vs_weight_con[,x])))
sapply(2:3,function(x)sum(abs(true_vs_weight[,1]+true_vs_weight_con[,1]-true_vs_weight[,x]-true_vs_weight_con[,x])))
colSums(true_vs_weight)
colSums(true_vs_weight_con)
colSums(true_vs_weight+true_vs_weight_con)
sum(true_vs_weight[,1]>true_vs_weight[,2])
sum(true_vs_weight[,1]>true_vs_weight[,3])
sum(true_vs_weight[,2]>true_vs_weight[,1])
sum(true_vs_weight[,3]>true_vs_weight[,1])
sum(true_vs_weight_con[,1]>true_vs_weight_con[,3])
sum(true_vs_weight_con[,3]>true_vs_weight_con[,1])
nrow(true_vs_weight)
pdf('figures/count_v_estimate.pdf'); par(mfrow=c(1,1),mar=c(5,5,2,2),cex=1.5,pch=20)
plot(true_vs_weight[,1],true_vs_weight[,3],col='navyblue',ylab='Estimated count',xlab='True count',frame=F,ylim=c(0,max(true_vs_weight)),xlim=c(0,max(true_vs_weight)))
points(true_vs_weight[,1],true_vs_weight[,2],cex=0.8,col='hotpink')
lines(c(0,max(true_vs_weight)),c(0,max(true_vs_weight)),col='grey',lty=2,lwd=2)
legend(x=0,y=22,legend=c('Binary weights','Continuous weights'),col=c('navyblue','hotpink'),pch=20,cex=0.75,bty='n')
dev.off()
1