rm(list=ls())
setwd('~/overflow_dropbox/ADAGIO/')
source('Code/functions_network.R')
library(igraph)
library(truncnorm)
library(infotheo)
library(xtable)
library(RColorBrewer)
library(plotrix)
library(profvis)
library(funique)
library(doParallel)
library(foreach)

## build network ###############################################################

number_of_households <- 100
household_sizes <- rnbinom(number_of_households,12.99,0.66)
while(any(household_sizes==0)){
  household_sizes[household_sizes==0] <- rnbinom(sum(household_sizes==0),3.45,1-0.66)
}

# assign individuals to households
label_start <- 0
hh <- list()
for(i in 1:number_of_households) {
  hh[[i]] <- make_full_graph(household_sizes[i]) %>%
    set_vertex_attr("name", value = label_start+1:household_sizes[i])
  label_start <- label_start + household_sizes[i]
}

# extract household data frames
attrs <- do.call(rbind,lapply(hh,function(x)as_data_frame(x,'vertices')))
# combine all
el <- do.call(rbind,lapply(hh,function(x)as_data_frame(x)))
# convert to network
new_g <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)
# save layout for plotting
save_layout <- layout_nicely(new_g)
# add household labels
hh_labels <- rep(1:number_of_households,household_sizes)
new_g <- set_vertex_attr(new_g,'hh',value=hh_labels)

## add index connections
index_person <- sapply(1:number_of_households,function(x)which(hh_labels==x)[1])
index_label <- rep(0,length(V(new_g)))
index_label[index_person] <- 1
new_g <- set_vertex_attr(new_g,'index',value=index_label)
for(i in 1:300) {
  first_person <- sample(index_person,1,replace=F)
  second_person <- sample(index_person[!index_person%in%c(first_person,ego(new_g,order=1,nodes=first_person)[[1]])],1,replace=F)
  new_g <- add_edges(new_g,edges=c(first_person,second_person))
}

## add child connections
young_person <- t(sapply(1:number_of_households,function(x)which(hh_labels==x)[2:4]))
# remove NA from smaller hhs
young_person <- young_person[!is.na(young_person)]
child_label <- rep(0,length(V(new_g)))
child_label[young_person] <- 1
new_g <- set_vertex_attr(new_g,'child',value=child_label)
class_size <- 25
for(i in 1:150) {
  for(j in 1:round(length(young_person)/class_size)){
    max_index <- min(class_size+class_size*(j-1),length(young_person))
    min_index <- 1+class_size*(j-1)
    young_people <- young_person[min_index:max_index]
    first_person <- sample(young_people,1,replace=F)
    second_person <- sample(young_people[!young_people%in%c(first_person,ego(new_g,order=1,nodes=first_person)[[1]])],1,replace=F)
    new_g <- add_edges(new_g,edges=c(first_person,second_person))
  }
}

## add random connections
for(i in 1:1000) {
  first_person <- sample(V(new_g),1)
  first_hh <- V(new_g)$hh[first_person]
  second_person <- sample(V(new_g)[V(new_g)$hh!=first_hh&!V(new_g)$name%in%ego(new_g,order=1,nodes=first_person)[[1]]],1)
  new_g <- add_edges(new_g,edges=c(first_person,second_person))
}
#plot.igraph(new_g,vertex.label=NA,vertex.size=1,layout=save_layout)
#cluster_sizes <- sapply(V(new_g),function(x)ego_size(new_g,order=2,nodes=x))
#hist(cluster_sizes,main='',xlab='Cluster size')
#c(mean(cluster_sizes),quantile(cluster_sizes,c(0.25,0.5,0.75)))

# plot degree distribution - aiming for mean=17.5
degreedistribution <- degree.distribution(new_g)*length(E(new_g))
barplot(degreedistribution,ylab='Number of people', xlab='Number of connections',names.arg=0:(length(degreedistribution)-1),main='')
average_contacts <- sum(degreedistribution*c(1:length(degreedistribution)-1)/length(E(new_g)))
length(E(new_g))/length(V(new_g))*2

# get list of neighbours
contact_list <- lapply(V(new_g),function(x) {cs <- as.vector(unlist(ego(new_g,order=1,nodes=x))); cs[cs!=x]})
mean(sapply(contact_list,length))

## get neighbourhood network
# assume 5 hh per neighbourhood
neighbourhood_sizes <- rep(5,length=floor(number_of_households/5)-1)
neighbourhood_sizes <- c(neighbourhood_sizes,number_of_households-sum(neighbourhood_sizes))
number_of_neighbourhoods <- length(neighbourhood_sizes)
# assume all hh within neighbourhood are connected
rate_within <- 1
within_rates <- diag(nrow=number_of_neighbourhoods,ncol=number_of_neighbourhoods,x=rate_within)
# make connections between hh across neighbourhoods to represent extended family
rate_between <- 0.045
between_rates <- matrix(rate_between,nrow=number_of_neighbourhoods,ncol=number_of_neighbourhoods) -
  diag(nrow=number_of_neighbourhoods,ncol=number_of_neighbourhoods,x=rate_between)
rates <- within_rates+between_rates
# create network
g2 <- sample_sbm(sum(neighbourhood_sizes),rates,neighbourhood_sizes)
median(degree(g2)*8)

## translate into individual-level network with connections between all hh members
neighbour_adjacency_matrix <- matrix(0,nrow=length(V(new_g)),ncol=length(V(new_g)))
# populate adjacency matrix edge by edge
for(i in 1:length(E(g2))) {
  hh_edge <- ends(g2, i, names = F)
  hh1 <- hh_edge[1]
  hh2 <- hh_edge[2]
  hh1_occupants <- V(new_g)$hh==hh1
  hh2_occupants <- V(new_g)$hh==hh2
  neighbour_adjacency_matrix[hh1_occupants,hh2_occupants] <- neighbour_adjacency_matrix[hh2_occupants,hh1_occupants] <- 1
}
neighbourhood_g <- graph_from_adjacency_matrix(neighbour_adjacency_matrix,mode='undirected')
degreedistribution <- degree.distribution(neighbourhood_g)*length(E(neighbourhood_g))
average_contacts <- sum(degreedistribution*c(1:length(degreedistribution)-1)/length(E(neighbourhood_g)))
rm(neighbour_adjacency_matrix)
# aiming for average contacts approx 60
##!! there are almost certainly duplicate edges here, so some people might get two tries to infect someone. Is that what we want?

# get list of neighbours
contact_of_contact_list <- lapply(V(neighbourhood_g),function(x) {cofc <- as.vector(unlist(ego(neighbourhood_g,order=1,nodes=x))); cofc[cofc!=x]})

household_list <- lapply(V(new_g),function(x){hh_members <- which(hh_labels==hh_labels[x]); as.vector(hh_members[hh_members!=x])})

# add high-risk labels, to be used for ring vaccination, could be used to increase disease spread
# assume high risk rate is constant across contacts and contacts of contacts
high_risk_rate <- sum(c(330,171,58,246,574,231))/sum(2151,1435,1104,1678,3796,2572)
high_risk_list <- lapply(V(new_g),function(x){
  sz <- length(unique(c(contact_list[[x]],contact_of_contact_list[[x]])))
  nhr <- rbinom(1,sz,high_risk_rate)
  ct <- contact_list[[x]]
  non_hh_ct <- ct[!ct%in%household_list[[x]]]
  if(nhr>length(contact_list[[x]])&length(non_hh_ct)>0){
    hr <- sample(non_hh_ct,min(nhr-length(ct),length(non_hh_ct)))
  }else{
    hr <- c()
  }
  as.vector(hr)
  })

average_cluster_size <- mean(sapply(1:length(contact_list),
                                    function(x)length(contact_list[[x]])+
                                      length(contact_of_contact_list[[x]])+
                                      ifelse(length(high_risk_list[[x]])==0,0,sapply(high_risk_list[[x]],function(y)length(household_list[[y]])))))
average_cluster_size

## functions #######################################################

source('network_functions.R')
source('evaluation_functions.R')

probability_by_lag <<- readRDS(paste0('probability_by_lag_8080.Rds'))
probability_after_day_0 <<- readRDS(paste0('probability_after_day_0_8080.Rds'))
probability_by_lag_given_removal <<- readRDS(paste0('probability_by_lag_given_removal_8080.Rds'))
probability_after_day_0_given_removal <<- readRDS(paste0('probability_after_day_0_given_removal_8080.Rds'))

## set up #######################################################

# Per-time-step hazard of infection for a susceptible nodes from an infectious
# neighbour
beta <- 0.0065
high_risk_scalar <- 2.17
# fraction of beta applied to neighbours ("contacts of contacts")
neighbour_scalar <- 0.39
# Gamma-distribution parameters of incubation and infectious period and wait times
incperiod_shape <- 3.11
incperiod_rate <- 0.32
infperiod_shape <- 1.13
infperiod_rate <- 0.226
hosp_shape_index <- 2
hosp_rate_index <- 0.5
hosp_shape <- 2
hosp_rate <- 1.5
recruit_shape <- 5.4
recruit_rate <- 0.47
hosp_mean_index <- 3.85
hosp_sd_index <- 2.76
hosp_mean <- 2.8
hosp_sd <- 1.5
vacc_mean <- 1.5
vacc_sd <- 1
recruit_mean <- 10.32
recruit_sd <- 4.79
direct_VE <- 0.0

g <<- new_g

g_name <- as.numeric(as.vector(V(g)$name))
vertices <- V(g)
cluster_size <- hosp_times <- recruit_times <- c()
results_list <- list()

## ring vaccination trial ##################################################
nClusters <- 100
nTrials <- 1000
ves <- c(0,0.9)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- c(0,1)
trial_designs <- expand.grid(VE=ves,cluster=cluster_flags,adapt=adaptations)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
trial_designs <- rbind(trial_designs,trial_designs[trial_designs$adapt=='',])
trial_designs$weight[(nCombAdapt+1):(nComb*(length(adaptations)+1))] <- 'binary'
trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- 0
ref_recruit_day <- 30
pval_binary_mle2 <- pval_binary_mle <- ve_est2 <- ve_est <- c()
registerDoParallel(cores=8)

trial_results <- foreach(des = 1:nCombAdapt) %dopar% {
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  vaccinated_count <- infectious_count <- vaccinated_countb <- infectious_countb <- 0
  for(tr in 1:nTrials){
    vaccinees <- trial_participants <- c()
    infectious_by_vaccine <- weight_hh_rem <- excluded <- matrix(0,nrow=nClusters,ncol=2)
    results_list <- list()
    allocation_ratio <- 0.5
    for(iter in 1:nClusters){
      ## select random person to start
      first_infected <- sample(g_name,1)
      inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
      #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
      hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
      inf_time <- min(inf_period,hosp_time)
      netwk <- simulate_contact_network(beta,neighbour_scalar,high_risk_scalar,first_infected,inf_time,cluster_flag=cluster_flag)
      results_list[[iter]] <- netwk[[1]]
      results <- results_list[[iter]]
      infectious_by_vaccine[iter,] <- c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+9),sum(!results$vaccinated&results$inTrial&results$DayInfectious>results$RecruitmentDay+9))
      excluded[iter,] <- c(sum(results$vaccinated&results$DayInfectious<results$RecruitmentDay+10),sum(!results$vaccinated&results$inTrial&results$DayInfectious<results$RecruitmentDay+10))
      vaccinees[iter] <- netwk[[4]]
      trial_participants[iter] <- netwk[[5]]
      ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
      if(adaptation!=''&&iter %% 31 == 0){
        allocation_ratio <- response_adapt(results_list,vaccinees,trial_participants,adaptation)
      }
    }
    eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants)
    pval_binary_mle2[tr]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
    ve_est2[tr]  <- eval_list[[1]]
    vaccinated_count <- vaccinated_count + sum(vaccinees)/nTrials
    infectious_count <- infectious_count + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
    if(adaptation==''){
      pop_sizes <- c(sum(vaccinees),sum(trial_participants) - sum(vaccinees)) - colSums(excluded)
      pval_binary_mle[tr] <- calculate_pval(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      ve_est[tr]  <- calculate_ve(colSums(infectious_by_vaccine,na.rm=T),pop_sizes)
      vaccinated_countb <- vaccinated_countb + sum(vaccinees)/nTrials
      infectious_countb <- infectious_countb + (sum(sapply(results_list,nrow))-length(results_list))/nTrials
    }
    ## ICC without weighting
    if(cluster_flag==1){
      vax <- vaccinees
      non_vax <- trial_participants - vax
      trial_case <- sapply(results_list,function(x)sum(x$inTrial==T))
      vax_case <- sapply(results_list,function(x)sum(x$vaccinated==T))
      non_vax_case <- trial_case - vax_case
      cid <- rep(1:length(trial_participants),times=trial_participants)
      non_cases <- trial_participants - trial_case
      y <- unlist(sapply(1:length(trial_case),function(x) c(rep(1,times=trial_case[x]),rep(0,times=non_cases[x]))))
      #icc <- iccbin(cid,y,data=data.frame(cid=factor(cid),y=y),method='aov',ci.type='aov')
    }
  }
  power <- sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
  VE_est <- mean(ve_est2,na.rm=T)
  VE_sd <- sd(ve_est2,na.rm=T)
  powerb <- VE_estb <- VE_sdb <- c()
  if(adaptation==''){
    powerb <- sum(pval_binary_mle<0.05,na.rm=T)/sum(!is.na(pval_binary_mle))
    VE_estb <- mean(ve_est,na.rm=T)
    VE_sdb <- sd(ve_est,na.rm=T)
  }
  return(list(power, VE_est, VE_sd,powerb, VE_estb, VE_sdb,vaccinated_count, infectious_count, vaccinated_countb, infectious_countb))
}
for(des in 1:nCombAdapt){
  cluster_flag <- trial_designs$cluster[des]
  direct_VE <- trial_designs$VE[des]
  adaptation <- trial_designs$adapt[des]
  trial_designs$vaccinated[des] <- trial_results[[des]][[7]]
  trial_designs$infectious[des] <- trial_results[[des]][[8]]
  if(adaptation==''){
    trial_designs$vaccinated[des+nComb] <- trial_results[[des]][[9]]
    trial_designs$infectious[des+nComb] <- trial_results[[des]][[10]]
  }
  trial_designs$power[des] <- trial_results[[des]][[1]]
  trial_designs$VE_est[des] <- trial_results[[des]][[2]]
  trial_designs$VE_sd[des] <- trial_results[[des]][[3]]
  if(adaptation==''){
    trial_designs$power[des+nComb] <- trial_results[[des]][[4]]
    trial_designs$VE_est[des+nComb] <- trial_results[[des]][[5]]
    trial_designs$VE_sd[des+nComb] <- trial_results[[des]][[6]]
  }
}
subset(trial_designs,VE==0)
subset(trial_designs,VE>0)

result_table <- subset(trial_designs,VE>0)[,c(2,3,4,5,6,8,7,9)]
result_table$t1e <- subset(trial_designs,VE==0)$power
result_table$VE <- paste0(round(result_table$VE_est,2),' (',round(result_table$VE_sd,2),')')
result_table <- result_table[,-c(6:7)]
result_table$adapt <- as.character(result_table$adapt)
result_table$adapt[result_table$adapt==''] <- 'None'
result_table$cluster[result_table$cluster==0] <- 'Individual'
result_table$cluster[result_table$cluster==1] <- 'Cluster'
colnames(result_table) <- c('Randomisation','Adaptation','Weighting','Infectious','Vaccinated','Power','Type 1 error','VE estimate')
print(xtable(result_table), include.rownames = FALSE)


