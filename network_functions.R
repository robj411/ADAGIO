## extracted from hitchings
infect_contacts <- function(potential_contacts,beta_value){
  num_neighbours_susc <- length(potential_contacts)
  # Sample from each group of neighbours in turn
  # First choose how many neighbours each node infects
  num_contacts_inf <- rbinom(1,num_neighbours_susc,1-exp(-beta_value))
  # Then sample from the neighbours
  # If one node gets picked twice by different nodes, just discard the duplicate.
  infectees_susc <- c()
  if(num_contacts_inf>0){
    sample_indices <- ceiling((runif(1) + 1:num_contacts_inf-1)*num_neighbours_susc/num_contacts_inf)
    infectees_temp <- potential_contacts[sample_indices] 
    infectees_susc <- funique(infectees_temp)
  }
  return(infectees_susc)
}

infect_from_source <- function( s_nodes, v_nodes, e_nodes_info, direct_VE,incperiod_shape, incperiod_rate,rate_from_source){
  # identify susceptibles
  infectees_susc <- c()
  s_hr_l <- s_nodes==1
  ##!! not trial ppts for now
  s_hr <- g_name[s_hr_l & v_nodes==0 & t_nodes==0]
  # infect
  if(length(s_hr)>0)
    infectees_susc <- funique(infect_contacts(s_hr,beta_value=rate_from_source))
  newinfected <- infectees_susc
  # infect vaccinated
  if(sum(v_nodes)>0){
    infectees_susc <- c()
    beta_v <- rate_from_source*(1-direct_VE)
    s_hr <- g_name[s_hr_l & v_nodes==1 & t_nodes==0]
    if(length(s_hr)>0)
      infectees_susc <- funique(infect_contacts(s_hr,beta_value=beta_v))
    newinfected <- c(newinfected,infectees_susc)
  }
  
  if (length(newinfected)>0) {
    # Give each newly exposed node an incubation/latent period
    inc_periods <- incperiod_const + rgamma(length(newinfected),shape=incperiod_shape,rate=incperiod_rate)
    # Add them to e_nodes and remove from s_nodes and v_nodes
    for(i in 1:length(newinfected))
      e_nodes_info <- rbind(e_nodes_info,c(newinfected[i],0,inc_periods[i]))
  }
  return(e_nodes_info)
}

spread <- function( s_nodes, v_nodes, e_nodes_info, current_infectious, direct_VE,incperiod_shape, incperiod_rate,susc_list=contact_list,beta_scalar=1){
  # Spread will create new infected nodes from two sources: infectious nodes within the the study
  # population, and external pressure from the source population
  # Inputs:
  # s_nodes, e_nodes and i_nodes are susceptible, exposed and infected nodes
  # beta_base is the hazard of infection for one contact
  # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
  # length, currently drawn from a gamma distribution
  scaled_beta <- beta_scalar * beta_base
  # Process: go through list of i_nodes, and choose a random number of its susceptible neighbours to
  # be infected, according to beta_base and choose a random number of its susceptible vaccinated neighbours to
  # be infected, according to beta_base and direct_VE
  infectees_susc <- c()
  # Get a list of all neighbours of all infected nodes
  potential_contacts <- c()
  #potential_hr_contacts <- c()
  #potential_neighbours <- c()
  for(i in current_infectious) {
    #potential_hr_contacts <- c(potential_hr_contacts,high_risk_list[[i]],household_list[[i]])
    #potential_neighbours <- c(potential_neighbours,contact_of_contact_list[[i]])
    potential_contacts <- c(potential_contacts,susc_list[[i]])
  }
  # remove duplications
  #potential_contacts <- potential_contacts[!potential_contacts%in%potential_hr_contacts]
  # infect high risk
  #s_hr_l <- s_nodes[potential_hr_contacts]==1
  #v_hr_l <- v_nodes[potential_hr_contacts]
  #s_hr <- potential_hr_contacts[s_hr_l & v_hr_l==0]
  #if(length(s_hr)>0)
  #  infectees_hr_susc <- funique(infect_contacts(s_hr,beta_value=scaled_beta*high_risk_scalar))
  # infect neighbours
  #s_nb_l <-  s_nodes[potential_neighbours]==1
  #v_nb_l <-  v_nodes[potential_neighbours]
  #s_nb <- potential_neighbours[s_nb_l & v_nb_l==0]
  #if(length(s_nb)>0)
  #  infectees_n_susc <- funique(infect_contacts(s_nb,beta_value=scaled_beta*neighbour_scalar))
  # infect other contacts
  if(length(potential_contacts)>0){
    s_contacts_l <- s_nodes[potential_contacts]==1
    v_contacts_l <- v_nodes[potential_contacts]
    s_contacts <- potential_contacts[s_contacts_l & v_contacts_l==0]
    if(length(s_contacts)>0)
      infectees_susc <- funique(infect_contacts(s_contacts,beta_value=scaled_beta))
  }
  newinfected <- funique(infectees_susc) #funique(c(infectees_susc,infectees_hr_susc,infectees_n_susc))
  # infect vaccinated
  if(sum(v_nodes)>0){
    infectees_susc <- c()
    beta_v <- scaled_beta*(1-direct_VE)
    # infect high risk
    #s_hr <- potential_hr_contacts[s_hr_l & v_hr_l==1]
    #if(length(s_hr)>0)
    #  infectees_hr_susc <- funique(infect_contacts(s_hr,beta_value=beta_v*high_risk_scalar))
    # infect neighbours
    #s_nb <- potential_neighbours[s_nb_l & v_nb_l==1]
    #if(length(s_nb)>0)
    #  infectees_n_susc <- funique(infect_contacts(s_nb,beta_value=beta_v*neighbour_scalar))
    # infect other contacts
    if(length(potential_contacts)>0){
      s_contacts <- potential_contacts[s_contacts_l & v_contacts_l==1]
      if(length(s_contacts)>0)
        infectees_susc <- funique(infect_contacts(s_contacts,beta_value=beta_v))
    }
    newinfected <- c(newinfected,funique(infectees_susc))
  }
  
  if (length(newinfected)>0) {
    # Give each newly exposed node an incubation/latent period
    inc_periods <- incperiod_const + rgamma(length(newinfected),shape=incperiod_shape,rate=incperiod_rate)
    # Add them to e_nodes and remove from s_nodes and v_nodes
    for(i in 1:length(newinfected))
      e_nodes_info <- rbind(e_nodes_info,c(newinfected[i],0,inc_periods[i]))
  }
  return(e_nodes_info)
}

recover <- function(e_nodes_info,i_nodes_info,infperiod_shape,infperiod_rate,cluster_people_index,time_diff=NULL) {
  # Input is a list of the exposed nodes, 
  # with number of days since infection and total incubation/latent
  # period, and equivalently for the infectious nodes.
  # For each of these nodes, we will add it to newinfectious if the number of days exposed has
  # reached the total length of the incubation period, and equivalently for the infectious nodes.
  
  # Advance infectious nodes first, otherwise you will doubly advance any nodes switching from
  # exposed to infectious at this time step
  indices_to_remove <- i_nodes_info[,2]>=i_nodes_info[,3]
  newremoved <- as.vector(i_nodes_info[indices_to_remove,1])
  
  # Add one day to the length of each infected individual's time infected
  i_nodes_info[,2] <- i_nodes_info[,2]+1
  
  # Remove any recovered from i_nodes and add to r_nodes
  i_nodes_info <- i_nodes_info[!indices_to_remove,,drop=FALSE]
  
  # Now advance exposed nodes
  indices_to_remove <- e_nodes_info[,2]>=e_nodes_info[,3]
  newinfectious <- as.vector(e_nodes_info[indices_to_remove,1])
  incubation_days <- as.vector(e_nodes_info[indices_to_remove,2])+1
  
  # Add one day to the length of each infected individual's time infected
  e_nodes_info[,2] <- e_nodes_info[,2]+1
  
  # Remove any progressing from e_nodes_info and add to i_nodes
  e_nodes_info <- e_nodes_info[!indices_to_remove,,drop=FALSE]
  if(length(newinfectious)>0){
    inf_periods <- infperiod_const + rgamma(length(newinfectious),infperiod_shape,rate=infperiod_rate)
    #hosp_time <- rgamma(length(newinfectious),shape=hosp_shape,rate=hosp_rate)
    if(!is.null(time_diff)){
      hosp_time <- c()
      for(i in 1:length(newinfectious))
        if(cluster_people_index[newinfectious[i]]==1){
          hosp_time[i] <- rtruncnorm(1,a=0,mean=hosp_mean,sd=hosp_sd)
        }else{
          hosp_time[i] <- rtruncnorm(1,a=0,mean=hosp_mean_index,sd=hosp_sd_index)
        }
      # person is infectious until the minimal time that they are hospitalised or removed otherwise
      for(i in 1:length(inf_periods)) if(hosp_time[i] < inf_periods[i]) inf_periods[i] <- hosp_time[i] 
      #min_time <- pmin(inf_periods,hosp_time)
      # if person is enrolled soon, they will be hospitalised then
      if(0<=time_diff) for(i in 1:length(inf_periods)) if(time_diff < inf_periods[i]) inf_periods[i] <- time_diff 
    }
    #i_nodes <- rbind(i_nodes,cbind(newinfectious,rep(0,length(newinfectious)),min_time,incubation_days))
    for(i in 1:length(newinfectious)) 
      i_nodes_info <- rbind(i_nodes_info,c(newinfectious[i],0,inf_periods[i],incubation_days[i],runif(1)<observed))
  }
  
  list(e_nodes_info, i_nodes_info, newremoved, newinfectious)
}

ebola_spread_wrapper <- function(i_nodes_info,s_nodes,v_nodes,e_nodes_info,direct_VE){
  # to contacts
  current_infectious <- i_nodes_info[,1]
  if(length(current_infectious)>0){
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_list,beta_scalar=1)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=hr_and_hh_list,beta_scalar=high_risk_scalar-1)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes_info <- spread(s_nodes,v_nodes,e_nodes_info,current_infectious,direct_VE,incperiod_shape,incperiod_rate,susc_list=contact_of_contact_list,beta_scalar=neighbour_scalar)
  }
  return(e_nodes_info)
}

simulate_contact_network <- function(first_infected,individual_recruitment_times=F,end_time=31,start_day=0,from_source=0,cluster_flag=0,allocation_ratio=0.5,
                                     direct_VE=0,base_rate=0,spread_wrapper=ebola_spread_wrapper){
  # set up info to store
  trajectories <- list()
  trajectories$S <- length(vertices) - 1
  trajectories$E <- 0
  trajectories$I <- 0
  trajectories$R <- 0
  e_nodes_info <- matrix(nrow=0,ncol=3)
  i_nodes_info <- matrix(nrow=0,ncol=5)
  e_nodes <- rep(0,length(g_name))
  i_nodes <- rep(0,length(g_name))
  v_nodes <- rep(0,length(g_name))
  c_nodes <- rep(0,length(g_name))
  s_nodes <- rep(1,length(g_name))
  r_nodes <- rep(0,length(g_name))
  t_nodes <- rep(0,length(g_name))
  
  # generate info for index case
  inc_time <- incperiod_const + rgamma(length(first_infected),shape=incperiod_shape,rate=incperiod_rate)
  ceil_inc_time <- ceiling(inc_time)
  #i_nodes_info <- rbind(i_nodes_info,c(first_infected,rep(0,length(first_infected)),inf_time,inc_time))
  e_nodes_info <- rbind(e_nodes_info,c(first_infected,0,inc_time))
  s_nodes[first_infected] <- 0
  e_nodes[first_infected] <- 1
  
  #recruitment_time <- round(rgamma(1,shape=recruit_shape,rate=recruit_rate))
  recruitment_time <- ceiling(rtruncnorm(1,a=0,mean=recruit_mean,sd=recruit_sd))
  results <- matrix(nrow=0,ncol=5)#c(first_infected,0,-inc_time,NA),nrow=1)
  numinfectious <- 1
  ##!! add in additional infectious people?
  
  # identify contacts of index case
  contacts <- contact_list[[first_infected]]
  contacts <- contacts[contacts!=first_infected]
  order_infected <- first_infected
  ## identify high-risk people
  ##!! all household members are high risk. 
  if(!exists('high_risk_list')) high_risk_list <- lapply(g_name,function(x)c())
  high_risk <- high_risk_list[[first_infected]]
  ## contacts of contacts
  contacts_of_contacts <- contact_of_contact_list[[first_infected]]
  contacts_of_contacts <- contacts_of_contacts[contacts_of_contacts!=first_infected]
  ## add households of high-risk contacts to contacts of contacts
  if(length(high_risk)>0) 
    for(hr in high_risk)
      contacts_of_contacts <- c(contacts_of_contacts,household_list[[hr]])
  #high_risk <- c(high_risk,household_list[[first_infected]])
  cluster_people <- funique(c(contacts,contacts_of_contacts))
  cluster_people_index <- g_name%in%cluster_people
  
  # enroll trial participants
  n_trial_participants <- rbinom(1,length(cluster_people),enrollment_rate)
  trial_participants <- sample(cluster_people,n_trial_participants,replace=F)
  t_nodes[trial_participants] <- 1
  t_nodes <<- t_nodes
  vaccinees <- c()
  if(cluster_flag==0){
    nvacc <- round(length(trial_participants)*allocation_ratio)
    vaccinees <- trial_participants[sort(sample.int(length(trial_participants),nvacc,replace=F))]
  }else{
    if(runif(1)<allocation_ratio)
      vaccinees <- trial_participants
  }
  vaccine_incubation_times <- 0
  if(length(vaccinees)>0)
    vaccine_incubation_times <- rgamma(length(vaccinees),shape=vacc_shape,rate=vacc_rate)
  if(individual_recruitment_times==F){
    recruitment_times <- rep(recruitment_time,n_trial_participants) + ceil_inc_time
  }else{
    recruitment_times <- sample(1:recruitment_time,n_trial_participants,replace=T) + ceil_inc_time
  }
  # roll epidemic forward one day at a time
  sim_time <- recruitment_time + end_time + ceil_inc_time
  for(time_step in 1:sim_time){
    ## vaccination given time to develop immunity
    if(length(vaccinees)>0) {
      developed <- vaccine_incubation_times<=time_step-recruitment_times[trial_participants%in%vaccinees]
      v_nodes[vaccinees[developed]] <- 1
    }
    
    # update everyone's internal clock
    newinfectious <- newremoved <- c()
    if ((nrow(e_nodes_info)>0)||(nrow(i_nodes_info)>0)) {
      time_diff <- NULL
      if(!individual_recruitment_times) time_diff <- recruitment_time-time_step
      rec_list <- recover(e_nodes_info,i_nodes_info,infperiod_shape,infperiod_rate,cluster_people_index=cluster_people_index,time_diff=time_diff)
      e_nodes_info <- rec_list[[1]]
      i_nodes_info <- rec_list[[2]]
      newremoved <- rec_list[[3]]
      i_nodes[newremoved] <- 0
      r_nodes[newremoved] <- 1 
      newinfectious <- rec_list[[4]]
      e_nodes[newinfectious] <- 0
      i_nodes[newinfectious] <- 1
    }
    # store new cases of infectiousness
    numnewinfectious <- length(newinfectious)
    if (numnewinfectious>0) {
      # Update results
      ord <- match(newinfectious,i_nodes_info[,1])
      results <- rbind(results,cbind(newinfectious,time_step,time_step-as.numeric(i_nodes_info[ord,4]),NA,i_nodes_info[ord,5]))
      numinfectious <- numinfectious+numnewinfectious
    }
    if(length(newremoved)>0){
      results[results[,1]%in%newremoved,4] <- time_step
    }
    
    ## spread infection
    e_nodes_info <- spread_wrapper(i_nodes_info,s_nodes,v_nodes,e_nodes_info,direct_VE)
    s_nodes[e_nodes_info[,1]] <- 0
    e_nodes[e_nodes_info[,1]] <- 1
    order_infected <- c(order_infected,e_nodes_info[,1])
    
    # infect from source
    rate_from_source <- max((start_day + time_step)*from_source + base_rate, 0)
    if(rate_from_source>0){
      e_nodes_info <- infect_from_source(s_nodes,v_nodes,e_nodes_info,direct_VE,incperiod_shape,incperiod_rate,rate_from_source)
      s_nodes[e_nodes_info[,1]] <- 0
      e_nodes[e_nodes_info[,1]] <- 1
      order_infected <- c(order_infected,e_nodes_info[,1])
    }
    
    # store information
    trajectories$S[time_step+1] <- sum(s_nodes)# + sum(v_nodes) + sum(c_nodes)
    trajectories$E[time_step+1] <- trajectories$S[time_step] - trajectories$S[time_step+1]
    trajectories$I[time_step+1] <- numnewinfectious
    trajectories$R[time_step+1] <- sum(r_nodes)-sum(trajectories$R)
    
  }
  
  # store information and format to return
  results <- as.data.frame(results)
  colnames(results) <- c('InfectedNode', 'DayInfectious', 'DayInfected', 'DayRemoved','Observed')
  results$inCluster <- results$InfectedNode%in%cluster_people
  results$contact <- results$InfectedNode%in%contacts
  #results$highrisk <- results$InfectedNode%in%high_risk
  results$inTrial <- results$InfectedNode%in%trial_participants
  results$vaccinated <- results$InfectedNode%in%vaccinees
  results$RecruitmentDay <- recruitment_times[match(results$InfectedNode,trial_participants)]
  
  return(list(results,length(cluster_people),recruitment_times,length(vaccinees),length(trial_participants),vaccinees,trial_participants,order_infected,vaccine_incubation_times+recruitment_times[trial_participants%in%vaccinees]))
}
