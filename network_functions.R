## extracted from hitchings
infect_contacts <- function(potential_contacts,excluded_nodes,beta_value){
  contact_index <- susc_index <- rep(0,length(g_name))
  contact_index[potential_contacts] <- 1
  contact_index[excluded_nodes] <- 0
  susc_contacts <- g_name[contact_index] # potential_contacts[potential_contacts%in%node_class]#intersect(potential_contacts,node_class)
  num_neighbours_susc <- length(susc_contacts)
  # Sample from each group of neighbours in turn
  # First choose how many neighbours each node infects
  num_contacts_susc <- rbinom(1,num_neighbours_susc,1-exp(-beta_value))
  # Then sample from the neighbours
  # If one node gets picked twice by different nodes, just discard the duplicate.
  infectees_susc <- c()
  if(num_contacts_susc>0){
    sample_indices <- ceiling((runif(1) + 1:num_contacts_susc-1)*length(susc_contacts)/num_contacts_susc)
    infectees_temp <- susc_contacts[sample_indices] # sample(susc_contacts,num_contacts_susc)
    infectees_susc <- funique(infectees_temp)
  }
  return(infectees_susc)
}

spread <- function( s_nodes, v_nodes, e_nodes, i_nodes,c_nodes, beta, direct_VE,
                    incperiod_shape, incperiod_rate){
  # Spread will create new infected nodes from two sources: infectious nodes within the the study
  # population, and external pressure from the source population
  # Inputs:
  # s_nodes, e_nodes and i_nodes are susceptible, exposed and infected nodes
  # beta is the hazard of infection for one contact
  # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
  # length, currently drawn from a gamma distribution
  
  # Process: go through list of i_nodes, and choose a random number of its susceptible neighbours to
  # be infected, according to beta and choose a random number of its susceptible vaccinated neighbours to
  # be infected, according to beta and direct_VE
  infectees_susc <- c()
  infectees_vacc <- c()
  infectees_cont <- c()
  if (nrow(i_nodes)>0) {
    # Make a beta vector
    beta_v <- beta*(1-direct_VE)
    # Get a list of all neighbours of all infected nodes
    #potential_contacts <- unlist(ego(g,order=1,nodes=i_nodes[1,]))
    #neighbour_hh <- unlist(ego(g2,order=1,nodes=V(g)$hh[i_nodes[1,]]))
    #neighbours <- unlist(ego(neighbourhood_g,order=1,nodes=i_nodes[1,]))
    potential_contacts <- c()
    hr_contacts <- c()
    neighbours <- c()
    for(i in i_nodes[,1]) {
      potential_contacts <- c(potential_contacts,contact_list[[i]])
      hr_contacts <- c(hr_contacts,high_risk_list[[i]],household_list[[i]])
      neighbours <- c(neighbours,contact_of_contact_list[[i]])
    }
    excluded_nodes <- c(e_nodes[,1], i_nodes[,1], r_nodes)
    infectees_susc <- infect_contacts(potential_contacts,excluded_nodes,beta_value=beta)
    infectees_hr_susc <- infect_contacts(hr_contacts,excluded_nodes,beta_value=beta*(high_risk_scalar-1))
    infectees_n_susc <- infect_contacts(neighbours,excluded_nodes,beta_value=beta*neighbour_scalar)
    
    infectees_susc <- funique(c(infectees_susc,infectees_hr_susc,infectees_n_susc))
  } 
  
  newinfected_susc <- infectees_susc
  newinfected_vacc <- infectees_vacc
  newinfected_cont <- infectees_cont
  newinfected <- c(newinfected_susc, newinfected_vacc,newinfected_cont)
  newinfected <- funique(newinfected)
  
  if (length(newinfected)>0) {
    # Give each newly exposed node an incubation/latent period
    inc_periods <- rgamma(length(newinfected),shape=incperiod_shape,rate=incperiod_rate)
    # Add them to e_nodes and remove from s_nodes and v_nodes
    for(i in 1:length(newinfected))
      e_nodes <- rbind(e_nodes,c(newinfected[i],0,inc_periods[i]))
    s_nodes <- s_nodes[!s_nodes%in%newinfected_susc] # setdiff(s_nodes,newinfected_susc)
    v_nodes <- v_nodes[!v_nodes%in%newinfected_susc]
    c_nodes <- c_nodes[!c_nodes%in%newinfected_susc]
  }
  list(s_nodes, v_nodes, e_nodes,c_nodes)
}

recover <- function(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate) {
  # Input is a list of the exposed nodes, 
  # with number of days since infection and total incubation/latent
  # period, and equivalently for the infectious nodes.
  # For each of these nodes, we will add it to newinfectious if the number of days exposed has
  # reached the total length of the incubation period, and equivalently for the infectious nodes.
  
  # Advance infectious nodes first, otherwise you will doubly advance any nodes switching from
  # exposed to infectious at this time step
  indices_to_remove <- i_nodes[,2]>=i_nodes[,3]
  newremoved <- as.vector(i_nodes[indices_to_remove,1])
  
  # Add one day to the length of each infected individual's time infected
  i_nodes[,2] <- i_nodes[,2]+1
  
  # Remove any recovered from i_nodes and add to r_nodes
  i_nodes <- i_nodes[!(i_nodes[,1] %in% newremoved),,drop=FALSE]
  r_nodes <- c(r_nodes,newremoved)
  
  # Now advance exposed nodes
  indices_to_remove <- e_nodes[,2]>=e_nodes[,3]
  newinfectious <- as.vector(e_nodes[indices_to_remove,1])
  incubation_days <- as.vector(e_nodes[indices_to_remove,2])+1
  
  # Add one day to the length of each infected individual's time infected
  e_nodes[,2] <- e_nodes[,2]+1
  
  # Remove any progressing from e_nodes and add to i_nodes
  e_nodes <- e_nodes[!indices_to_remove,,drop=FALSE]
  if(length(newinfectious)>0){
    inf_periods <- rgamma(length(newinfectious),infperiod_shape,rate=infperiod_rate)
    #hosp_time <- rgamma(length(newinfectious),shape=hosp_shape,rate=hosp_rate)
    hosp_time <- rtruncnorm(length(newinfectious),a=0,mean=hosp_mean,sd=hosp_sd)
    min_time <- inf_periods
    for(i in 1:length(min_time)) if(hosp_time[i] < min_time[i]) min_time[i] <- hosp_time[i] 
    #min_time <- pmin(inf_periods,hosp_time)
    if(t<=recruitment_time) for(i in 1:length(min_time)) if(recruitment_time-t < min_time[i]) min_time[i] <- recruitment_time-t # min_time <- pmin(min_time,recruitment_time-t)
    #i_nodes <- rbind(i_nodes,cbind(newinfectious,rep(0,length(newinfectious)),min_time,incubation_days))
    for(i in 1:length(newinfectious)) i_nodes <- rbind(i_nodes,c(newinfectious[i],0,min_time[i],incubation_days[i]))
  }
  
  list(e_nodes, i_nodes, r_nodes, newinfectious)
}

simulate_contact_network <- function(beta,neighbour_scalar,high_risk_scalar,first_infected,inf_time,end_time=40){
  trajectories <- list()
  trajectories$S <- length(vertices) - 1
  trajectories$E <- 0
  trajectories$I <- 0
  trajectories$R <- 0
  
  e_nodes <- matrix(nrow=0,ncol=3)
  i_nodes <- matrix(nrow=0,ncol=4)
  v_nodes <- c()
  c_nodes <- c()
  s_nodes <- g_name
  r_nodes <- c()
  
  inc_time <- rgamma(length(first_infected),shape=incperiod_shape,rate=incperiod_rate)
  i_nodes <- rbind(i_nodes,c(first_infected,rep(0,length(first_infected)),inf_time,inc_time))
  s_nodes <- setdiff(s_nodes,first_infected)
  
  #recruitment_time <- round(rgamma(1,shape=recruit_shape,rate=recruit_rate))
  recruitment_time <- ceiling(rtruncnorm(1,a=0,mean=recruit_mean,sd=recruit_sd))
  results <- matrix(c(first_infected,0,recruitment_time,-inc_time),nrow=1)
  numinfectious <- 1
  ##!! add in additional infectious people?
  
  contacts <- contact_list[[first_infected]]
  contacts <- contacts[contacts!=first_infected]
  ## identify high-risk people
  ##!! all household members are high risk. 
  ##!! does it mean anything for a contact_of_contact (i.e. neighbour) to be high risk?
  high_risk <- high_risk_list[[first_infected]]
  ## contacts of contacts
  contacts_of_contacts <- contact_of_contact_list[[first_infected]]
  contacts_of_contacts <- contacts_of_contacts[contacts_of_contacts!=first_infected]
  ## add households of high-risk contacts to contacts of contacts
  if(length(high_risk)>0) 
    for(hr in high_risk)
      contacts_of_contacts <- c(contacts_of_contacts,household_list[[hr]])
  high_risk <- c(high_risk,household_list[[first_infected]])
  cluster_people <- funique(c(contacts,contacts_of_contacts))
  
  sim_time <- recruitment_time + end_time
  for(t in 1:sim_time){
    
    newinfectious <- c()
    if ((nrow(i_nodes)>0)||(nrow(e_nodes)>0)) {
      rec_list <- recover(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate)
      e_nodes <- rec_list[[1]]
      i_nodes <- rec_list[[2]]
      r_nodes <- rec_list[[3]]
      newinfectious <- rec_list[[4]]
    }
    
    #sub_g <- if(ncol(i_nodes)>0) induced_subgraph(g,c(i_nodes[1,],s_nodes,v_nodes)) else NULL
    spread_list <- spread(s_nodes,v_nodes,e_nodes,i_nodes,c_nodes,beta,direct_VE,incperiod_shape,incperiod_rate)
    s_nodes <- spread_list[[1]]
    v_nodes <- spread_list[[2]]
    e_nodes <- spread_list[[3]]
    c_nodes <- spread_list[[4]]
    
    numnewinfectious <- length(newinfectious)
    if (numnewinfectious>0) {
      # Update results
      results <- rbind(results,cbind(newinfectious,t,recruitment_time,t-as.numeric(i_nodes[match(newinfectious,i_nodes[,1]),4])))
      
      numinfectious <- numinfectious+numnewinfectious
    }
    
    trajectories$S[t+1] <- length(s_nodes) + length(v_nodes) + length(c_nodes)
    trajectories$E[t+1] <- trajectories$S[t] - trajectories$S[t+1]
    trajectories$I[t+1] <- numnewinfectious
    trajectories$R[t+1] <- length(r_nodes)-sum(trajectories$R)
    
  }
  
  
  results <- as.data.frame(results)
  colnames(results) <- c('InfectedNode', 'DayInfectious', 'RecruitmentDay', 'DayInfected')
  results$inCluster <- results$InfectedNode%in%cluster_people
  results$contact <- results$InfectedNode%in%contacts
  results$highrisk <- results$InfectedNode%in%high_risk
  
  return(list(results,length(cluster_people),recruitment_time))
}