
## from hitchings
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

## script as function
core_trial_script <- function(trial_design,g){
  trial_startday <- trial_design$trial_startday#100#50 + (sday-1)*100
  trial_length <- trial_design$trial_length#500 - trial_startday
  bTrial <- trial_design$bTrial
  bCluster <- trial_design$bCluster
  #adaptation <- trial_designs[[tr]]$adaptation
  #vaccination_gap <- trial_designs[[tr]]$vaccination_gap
  follow_up <- trial_design$follow_up
  #adaptation_day <- trial_designs[[tr]]$adaptation_day
  #profvis(
  list[results,trial_nodes,trajectories,allocation_rates]<-
    network_epidemic(g,disease_dynamics,direct_VE,infected_trajectory,trial_design)
  #)
  
  list[VE,pval,events_vacc,events_cont,analysed_trialsize] <- 
    analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,bCluster,follow_up,trial_design$revisit)
  VE <- VE[1]
  # if(tr!=2){
  ## add analysis for ring-end
  if(trial_design$reevaluate==1){
    pvals$DiRCT[simnum] <- pval
    VaccineEfficacy$DiRCT[simnum] <- VE[1]
    # duplicate results with different end point
    list[VE2,pval2,events_vacc,events_cont,analysed_trialsize] <- 
      analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,bCluster,follow_up,revisit=1)
    pval <- c(pval2,pval)
    VE <- c(VE2[1],VE[1])
  }
  # }else{
  # list[VE_gaussian_coxme,pval_gaussian_coxme,VE_gaussian_coxph,pval_gaussian_coxph,
  #      VE_gamma_coxph,pval_gamma_coxph,VE_gee,pval_gee,events_vacc,events_cont,analysed_trialsize,
  #      ICC,deff,prop_zeros] <- analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,bCluster,vaccination_gap)
  # 
  # pvals$cRCT_gaussian_coxme[simnum] <- pval_gaussian_coxme
  # VaccineEfficacy$cRCT_gaussian_coxme[simnum] <- VE_gaussian_coxme[1]
  # 
  # pvals$cRCT_gaussian_coxph[simnum] <- pval_gaussian_coxph
  # VaccineEfficacy$cRCT_gaussian_coxph[simnum] <- VE_gaussian_coxph[1]
  # 
  # pvals$cRCT_gamma_coxph[simnum] <- pval_gamma_coxph
  # VaccineEfficacy$cRCT_gamma_coxph[simnum] <- VE_gamma_coxph[1]
  # 
  # pvals$cRCT_gee[simnum] <- pval_gee
  # VaccineEfficacy$cRCT_gee[simnum] <- VE_gee[1]
  # 
  # ICCs[simnum] <- ICC
  # deffs[simnum] <- deff
  # props_zeros[simnum] <- prop_zeros
  # 
  # pval <- pval_gee
  # VE <- VE_gee[1]
  #}
  num_enrolled <- nrow(trial_nodes)#analysed_trialsize
  num_vacc <- sum(trial_nodes$TrialStatus==1)
  numevents <- nrow(results)
  list(num_enrolled=num_enrolled,num_vacc=num_vacc,events_vacc=events_vacc,events_cont=events_cont,numevents=numevents,pval=pval,VaccineEfficacy=VE,
       trajectories=trajectories,vaccinationDays=trial_nodes$DayVaccinated,allocation_rate=allocation_rates)
}

## script as function
outer_trial_script <- function(nsim,direct_VE,trial_indicies,trial_designs){
  trials <- length(trial_designs)
  extra_trials <- sum(sapply(trial_designs,function(x) x$reevaluate))
  allocation_rates <- numevents_cont <- numevents <- numevents_vacc <- num_vacc <-  num_enrolled <- matrix(NA,nrow=trials+extra_trials,ncol=nsim)
  allocation_rate_list <- trajectory_list <- list()
  for (simnum in 1:nsim) {
    if(direct_VE==0.6) trajectory_list[[simnum]] <- list() 
    g <<- make_network(ave_community_size, community_size_range, num_communities,rate_within, rate_between)
    
    trial_outcomes <- foreach(tr = trial_indicies) %dopar% {
      core_trial_script(trial_design=trial_designs[[tr]],g)
    }
    allocation_rate_list[[simnum]] <- sapply(trial_outcomes,function(x)x$allocation_rate)
    index <- 0
    for(tr in trial_indicies){
      
      ## add analysis for ring-end
      index <- index + 1
      if(trial_designs[[tr]]$reevaluate==1) index <- c(index,index+1)
      allocation_rates[index,simnum] <- trial_outcomes[[tr]]$allocation_rate[length(trial_outcomes[[tr]]$allocation_rate)]
      num_enrolled[index,simnum] <- trial_outcomes[[tr]]$num_enrolled
      num_vacc[index,simnum] <- trial_outcomes[[tr]]$num_vacc
      numevents_vacc[index,simnum] <- trial_outcomes[[tr]]$events_vacc
      numevents_cont[index,simnum] <- trial_outcomes[[tr]]$events_cont
      numevents[index,simnum] <- trial_outcomes[[tr]]$numevents
      mle_pvals[simnum,index] <- trial_outcomes[[tr]]$pval
      mle_VaccineEfficacy[simnum,index] <- trial_outcomes[[tr]]$VaccineEfficacy
      if(direct_VE==0.6) trajectory_list[[simnum]][[tr]] <- trial_outcomes[[tr]]$trajectories
      index <- max(index)
    }
    cat("Simulation ",simnum,"\n")
  }
  return(list(allocation_rates,num_enrolled,num_vacc,numevents_vacc,numevents_cont,numevents,mle_pvals,mle_VaccineEfficacy,trajectory_list,allocation_rate_list))
}

## adapted from hitchings
make_network <- function(ave_community_size, community_size_range, 
                         num_communities, rate_within, rate_between) {
  # Function to make a network of closely-connected communities that are more sparsely
  # connected with each other. I use a stochastic block model.
  # Inputs:
  # Ave_community_size is the average number of people in a community
  # Community_size_range is the range of community sizes. Currently size of a community
  # is being chosen from a uniform distribution on (ave-range/2, ave+range/2)
  # Num_communities is the number of communities in the study population
  # rate_within is the probability of an edge between any two nodes within the same community
  # rate_between is the probability of an edge between any two nodes in different communities
  
  # Create the network, and assign all members a community number
  # The community number will be output by the epidemic function and used in the 
  # gamma frailty model and calculation of the ICC
  community_sizes <- ave_community_size + round(runif(num_communities,-community_size_range/2,community_size_range/2))
  studypop_size <- sum(community_sizes)
  # All communities have the same connectedness, and all communities are equally
  # connected to each other
  within_rates <- diag(nrow=num_communities,ncol=num_communities,x=rate_within)
  between_rates <- matrix(rate_between,nrow=num_communities,ncol=num_communities) -
    diag(nrow=num_communities,ncol=num_communities,x=rate_between)
  rates <- within_rates+between_rates
  g <- sample_sbm(studypop_size,rates,community_sizes)
  # Give the nodes a name so that igraph remembers them
  V(g)$name <- 1:studypop_size
  V(g)$community <- rep(1:num_communities,community_sizes)
  # Trial status will track whether a node is not in the trial (NA), in the control arm (0) or
  # in the vaccine arm (1)
  V(g)$trialstatus <- NA
  V(g)$enrollmentday <- NA
  V(g)$treatmentday <- NA
  
  return(g)
  
}
## adapted from hitchings
network_epidemic<-function(g,disease_dynamics,direct_VE,infected_trajectory,trial_design) {
  # Inputs:
  # g - the graph to run the epidemic on
  # beta - Every infectious individual contacts all their neighbours in a time step
  # and infects each susceptible with hazard beta. So beta represents the hazard of infection from one
  # contact between an infectious and a susceptible.
  # num_introductions - how many separate introductions we expect on average from the main epidemic. This is used
  # to calibrate the external force of infection
  # direct_VE - direct leaky efficacy of the vaccine
  # bTrial - whether we are recruiting ad hoc (1) or with a ring (2)
  # bCluster - indicator of whether we are running the cRCT (1) or the iRCT (0)
  # trial_startday - first day of trial enrollment
  # trial_length - end of follow-up of trial partcipants, counting from the first day of enrollment
  # num_enrolled_per_day - number of individuals/clusters enrolled per day
  # enrollment_period - length of enrollment period
  # cluster_coverage - The proportion of each cluster we expect to enroll
  {
    for(i in 1:length(trial_design)) assign(names(trial_design)[i],trial_design[[names(trial_design)[i]]])
    for(i in 1:length(disease_dynamics)) assign(names(disease_dynamics)[i],disease_dynamics[[names(disease_dynamics)[i]]])
    
    trajectories <- list()
    trajectories$S <- c()
    trajectories$E <- c()
    trajectories$I <- c()
    trajectories$R <- c()
    num_timesteps <- trial_startday + trial_length + enrollment_period - 1
    times <- seq(0,num_timesteps,1)
  # Define how the study population is linked to the source population
  # Connect all individuals to source population at same hazard
  # Constant of proportionality varies by community
  connected_to_source <- V(g)$name
  g_name <- V(g)$name
  g_community <- V(g)$community
  
  # Calibrate extF to the number of introductions, given the progression of the epidemic in the source population
  num_communities <- max(V(g)$community)
  #nodes_by_community <- lapply(1:num_communities,function(x) V(g)[V(g)$community==x])
  comm_sizes <- sapply(1:num_communities,function(x) sum(g_community==x))
  sumsqrt <- sum(sqrt(comm_sizes))
  # each community has a constant source of spontaneous infection
  extF <- -log(1-num_introductions/(sqrt(comm_sizes)*sumsqrt))/trapz(times,infected_trajectory)
  
  # Parameters to do with trial recruitment
  # Enrollment per day is number of clusters enrolled per day
  enrollment_schedule <- rep(num_enrolled_per_day,enrollment_period)
  trial_days <- seq(trial_startday,trial_startday+enrollment_period*enrollment_gap-1,enrollment_gap)
  vaccination_days <- seq(trial_startday,trial_startday+trial_length-1,vaccination_gap)
  non_trial_clusters <- 1:max(g_community)
  
  if (adaptation_day>0) {
    # Parameters to do with trial recruitment
    # Enrollment per day is number of clusters enrolled per day
    adaptation_days <- seq(trial_startday,trial_startday+trial_length-1,adaptation_day)
    non_trial_clusters <- 1:max(g_community)
  }
  number_to_treat <- cluster_coverage*sum(comm_sizes)/length(vaccination_days)
  
  # Initialize the S, E, I, and R nodes. I seed the epidemic from an SIR curve in a source population,
  # so initially all nodes in the study population are susceptible
  # i_nodes and e_nodes are matrices. The first row is the  identity of the node. The second row
  # is the number of days since infection/infectiousness. The third row is the total incubation/infectious period, 
  # drawn from a distribution when it becomes infected/infectious.
  # We are only going to consider a vaccine with effects on susceptibility, so only need one
  # vaccinated class
  e_nodes <- matrix(nrow=3,ncol=0)
  i_nodes <- matrix(nrow=3,ncol=0)
  v_nodes <- c()
  c_nodes <- c()
  s_nodes <- as.vector(V(g))
  r_nodes <- c()
  vertices <- V(g)
  
  # Initialize results.
  # Results will be, for each newly-infected node, the identity of the node, the day it was infected,
  # the community of which it is a member and its trial status at the time of infection. 
  # This is enough information to run a Cox PH with gamma frailty.
  # Make a data frame the size of the study pop and fill it in, then trim at the end
  results <- data.frame("InfectedNode"=numeric(),
                        "DayInfected"=numeric(),
                        "Community"=numeric(),
                        "TrialStatus"=numeric(),
                        "DayEnrolled"=numeric())
  numinfectious <- 0
  allocation_rate <- 0.5
  
  trial_nodes_info <- data.frame("Node"=numeric(),
                                 "Community"=numeric(),
                                 "TrialStatus"=numeric(),
                                 "DayEnrolled"=numeric(),
                                 "DayVaccinated"=numeric())
  allocation_rates <- rep(allocation_rate,num_timesteps)
  count=0
  }
  for (t in 1:num_timesteps) {
    
    
    if(adaptation!='' && 
       bCluster==0 && 
       sum(!is.na(results$TrialStatus))>0 && #&(results$DayInfected-results$DayVaccinated)>ave_inc_period
       any(results$DayInfected-results$DayVaccinated>ave_inc_period,na.rm=T) &&
       t%in%adaptation_days){
      ## subset those (a) in trial and (b) enrolled at least follow_up days ago
      recently_vaccinated <- !is.na(results$DayVaccinated) #& results$DayVaccinated < t - follow_up
      eligible_results <- results[recently_vaccinated,]
      recently_vaccinated <- !is.na(trial_nodes_info$DayVaccinated) #& trial_nodes_info$DayVaccinated < t - follow_up
      eligible_trial_nodes <- trial_nodes_info[recently_vaccinated,]
      ## select those to visit. Either all, or only those who haven't been visited before.
      results_to_visit <- nodes_to_visit <- T
      if(revisit==0){
        results_to_visit <- eligible_results$DayVaccinated > t - adaptation_day - follow_up
        nodes_to_visit <- eligible_trial_nodes$DayVaccinated > t - adaptation_day - follow_up
      }
      ## count successes and fails
      excluded0 <- sum(eligible_results$TrialStatus==0&eligible_results$DayInfected-eligible_results$DayVaccinated<=ave_inc_period&results_to_visit,na.rm=T)
      excluded1 <- sum(eligible_results$TrialStatus==1&eligible_results$DayInfected-eligible_results$DayVaccinated<=ave_inc_period&results_to_visit,na.rm=T)
      fail0 <- sum(eligible_results$TrialStatus==0&results_to_visit,na.rm=T) - excluded0
      fail1 <- sum(eligible_results$TrialStatus==1&results_to_visit,na.rm=T) - excluded1
      total0 <- sum(eligible_trial_nodes$TrialStatus==0&nodes_to_visit,na.rm=T) - excluded0
      total1 <- sum(eligible_trial_nodes$TrialStatus == 1&nodes_to_visit,na.rm=T) - excluded1
      success0 <- total0 - fail0
      success1 <- total1 - fail1
      if(adaptation%in%c('Ros','Ney')){
        p0 <- success0/total0
        p1 <- success1/total1
        if(adaptation=='Ros'){
          R_val <- sqrt(p1/p0)  # ros
          allocation_rate <- R_val / (1+R_val)
        }else if(adaptation=='Ney'){
          allocation_rate <- ifelse(p0*(1-p0)==0||p1*(1-p1)==0, 0.5, sqrt(p1*(1-p1)) / (sqrt(p0*(1-p0))+ sqrt(p1*(1-p1))) )# ney
          print(c(excluded0,total0,excluded1,total1,p0,p1,p0*(1-p0)==0,p1*(1-p1)==0,p0*(1-p0)==0||p1*(1-p1)==0,allocation_rate))
        }
      }else if(adaptation%in%c('TS','TST')){
        j <- t - trial_startday
        bigT <- trial_length
        tuning_c <- ifelse(adaptation=='TS',1,j/bigT)
        #print(tuning_c)
        p0 <- rbeta(1000,1+success0,1+fail0)
        p1 <- rbeta(1000,1+success1,1+fail1)
        prob1 <- sum(p1>p0)/1000
        allocation_rate <- prob1^tuning_c / (prob1^tuning_c + (1 - prob1)^tuning_c)
      }
      allocation_rates[t:num_timesteps] <- allocation_rate
    }
    
    # I'm recovering first, so I need to ensure that everyone has at least one chance to infect.
    # I do this by initializing an infectious node with 0 days since infection, seeing whether they
    # recover, then advancing them one day along their infectious period.
    if (bTrial==1) { ## fixed enrollment / cluster enrollment
      # Recruit and randomize if during the enrollment period
      if (t%in%trial_days) {
        
        if (bCluster == 1)  {  ## cRCT only
          num_to_enroll <- enrollment_schedule[(t-trial_startday+enrollment_gap)/enrollment_gap]
          # We try and enroll as many from the cluster as you can. I have set an
          # enrollment rate rather than cluster size, e.g. 70% enrolled in each cluster.
          # It means that every simulated trial would have slightly different numbers enrolled 
          # (=coverage*ave_community_size)
          
          # Need to choose from the clusters not already enrolled
          new_clusters <- sample(non_trial_clusters,num_to_enroll)
          new_clusters_v <- sample(new_clusters,num_to_enroll/2)
          new_clusters_c <- setdiff(new_clusters,new_clusters_v)
          # From the chosen clusters, a fraction of the non-infectious individual. That fraction is defined in the inputs
          # This looks complicated: for each new cluster, I sample from that cluster a proportion of the whole cluster,
          # but only from the susceptible or exposed individuals. If I'm trying to sample more than are available (because
          # there are lots of infectious/recovered individuals), just sample all of them.
          new_vacc <- unlist(lapply(new_clusters_v, function(x) {vertices <- g_name[g_community==x];
                                    possibles <- vertices[vertices%in%c(e_nodes[1,],s_nodes)]#intersect(vertices,c(e_nodes[1,],s_nodes))
                                    sample(possibles,min(round(cluster_coverage*length(vertices)),length(possibles)))
                                    }))
          new_controls <- unlist(lapply(new_clusters_c,
                                        function(x) {vertices <- g_name[g_community==x];
                                        possibles <- vertices[vertices%in%c(e_nodes[1,],s_nodes)]#intersect(vertices,c(e_nodes[1,],s_nodes))
                                        sample(possibles,min(round(cluster_coverage*length(vertices)),length(possibles)))
                                        }))
          
          enrolled_so_far <- nrow(trial_nodes_info)
          len_new_recruits <- length(new_controls)+length(new_vacc)
          trial_nodes_info[1:len_new_recruits+enrolled_so_far,] <- cbind(c(new_controls,new_vacc),
                                                                         g_community[c(new_controls,new_vacc)],
                                                                         c(rep(0,length(new_controls)),rep(1,length(new_vacc))),
                                                                         t, t)
          
          # Move the vaccinated susceptibles to from s_nodes to v_nodes
          vacc_susc <- intersect(s_nodes,new_vacc)
          s_nodes <- setdiff(s_nodes,vacc_susc)
          v_nodes <- c(v_nodes,vacc_susc)
          non_trial_clusters <- setdiff(non_trial_clusters,new_clusters)
        }
      }
      
      if (t%in%vaccination_days && bCluster==0) { ## staggered/instantaneous recruitment
        possibles <- vertices[vertices%in%c(setdiff(e_nodes[1,],c(trial_nodes_info$Node)),s_nodes)] # excludes v and c
        number_to_treat_rpois <- rpois(1,number_to_treat)
        new_recruits <- sample(possibles,number_to_treat_rpois)
        #V(g)[name %in% unlist(new_recruits)]$enrollmentday <- t
        new_vacc <- sample(new_recruits,allocation_rate*length(new_recruits))
        new_controls <- setdiff(new_recruits,new_vacc)
        
        enrolled_so_far <- nrow(trial_nodes_info)
        len_new_recruits <- length(new_controls)+length(new_vacc)
        trial_nodes_info[1:len_new_recruits+enrolled_so_far,] <- cbind(c(new_controls,new_vacc),
                                                                       g_community[c(new_controls,new_vacc)],
                                                                       c(rep(0,length(new_controls)),rep(1,length(new_vacc))),
                                                                       t, t)
        
        # Move the vaccinated susceptibles to from s_nodes to v_nodes
        vacc_susc <- intersect(s_nodes,new_vacc)
        s_nodes <- setdiff(s_nodes,vacc_susc)
        v_nodes <- c(v_nodes,vacc_susc)
        cont_susc <- intersect(s_nodes,new_controls)
        s_nodes <- setdiff(s_nodes,cont_susc)
        c_nodes <- c(c_nodes,cont_susc)
      }
    }
      
    # Only need to recover if there are any infected or exposed
    newinfectious <- c()
    if ((ncol(i_nodes)>0)||(ncol(e_nodes)>0)) 
      list[e_nodes,i_nodes,r_nodes,newinfectious] <- recover(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate)
    #sub_g <- if(ncol(i_nodes)>0) induced_subgraph(g,c(i_nodes[1,],s_nodes,v_nodes)) else NULL
    list[s_nodes,v_nodes,e_nodes,c_nodes] <- spread(sub_g=g,g_community,s_nodes,v_nodes,e_nodes,i_nodes,c_nodes,
             beta,direct_VE,incperiod_shape,incperiod_rate,connected_nodes=connected_to_source,external_inf_F=extF,source_num_inf=infected_trajectory[t])
    numnewinfectious <- length(newinfectious)
    if (numnewinfectious>0) {
      new_indices <- g_name %in% newinfectious
      v_subset <- V(g)[new_indices]
      match_indices <- match(newinfectious,trial_nodes_info$Node)
      if(length(newinfectious)!=length(v_subset)||length(match_indices)!=length(v_subset)||length(match_indices)!=length(newinfectious)){
        browser()
      }
      # Update results
      results <- rbind(results,data.frame("InfectedNode"=newinfectious,
                            "DayInfected"=t,
                            "Community"=v_subset$community,
                            "TrialStatus"=trial_nodes_info$TrialStatus[match_indices],
                            "DayEnrolled"=trial_nodes_info$DayEnrolled[match_indices],
                            "DayVaccinated"=trial_nodes_info$DayVaccinated[match_indices]))
      
      numinfectious <- numinfectious+numnewinfectious
    }
    
    if(bTrial==2){ ## ring vaccination
      if(t>=trial_startday&numnewinfectious>0){
        # get all infected people's neighbours
        all_contacts <- unlist(ego(g,order=1,nodes=newinfectious))
        untrialled_contacts <- all_contacts[!all_contacts%in%trial_nodes_info$Node]
        susc_contacts <- unique(untrialled_contacts[untrialled_contacts%in%c(s_nodes,e_nodes[1,])])
        if(length(susc_contacts)>0){
          new_recruits <- if(length(susc_contacts)==1){susc_contacts}else{base::sample(susc_contacts,round(length(susc_contacts)*cluster_coverage))}
          new_vacc <- if(length(new_recruits)==1&runif(1)<allocation_rate){new_recruits}else{base::sample(new_recruits,round(length(new_recruits)*allocation_rate))}
          new_controls <- setdiff(new_recruits,new_vacc)
          
          enrolled_so_far <- nrow(trial_nodes_info)
          len_new_recruits <- length(new_recruits)
          trial_nodes_info[1:len_new_recruits+enrolled_so_far,] <- cbind(c(new_controls,new_vacc),
                                                                         g_community[c(new_controls,new_vacc)],
                                                                         c(rep(0,length(new_controls)),rep(1,length(new_vacc))),
                                                                         t, t)
          # Move the vaccinated susceptibles to from s_nodes to v_nodes
          vacc_susc <- intersect(s_nodes,new_vacc)
          s_nodes <- setdiff(s_nodes,vacc_susc)
          v_nodes <- c(v_nodes,vacc_susc)
          cont_susc <- intersect(s_nodes,new_controls)
          s_nodes <- setdiff(s_nodes,cont_susc)
          c_nodes <- c(c_nodes,cont_susc)
        }
      }
    }
    
    trajectories$S <- c(trajectories$S,length(s_nodes) + length(v_nodes) + length(c_nodes))
    trajectories$E <- c(trajectories$E,ifelse(length(trajectories$E)==0,0,-diff(tail(trajectories$S,2))))
    trajectories$I <- c(trajectories$I,numnewinfectious)
    trajectories$R <- c(trajectories$R,ifelse(length(trajectories$R)==0,0,length(r_nodes)-sum(trajectories$R)))
  }
  
  # Tidy up results
  if (numinfectious>0) {
    results<-results[1:numinfectious,]
  } else {
    results <- results[1,]
  }
  
  list(results,trial_nodes_info,trajectories,allocation_rates)
}

## extracted from hitchings
get_infected_trajectory <- function(times){
  N <- 50000
  y <- c(S=N-1,E1=0,E2=0,E3=0,I1=1,I2=0,I3=0,R=0)
  parms <- c(betahat=0.94,a1=0.19,a2=0.6,atau=27.79,sigma=0.14,gamma=0.33)
  out <- as.data.frame(lsoda(y,times,source_population_model,parms))
  out$I1 + out$I2 + out$I3
}

# Recover and spread functions
## from hitchings
recover<-function(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate) {
  # Input is a list of the exposed nodes, 
  # with number of days since infection and total incubation/latent
  # period, and equivalently for the infectious nodes.
  # For each of these nodes, we will add it to newinfectious if the number of days exposed has
  # reached the total length of the incubation period, and equivalently for the infectious nodes.
  
  # Advance infectious nodes first, otherwise you will doubly advance any nodes switching from
  # exposed to infectious at this time step
  indices_to_remove <- i_nodes[2,]>=i_nodes[3,]
  newremoved <- as.vector(i_nodes[1,indices_to_remove])
  
  # Add one day to the length of each infected individual's time infected
  i_nodes[2,] <- i_nodes[2,]+1
  
  # Remove any recovered from i_nodes and add to r_nodes
  i_nodes <- i_nodes[,!(i_nodes[1,] %in% newremoved),drop=FALSE]
  r_nodes <- c(r_nodes,newremoved)
  
  # Now advance exposed nodes
  indices_to_remove <- e_nodes[2,]>=e_nodes[3,]
  newinfectious <- as.vector(e_nodes[1,indices_to_remove])
  
  # Add one day to the length of each infected individual's time infected
  e_nodes[2,] <- e_nodes[2,]+1
  
  # Remove any progressing from e_nodes and add to i_nodes
  e_nodes <- e_nodes[,!(e_nodes[1,] %in% newinfectious),drop=FALSE]
  inf_periods <- rgamma(length(newinfectious),infperiod_shape,infperiod_rate)
  i_nodes <- cbind(i_nodes,rbind(newinfectious,rep(0,length(newinfectious)),inf_periods))
  
  list(e_nodes, i_nodes, r_nodes, sort(newinfectious))
}

## from hitchings
spread<-function(sub_g,g_community, s_nodes, v_nodes, e_nodes, i_nodes,c_nodes, beta, direct_VE,
                 incperiod_shape, incperiod_rate,connected_nodes,external_inf_F,source_num_inf){
  # Spread will create new infected nodes from two sources: infectious nodes within the the study
  # population, and external pressure from the source population
  # Inputs:
  # g is the graph, used to find neighbours of infected nodes
  # s_nodes, e_nodes and i_nodes are susceptible, exposed and infected nodes
  # beta is the hazard of infection for one contact
  # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
  # length, currently drawn from a gamma distribution
  # connected_nodes is a list of nodes that are connected to the source population
  # external_inf_F is a constant of proportionality that defines infectious pressure from source population 
  # to an individual
  # source_num_inf is the number of infectious individuals in the source population
  
  # Process: go through list of i_nodes, and choose a random number of its susceptible neighbours to
  # be infected, according to beta and choose a random number of its susceptible vaccinated neighbours to
  # be infected, according to beta and direct_VE
  # Then go through list of nodes that are connected to source population and infect each susceptible
  # one with probability 1-exp(-FI), where I is the number/proportion of infectious, and F is a constant
  infectees_susc <- c()
  infectees_vacc <- c()
  infectees_cont <- c()
  if (ncol(i_nodes)>0) {
    # Make a beta vector
    beta_v <- beta*(1-direct_VE)
    # Get a list of all neighbours of all infected nodes
    potential_contacts <- unlist(ego(g,order=1,nodes=i_nodes[1,]))
    infectees_susc <- infect_neighbours(potential_contacts,node_class=s_nodes,beta_value=beta)
    if (length(v_nodes)>0) 
      infectees_vacc <- infect_neighbours(potential_contacts,node_class=v_nodes,beta_value=beta_v)
    if (length(c_nodes)>0) 
      infectees_cont <- infect_neighbours(potential_contacts,node_class=c_nodes,beta_value=beta)
  } 
  potential_connected_nodes <- connected_nodes[!connected_nodes%in%c(infectees_susc,infectees_vacc,infectees_cont)]
  conn_inf_susc <- infect_from_source(g_community,num_communities,connected_nodes=potential_connected_nodes,target_nodes=s_nodes,direct_VE=0,
                                      source_num_inf,external_inf_F)
  conn_inf_vacc <- infect_from_source(g_community,num_communities,connected_nodes=potential_connected_nodes,target_nodes=v_nodes,direct_VE=direct_VE,
                                      source_num_inf,external_inf_F)
  conn_inf_cont <- infect_from_source(g_community,num_communities,connected_nodes=potential_connected_nodes,target_nodes=c_nodes,direct_VE=0,
                                      source_num_inf,external_inf_F)
  
  newinfected_susc <- c(infectees_susc,conn_inf_susc)
  newinfected_vacc <- c(infectees_vacc,conn_inf_vacc)
  newinfected_cont <- c(infectees_cont,conn_inf_cont)
  newinfected <- c(newinfected_susc, newinfected_vacc,newinfected_cont)
  newinfected <- unique(newinfected)
  
  if (length(newinfected)>0) {
    # Give each newly exposed node an incubation/latent period
    inc_periods <- rgamma(length(newinfected),incperiod_shape,incperiod_rate)
    # Add them to e_nodes and remove from s_nodes and v_nodes
    e_nodes <- cbind(e_nodes,rbind(newinfected,rep(0,length(newinfected)),inc_periods))
    s_nodes<-setdiff(s_nodes,newinfected_susc)
    v_nodes <- setdiff(v_nodes,newinfected_vacc)
    c_nodes <- setdiff(c_nodes,newinfected_cont)
  }
  list(s_nodes, v_nodes, e_nodes,c_nodes)
}

## extracted from hitchings
infect_from_source <- function(g_community,num_communities,connected_nodes,target_nodes,direct_VE,source_num_inf,external_inf_F){
  # Pick out the nodes connected to the source that are still susceptible 
  # and haven't just been infected
  target_cnodes <- target_nodes[target_nodes%in%connected_nodes]
  conn_inf_susc <- c()
  if (length(target_cnodes)>0) {
    # Make a vector to represent external infection hazard for each individual
    communities <- g_community[target_cnodes]
    #comm_sizes <- sapply(1:num_communities,function(x) sum(communities==x))
    # Hazard of infection
    extFs <- external_inf_F[communities]#rep(external_inf_F,comm_sizes)
    # Probability of infection
    prob_inf_fromsource <- 1 - exp(-(1-direct_VE)*mean(extFs)*source_num_inf)
    
    ## 
    #pdf('prob_inf_fromsource.pdf',width=6,height=4); par(mar=c(5,5,2,2))
    #plot(1:length(source_num_inf),1 - exp(-(1-direct_VE)*external_inf_F[1]*source_num_inf),typ='l',col='navyblue',ylim=c(0,3e-5),
    #     frame=F,cex.lab=1.5,cex.axis=1.5,xlab='Day',ylab=TeX('$G_i(t)$'))
    #for(i in 2:length(external_inf_F)) lines(1:length(source_num_inf),1 - exp(-(1-direct_VE)*external_inf_F[i]*source_num_inf),typ='l',col='navyblue')
    #dev.off()
    
    #pdf('sum_prob_inf_fromsource.pdf',width=4,height=4); par(mar=c(5,5,2,2))
    #plot(sapply(1:length(external_inf_F),function(i)sum(communities==i)),
    #     sapply(1:length(external_inf_F),function(i)sum(1 - exp(-(1-direct_VE)*external_inf_F[i]*source_num_inf))*sum(communities==i)),
    #     frame=F,cex.lab=1.5,cex.axis=1.5,xlab='Community size',ylab='Mean infections from source',col='navyblue')
    #dev.off()
    ##
    
    # Choose a number of individuals to be infected, then sample those individuals
    num_conn_inf_susc <- rbinom(1,length(target_cnodes),prob_inf_fromsource)
    if(num_conn_inf_susc>0) conn_inf_susc <- sample(target_cnodes,num_conn_inf_susc,prob=extFs)
  }
  conn_inf_susc
}

## extracted from hitchings
infect_neighbours <- function(potential_contacts,node_class,beta_value){
  susc_contacts <- potential_contacts[potential_contacts%in%node_class]#intersect(potential_contacts,node_class)
  num_neighbours_susc <- length(susc_contacts)
  # Sample from each group of neighbours in turn
  # First choose how many neighbours each node infects
  num_contacts_susc <- rbinom(1,num_neighbours_susc,1-exp(-beta_value))
  # Then sample from the neighbours
  # If one node gets picked twice by different nodes, just discard the duplicate.
  infectees_temp <- sample(susc_contacts,num_contacts_susc)
  infectees_susc <- unique(infectees_temp)
  return(infectees_susc)
}

# function to run THE EPIDEMIC IN THE SOURCE POPULATION
# This is to define external infectious pressure to the network
## from hitchings
source_population_model <- function(t, y, parms) {
  with(as.list(c(y,parms)), {
    
    beta <- betahat * (1 - a2/(1 + exp(-a1 * (t - atau))))
    
    dS <- -beta * S * (I1+I2+I3) / (S+E1+E2+E3+I1+I2+I3+R)
    dE1 <- beta * S * (I1+I2+I3) / (S+E1+E2+E3+I1+I2+I3+R) - sigma * 3 * E1
    dE2 <- sigma * 3 * E1 - sigma * 3 * E2
    dE3 <- sigma * 3 * E2 - sigma * 3 * E3
    dI1 <- sigma * 3 * E3 - gamma * 3 * I1
    dI2 <- gamma * 3 * I1 - gamma * 3 * I2
    dI3 <- gamma * 3 * I2 - gamma * 3 * I3
    dR <- gamma * 3 * I3
    list(c(dS,dE1,dE2,dE3,dI1,dI2,dI3,dR))
  })
}

## adapted from hitchings
analyse_data <- function(results,trial_nodes,trial_startday,trial_length,ave_inc_period,
                         bCluster,follow_up,revisit) {
  if(revisit==1) follow_up <- trial_length
  excluded0 <- sum(results$TrialStatus==0&results$DayInfected-results$DayVaccinated<=ave_inc_period,na.rm=T)
  excluded1 <- sum(results$TrialStatus==1&results$DayInfected-results$DayVaccinated<=ave_inc_period,na.rm=T)
  fail0 <- sum(results$TrialStatus==0&(results$DayInfected<follow_up+results$DayVaccinated),na.rm=T) - excluded0 #
  fail1 <- sum(results$TrialStatus==1&(results$DayInfected<follow_up+results$DayVaccinated),na.rm=T) - excluded1 #
  n0 <- ifelse(nrow(trial_nodes)==0,0,sum(trial_nodes==0,na.rm=T)) - excluded0
  n1 <- ifelse(nrow(trial_nodes)==0,0,sum(trial_nodes==1,na.rm=T)) - excluded1
  success0 <- n0 - fail0
  success1 <- n1 - fail1
  p0 <- success0/n0
  p1 <- success1/n1
  sigma0 <- p0 * ( 1 - p0 ) /n0
  sigma1 <- p1 * ( 1 - p1 ) /n1
  zval <- (p1-p0)/(sqrt(sigma0+sigma1))
  pval_binary_mle <- dnorm(zval)
  
  VE_pointest_binary_mle <- 1 - (fail1/n1)/(fail0/n0)
  
  results$DayInfected <- results$DayInfected - results$DayVaccinated
  
  # Get a list of nodes that were enrolled in the trial but never infected
  noninf <- setdiff(trial_nodes$Node,results$InfectedNode)
  # Get list of nodes that became infectious while they were in the trial
  # This is the step that excludes those who were R at time of enrollment
  results_analysis <- results[!is.na(results$TrialStatus),]
  # Get a list of nodes who were infected after their follow-up time was over
  # (i.e. those enrolled at the beginning but infected right at the end)
  #censored <- results_analysis[results_analysis$DayInfected>trial_length,]
  #results_analysis <- results_analysis[results_analysis$DayInfected<=trial_length,]
  # Assign them eventstatus=1 for the Cox analysis
  #results_analysis$eventstatus <- 1
  # Make data frame for those who were never infected (i.e. censored by end of study)
  #noninfdf <- data.frame(InfectedNode=noninf,DayInfected=rep(trial_length,length(noninf)),
  #                     Community=trial_nodes$Community[trial_nodes$Node %in% noninf],
  #                     TrialStatus=trial_nodes$TrialStatus[trial_nodes$Node %in% noninf],
  #                     DayEnrolled=trial_nodes$DayEnrolled[trial_nodes$Node %in% noninf],
  #                     DayVaccinated=trial_nodes$DayVaccinated[trial_nodes$Node %in% noninf],
  #                     eventstatus=0)
  # if (nrow(censored)>0) {
  #   censored$DayInfected <- trial_length
  #   censored$eventstatus <- 0
  # }
  # Remove column with simulation number so the columns match up
  #results_analysis$SimulationNumber <- NULL
  #results_analysis <- rbind(results_analysis,noninfdf)#,censored)
  
  # Finally, exclude any cases who were infected during the first n days of follow-up
  # This tries to rid of those who were already latently infected when enrolled
  results_analysis <- results_analysis[results_analysis$DayInfected>ave_inc_period,]
  
  numevents_vacc <- nrow(results_analysis[ (results_analysis$TrialStatus==1),])
  numevents_cont <- nrow(results_analysis[ (results_analysis$TrialStatus==0),])
  
  #total_vacc_pt <- sum(results_analysis$DayInfected[results_analysis$TrialStatus==1])
  #total_cont_pt <- sum(results_analysis$DayInfected[results_analysis$TrialStatus==0])
  #VE_pointest <- 1 - (numevents_vacc/total_vacc_pt)/(numevents_cont/total_cont_pt)
  
  sample_size <- nrow(results_analysis) + length(noninf)
  
  pval <- pval_binary_mle
  vaccEffEst <- VE_pointest_binary_mle
  return(list(vaccEffEst,pval,numevents_vacc,numevents_cont,sample_size))
  
  # Calculate ICC
  # Number of events in each cluster
  events_by_cluster <- aggregate(results_analysis$eventstatus,
                               by=list(Community=results_analysis$Community),FUN=sum)
  # Size of each cluster
  cluster_sizes <- aggregate(results_analysis$InfectedNode,by=list(Community=results_analysis$Community),FUN=length)
  # Overall trial size
  N <- sum(cluster_sizes$x)
  # Overall number of clusters
  K <- nrow(cluster_sizes)
  n0 <- 1/(K-1) * (N - sum(cluster_sizes$x^2)/N)
  n01 <- 1/(K-1) * ((K-1)*n0 - sum(cluster_sizes$x^2)/N)
  MSB <- 1/(K-1) * sum((events_by_cluster$x-mean(events_by_cluster$x))^2/cluster_sizes$x)
  MSW <- 1/(N-K-1) * sum(events_by_cluster$x-events_by_cluster$x^2/cluster_sizes$x)
  ICC <- (MSB - MSW) / (MSB + (n01-1) * MSW)
  deff <- 1+(mean(cluster_sizes$x)-1)*ICC
  
  # Proportion of clusters that have zero cases
  prop_zeros <- sum(events_by_cluster$x==0)/K
  
  if (bCluster == 0) {
    # If no events are observed in either arm, the trial has failed and no result can be obtained
    vaccEffEst <- c(NA,NA,NA)
    pval <- NA
    # Analysis for iRCT
    if ((numevents_vacc>0)&&(numevents_cont>0)) {
      # If we have events in both arms, can try a Cox PH. It can still be singular, so if it
      # throws an error, the trial has failed (not enough events)
      list[vaccEffEst,pval] <- coxmodel(results_analysis,VE_pointest)
      
    } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
      # If there are no events in the control arm but
      # events in the vaccine arm, VE estimate is -1 and p-value is 1
      vaccEffEst <- c(-1,-1,-1)
      pval <- 1
      
    } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
      # for the p-value we add one event to both arms and do a Cox regression on that data
      # I give both events the median time among control events.
      newevent_v_rownum <- min(which((results_analysis$eventstatus==0)&(results_analysis$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis$eventstatus==0)&(results_analysis$TrialStatus==0)))
      
      eventtime <- median(results_analysis$DayInfected[results_analysis$eventstatus==1])
      
      results_analysis$DayInfected[c(newevent_v_rownum,newevent_c_rownum)] <- eventtime
      results_analysis$eventstatus[c(newevent_v_rownum,newevent_c_rownum)] <- 1
      
      list[,pval] <- coxmodel(results_analysis,VE_pointest)
      vaccEffEst<-1
    } 
    return(list(vaccEffEst,pval,numevents_vacc,numevents_cont,sample_size))
    
  } else {
    # If no events are observed in either arm, the trial has failed and no result can be obtained
    vaccEffEst_gaussian_coxme <- c(NA,NA,NA)
    pval_gaussian_coxme <- NA
    vaccEffEst_gaussian_coxph <- c(NA,NA,NA)
    pval_gaussian_coxph <- NA
    vaccEffEst_gamma_coxph <- c(NA,NA,NA)
    pval_gamma_coxph <- NA
    vaccEffEst_gee <- c(NA,NA,NA)
    pval_gee <- NA
    
    if ((numevents_vacc>0)&&(numevents_cont>0)) {
      # Run the models for the clustered data
      list[vaccEffEst_gaussian_coxme,pval_gaussian_coxme,
           vaccEffEst_gaussian_coxph,pval_gaussian_coxph,
           vaccEffEst_gamma_coxph,pval_gamma_coxph,
           vaccEffEst_gee,pval_gee] <- clustermodels(results_analysis,VE_pointest)
      
    } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
      # If there are no events in the control arm but
      # events in the vaccine arm, VE estimate is -1 and p-value is 1
      vaccEffEst_gaussian_coxme <- c(-1,-1,-1)
      pval_gaussian_coxme <- 1
      vaccEffEst_gaussian_coxph <- c(-1,-1,-1)
      pval_gaussian_coxph <- 1
      vaccEffEst_gamma_coxph <- c(-1,-1,-1)
      pval_gamma_coxph <- 1
      vaccEffEst_gee <- c(-1,-1,-1)
      pval_gee <- 1
      
    } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1
      # For p-value, add one event to both arms. Give it the median event time among control events.
      # Cluster for event in vaccine arm doesn't matter
      # Choose cluster for event in control arm to be most conservative
      
      communities <- results_analysis$Community[((results_analysis$eventstatus==1)&(results_analysis$TrialStatus==0))]
      freq_table <- sort(table(communities),decreasing=TRUE)
      community <- as.numeric(names(freq_table[1]))
      
      newevent_v_rownum <- min(which((results_analysis$eventstatus==0)&(results_analysis$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis$eventstatus==0)&(results_analysis$TrialStatus==0)&(results_analysis$Community==community)))
      
      eventtime <- median(results_analysis$DayInfected[results_analysis$eventstatus==1])
      
      results_analysis$DayInfected[c(newevent_v_rownum,newevent_c_rownum)] <- eventtime
      results_analysis$eventstatus[c(newevent_v_rownum,newevent_c_rownum)] <- 1
      
      list[,pval_gaussian_coxme,
           ,pval_gaussian_coxph,
           ,pval_gamma_coxph,
           ,pval_gee] <- clustermodels(results_analysis,VE_pointest)
      vaccEffEst_gaussian_coxme <- 1
      vaccEffEst_gaussian_coxph <- 1
      vaccEffEst_gamma_coxph <- 1
      vaccEffEst_gee <- 1
    }
    pval_gee <- pval_binary_mle
    vaccEffEst_gee <- VE_pointest_binary_mle
    return(list(vaccEffEst_gaussian_coxme,pval_gaussian_coxme,
                vaccEffEst_gaussian_coxph,pval_gaussian_coxph,
                vaccEffEst_gamma_coxph,pval_gamma_coxph,
                vaccEffEst_gee,pval_gee,
                numevents_vacc,numevents_cont,sample_size,ICC,deff,prop_zeros))
  }
}

## from hitchings
coxmodel <- function(data,VEpointest) {
  vaccEffEst<-c(VEpointest,NA,NA)
  pval <- NA
  
  survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus+strata(Community),data),silent=T)
  usesurvmod <- !inherits(survmodel, 'try-error')
  
  if (usesurvmod && vcov(survmodel)>=0){
    # If no error was thrown and the variance is positive, use the results of the model
    vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
    zval <- survmodel$coefficient/sqrt(survmodel$var)
    pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
  } 
  list(vaccEffEst,pval)
}

## from hitchings
clustermodels <- function(data,VEpointest) {
  
  # Gaussian shared frailty, using coxme
  mod <- try(coxme(Surv(DayInfected, eventstatus) ~ TrialStatus + (1|Community), data = data),silent=T)
  usemod <- !inherits(mod, 'try-error')
  
  if (usemod && is.na(vcov(mod))) {
    vaccEffEst_gaussian_coxme <- c(VEpointest,NA,NA)
    pval_gaussian_coxme<-NA
  } else {
    if (usemod && vcov(mod)>=0){
      # If the model converges and variance is positive, use results of model
      
      vaccEffEst_gaussian_coxme <- 1-exp(mod$coefficients + c(0, 1.96, -1.96)*as.vector(sqrt(vcov(mod))))
      zval <- mod$coefficients/sqrt(vcov(mod))
      pval_gaussian_coxme <- pnorm(zval, lower.tail = vaccEffEst_gaussian_coxme[1]>0)*2
    } else {
      vaccEffEst_gaussian_coxme <- c(VEpointest,NA,NA)
      pval_gaussian_coxme<-NA
    }
  }
  
  # Try with a gaussian-distributed random effect.
  mod2 <- try(coxph(Surv(DayInfected, eventstatus) ~ TrialStatus + frailty(Community,distribution="gaussian",sparse=FALSE), data = data), silent=T)
  usemod2 <- !inherits(mod2, 'try-error')
  
  if (usemod2 && is.na(vcov(mod2))) {
    vaccEffEst_gaussian_coxph <- c(VEpointest,NA,NA)
    pval_gaussian_coxph<-NA
  } else {
    if (usemod2 && vcov(mod2)>=0){
      
      vaccEffEst_gaussian_coxph <- 1-exp(mod2$coefficients[1] + c(0, 1.96, -1.96)*as.vector(sqrt(vcov(mod2)[1])))
      zval2 <- mod2$coefficients[1]/sqrt(vcov(mod2)[1])
      pval_gaussian_coxph <- pnorm(zval2, lower.tail = vaccEffEst_gaussian_coxph[1]>0)*2
    } else {
      vaccEffEst_gaussian_coxph <- c(VEpointest,NA,NA)
      pval_gaussian_coxph<-NA
    }
  }
  
  # Try with a gamma-distributed random effect, using coxph and frailty
  mod3 <- try(coxph(Surv(DayInfected, eventstatus) ~ TrialStatus + frailty(Community,distribution="gamma",sparse=FALSE), data = data), silent=T)
  usemod3 <- !inherits(mod3, 'try-error')
  
  if (usemod3 && is.na(vcov(mod3))) {
    vaccEffEst_gamma_coxph <- c(VEpointest,NA,NA)
    pval_gamma_coxph<-NA
  } else {
    if (usemod3 && vcov(mod3)>=0){
      
      vaccEffEst_gamma_coxph <- 1-exp(mod3$coefficients[1] + c(0, 1.96, -1.96)*as.vector(sqrt(vcov(mod3)[1])))
      zval3 <- mod3$coefficients[1]/sqrt(vcov(mod3)[1])
      pval_gamma_coxph <- pnorm(zval3, lower.tail = vaccEffEst_gamma_coxph[1]>0)*2
    } else {
      vaccEffEst_gamma_coxph <- c(VEpointest,NA,NA)
      pval_gamma_coxph<-NA
    }
  }
  
  # Generalized estimating equations/robust variance
  mod5 <- try(coxph(Surv(DayInfected, eventstatus) ~ TrialStatus + cluster(Community), data = data), silent=T)
  usemod5 <- !inherits(mod5, 'try-error')
  
  if (usemod5 && is.na(vcov(mod5))) {
    vaccEffEst_gee <- c(VEpointest,NA,NA)
    pval_gee<-NA
  } else {
    if (usemod5 && vcov(mod5)>=0){
      
      vaccEffEst_gee <- 1-exp(mod5$coefficients[1] + c(0, 1.96, -1.96)*as.vector(sqrt(vcov(mod5)[1])))
      zval5 <- mod5$coefficients[1]/sqrt(vcov(mod5)[1])
      pval_gee <- pnorm(zval5, lower.tail = vaccEffEst_gee[1]>0)*2
    } else {
      vaccEffEst_gee <- c(VEpointest,NA,NA)
      pval_gee<-NA
    }
  }
  
  list(vaccEffEst_gaussian_coxme,pval_gaussian_coxme,
       vaccEffEst_gaussian_coxph,pval_gaussian_coxph,
       vaccEffEst_gamma_coxph,pval_gamma_coxph,
       vaccEffEst_gee,pval_gee)
  
}

# In vaccinated clusters, the CI among unvaccinated and vaccinated
## from hitchings
finalsizes_vacc <- function(CIs,parms) {
  finalsize_U <- 1 - exp(-parms$R0 * ((1-parms$enrollment_perc) * CIs[1] + parms$enrollment_perc * CIs[2])) - CIs[1]
  finalsize_V <- 1 - exp(-parms$R0 * (1-parms$direct_VE) * 
                           ((1-parms$enrollment_perc) * CIs[1] + parms$enrollment_perc * CIs[2])) - CIs[2]
  c(F1=finalsize_U,F2=finalsize_V)
}

# In unvaccinated clusters, the CI among the unvaccinated
## from hitchings
finalsize_unvacc <- function(CI,parms) {
  finalsize_U <- 1 - exp(-parms$R0 * CI) - CI
  return(finalsize_U)
}

# Function to solve for outbreak probability in vaccinated clusters
## from hitchings
outbreakprob_vacc <- function(x,parms) {
  prob <- exp(-parms$R0_vacc*(1-x)) - x
  return(prob)
}


# CI in iRCT
## from hitchings
finalsizes_iRCT <- function(CIs,parms) {
  finalsize_U <- 1 - exp(-parms$R0 * (parms$enrollment_perc*0.5 * CIs[2] + (1-parms$enrollment_perc*0.5) * CIs[1])) - CIs[1]
  finalsize_V <- 1 - exp(-parms$R0 * (1-parms$direct_VE) * 
                           (parms$enrollment_perc*0.5 * CIs[2] + (1-parms$enrollment_perc*0.5) * CIs[1])) - CIs[2]
  c(F1=finalsize_U,F2=finalsize_V)
}

# Probability of an outbreak in a cluster (only if R0>1)
## from hitchings
outbreakprob_iRCT <- function(x,parms) {
  prob <- exp(-parms$R0_iRCT*(1-x)) - x
  return(prob)
}
