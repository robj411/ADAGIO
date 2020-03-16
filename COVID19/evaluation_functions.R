get_weighted_results <- function(results,how_surprising=c()){
  # the day the cluster is recruited
  recruit_day <- results$RecruitmentDay[1]
  # the true nunber infected after the recruitment day
  true_positives <- sum(results$DayInfected>recruit_day)
  # the number we estimate with a binary weight
  binary <- sum(results$DayInfectious>recruit_day+9)
  # the day individuals became infectious
  days_infectious <- results$DayInfectious
  # the subset of those that might have infected others
  infectors <- days_infectious[1:(nrow(results)-1)]
  # the durations for which they were infectious
  infector_durations <- results$DayRemoved[1:(nrow(results)-1)] - infectors
  infector_durations[is.na(infector_durations)|infector_durations>21] <- 21
  # those who were infected by someone else
  infectees <- days_infectious[2:(nrow(results))]
  infector_names <- results$InfectedNode[1:(nrow(results)-1)]
  infectee_names <- results$InfectedNode[2:(nrow(results))]
  weight <- weight_hh <- weight_hh_rem <- 0
  for(j in 1:length(infectees)){
    # which infectors might have infected infectee j
    infectors_for_j <- infectors<infectees[j]
    # rows: the time lag between infector and infectee becoming infectious
    rows <- pmin(infectees[j]-infectors[infectors_for_j],nrow(probability_by_lag))
    # cols: the day the potential infector became infectious relative to recruitment day
    cols <- pmax(ref_recruit_day-recruit_day+infectors[infectors_for_j],1)
    # probabilities for infectors to infect infectee j
    prob_infectors <- probability_by_lag[cbind(rows,cols)]
    #if(results$inCluster[j+1]==T){
    # how surprising is it that infectee j became infectious?
    how_surprising[length(how_surprising)+1] <- prod(1-prob_infectors[sapply(infector_names[infectors_for_j],
                                                                             function(x)x%in%c(contact_list[[infectee_names[j]]],contact_of_contact_list[[infectee_names[j]]])
    )])
    #}
    # probabilities infected after recruitment day given infected by infector
    prob_after_0 <- probability_after_day_0[cbind(rows,cols)]
    # weights for all relationships given known network
    hh_weight <- sapply(infector_names[infectors_for_j],
                        function(x)as.numeric(x%in%contact_list[[infectee_names[j]]]) +
                          as.numeric(x%in%household_list[[infectee_names[j]]])*(high_risk_scalar-1) +
                          as.numeric(x%in%high_risk_list[[infectee_names[j]]])*(high_risk_scalar-1) +
                          as.numeric(x%in%contact_of_contact_list[[infectee_names[j]]])*neighbour_scalar
    )
    # weighted probability
    normalised_prob_infectors <- prob_infectors/sum(prob_infectors)
    weight <- weight + sum(prob_after_0*normalised_prob_infectors)
    # weighted probability using contact information
    normalised_prob_infectors <- prob_infectors*hh_weight/sum(prob_infectors*hh_weight)
    weight_hh <- weight_hh + sum(prob_after_0*normalised_prob_infectors)
    # weighted probability using contact and removal information
    #prob_infectors <- sapply(1:length(rows),function(x)probability_by_lag_given_removal[[max(infector_durations[x]-1,1)]][rows[x],cols[x]])
    prob_infectors <- sapply(1:length(rows),function(x)probability_by_lag_given_removal_mat[rows[x],max(infector_durations[x]-1,1)])
    prob_after_0 <- sapply(1:length(rows),function(x)probability_after_day_0_given_removal[[max(infector_durations[x]-1,1)]][rows[x],cols[x]])
    normalised_prob_infectors <- prob_infectors*hh_weight/sum(prob_infectors*hh_weight)
    weight_hh_rem <- weight_hh_rem + sum(prob_after_0*normalised_prob_infectors)
  }
  return(list(c(true_positives,weight,binary,weight_hh,weight_hh_rem),how_surprising))
}

get_weighted_results_given_ve <- function(results,ve_point_est,tested=F){
  weight_hh_rem <- colSums(get_infectee_weights(results,ve_point_est,tested)[[1]])
  return(weight_hh_rem)
}

get_infectee_weights <- function(results,ve_point_est,tested=F){
  
  # the day the cluster is recruited
  recruit_day <- results$RecruitmentDay
  # the day individuals became infectious
  days_infectious <- results$DayInfectious
  weight_hh_rem <- matrix(0,ncol=2,nrow=1)
  infectee_names <- c()
  # those who were infected by someone else
  infectee_index <- !is.na(recruit_day) & days_infectious>recruit_day
  if(sum(infectee_index)>0){
    weight_hh_rem <- matrix(0,ncol=2,nrow=sum(infectee_index))
    infectees <- days_infectious[infectee_index] - recruit_day[infectee_index]
    infectee_names <- results$InfectedNode[infectee_index]
    infectee_trial <- results$inTrial[infectee_index]
    infectee_vaccinated <- results$vaccinated[infectee_index]
    for(j in 1:length(infectees)){
      if(infectee_trial[j]){
        prob_after_0 <- pgamma(infectees[j],shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
        if(tested) {
          # if test positive, probability=0
          if(c(results$DayInfected[infectee_index])[j]<c(results$RecruitmentDay[infectee_index])[j]){
            prob_after_0 <- 0
          }else{
            prob_after_0 <- pgamma(infectees[j],shape=vacc_shape,rate=vacc_rate)
          }
        }
        # store complement
        yij <- 1 - prob_after_0
        # if vaccinated, adjust probability to be infected after day 0
        if(infectee_vaccinated[j]) prob_after_0 <- (1-ve_point_est)*prob_after_0/(yij+(1-ve_point_est)*prob_after_0)
        ## something about a test if infected on day 0
        # add to weight for vaccinated or unvaccinated
        if(infectee_vaccinated[j]){
          weight_hh_rem[j,1] <- prob_after_0
        }else{
          weight_hh_rem[j,2] <- prob_after_0
        }
      }
    }
  }
  return(list(weight_hh_rem,infectee_names))
}

calculate_pval <- function(fails,sizes){
  fail1 <- fails[1]
  fail0 <- fails[2]
  n1 <- sizes[1]
  n0 <- sizes[2]
  success0 <- n0 - fail0
  success1 <- n1 - fail1
  p0 <- success0/n0
  p1 <- success1/n1
  sigma0 <- p0 * ( 1 - p0 ) /n0
  sigma1 <- p1 * ( 1 - p1 ) /n1
  zval <- (p1-p0)/(sqrt(sigma0+sigma1))
  dnorm(zval)
}

calculate_ve <- function(fails,sizes){
  fail1 <- fails[1]
  fail0 <- fails[2]
  n1 <- sizes[1]
  n0 <- sizes[2]
  1 - (fail1/n1)/(fail0/n0)
}

response_adapt <- function(results_list,vaccinees,trial_participants, adaptation='TST',func=get_efficacious_probabilities){
  probs <- func(results_list,vaccinees,trial_participants,max_time=length(results_list))
  pop_sizes2 <- probs[[2]]
  fails <- probs[[3]]
  successes <- pop_sizes2 - fails
  
  if(adaptation%in%c('Ros','Ney')){
    ps <- successes/pop_sizes2
    if(adaptation=='Ros'){
      R_val <- sqrt(ps[1]/ps[2])  # ros
      allocation_rate <- R_val / (1+R_val)
    }else if(adaptation=='Ney'){
      allocation_rate <- ifelse(any(ps*(1-ps)==0), 0.5, sqrt(ps[1]*(1-ps[1])) / (sqrt(ps[2]*(1-ps[2]))+ sqrt(ps[1]*(1-ps[1]))) )# ney
    }
  }else if(adaptation%in%c('TS','TST')){
    j <- length(results_list) # t - trial_startday
    bigT <- nClusters # trial_length
    tuning_c <- ifelse(adaptation=='TS',1,(j/bigT))
    #print(tuning_c)
    p0 <- rbeta(1000,1+successes[2],1+fails[2])
    p1 <- rbeta(1000,1+successes[1],1+fails[1])
    prob1 <- sum(p1>p0)/1000
    allocation_rate <- prob1^tuning_c / (prob1^tuning_c + (1 - prob1)^tuning_c)
  }
  return(allocation_rate)
}

get_efficacious_probabilities <- function(results_list,vaccinees,trial_participants,max_time=10000,tested=F){
  ve_estimate <- c(0.6,1)
  weight_hh_rem <- matrix(0,ncol=2,nrow=length(results_list))
  break_count <- 0
  while(abs(ve_estimate[1]-ve_estimate[2])>0.005&&break_count<5){
    v_count <- 0
    c_count <- 0
    for(iter in 1:length(results_list)){
      results <- results_list[[iter]]
      results <- results[results$DayInfected<=max_time,]
      if(nrow(results)>1){
        weights_out <- get_weighted_results_given_ve(results,ve_point_est=ve_estimate[1],tested)
        weight_hh_rem[iter,] <- weights_out
        v_count <- v_count + sum(results$inTrial==T&results$vaccinated==T)
        c_count <- c_count + sum(results$inTrial==T&results$vaccinated==F)
      }
    }
    ve_estimate[2] <- ve_estimate[1]
    weight_sums <- colSums(weight_hh_rem,na.rm=T)
    pop_sizes2 <- c(sum(vaccinees)-v_count+weight_sums[1], sum(trial_participants) - sum(vaccinees) - c_count+weight_sums[2])
    #print(c(1,des,pop_sizes2,weight_sums[2]>0&&!any(pop_sizes2==0)))
    if(weight_sums[2]>0&&!any(pop_sizes2==0))
      ve_estimate[1] <- calculate_ve(weight_sums,pop_sizes2)
    break_count <- break_count + 1
  }
  return(list(ve_estimate[1],pop_sizes2,weight_sums))
}

get_weight_matrix <- function(infected_nodes,potential_infectees){
  weight_matrix <- matrix(0,nrow=length(infected_nodes),ncol=length(potential_infectees))
  for(i in 1:length(infected_nodes)){
    i_node <- infected_nodes[i]
    i_hr <- c(high_risk_list[[i_node]],household_list[[i_node]])
    i_contact <- contact_list[[i_node]]
    i_nb <- contact_of_contact_list[[i_node]]
    j_hr <- potential_infectees%in%i_hr
    j_contact <- potential_infectees%in%i_contact
    j_nb <- potential_infectees%in%i_nb
    for(j in 1:length(potential_infectees)){
      j_node <- potential_infectees[j]
      if(j_hr[j]){
        weight_matrix[i,j] <- high_risk_scalar
      }else if(j_contact[j]){
        weight_matrix[i,j] <- 1
      }else if(j_nb[j]){
        weight_matrix[i,j] <- neighbour_scalar
      }
    }
  }
  sumzero <- colSums(weight_matrix)>0
  if(sum(sumzero)==0) return(NULL)
  
  # excise nonrelevant nodes
  keep_participants <- potential_infectees[sumzero]
  weight_matrix <- weight_matrix[,sumzero,drop=F]
  return(list(keep_participants,weight_matrix))
}

get_exposures <- function(potential_infectors,removal_days,infectious_days,inc_plus_vacc_shape,inc_plus_vacc_rate,rec_day,weight_matrix){
  #potential_infectors$durations <- removal_days-infectious_days
  end_day <- potential_infectors$RecruitmentDay + 31
  # the total 'infectious force' exerted by each infectious person
  force_of_infection <- force_of_infection_after_0 <- c()
  for(x in 1:nrow(potential_infectors))
    force_of_infection[x] <- sum(pgamma_vector[end_day[x]-infectious_days[x]:removal_days[x]])
  
  get_infectious_force_after_0 <- function(x){
    days <- infectious_days[x]:removal_days[x]
    days <- days[days>rec_day]
    sum(pgamma_vector[end_day[x]-days])
  }
  # the total force after day 0
  for(x in 1:nrow(potential_infectors))
    force_of_infection_after_0[x] <- get_infectious_force_after_0(x)
  # the total force before day 0
  force_of_infection_before_0 <- force_of_infection - force_of_infection_after_0 
  # get exposures for each person
  exposures <- pretrial_exposures <- posttrial_exposures <- c()
  foi <- force_of_infection
  foib0 <- force_of_infection_before_0
  foia0 <- force_of_infection_after_0
  for(i in 1:ncol(weight_matrix)){
    exposures[i] <- sum(weight_matrix[,i] * foi)
    pretrial_exposures[i] <- sum(weight_matrix[,i] * foib0)
    posttrial_exposures[i] <- sum(weight_matrix[,i] * foia0)
  }
  return(list(exposures,pretrial_exposures,posttrial_exposures))
}

get_expected_infectious_exposures <- function(){
  expected_exposure <- expected_pre_exposure <- c()
  # for each person who is infected
  for(x in 1:length(infectee_columns)){
    i <- infectee_columns[x]
    x2 <- which(infected_nodes==keep_participants[i])
    prob_vals <- weight_vals <- vals <- c()
    pre_prob_vals <- pre_weight_vals <- pre_vals <- c()
    # for each person who might have infected them
    for(j  in 1:x2){
      # get all days
      start_day <- infectious_days[j]
      end_day <- min(removal_days[j],infectious_days[x2])
      days <- start_day:end_day
      # split into before and after
      ##!! the 1 is for vax development
      predays <- days[days<=rec_day]
      postdays <- days[days>rec_day]
      # get pre and post probabilities
      if(length(postdays)>0){
        inc_days <- infectious_days[x2]-postdays
        probs <- dgamma(inc_days,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
        weights <- rep(weight_matrix[j,i] ,length(inc_days))
        vals <- c(vals,inc_days)
        weight_vals <- c(weight_vals,weights)
        prob_vals <- c(prob_vals,probs)
      }
      if(length(predays)>0){
        inc_days <- infectious_days[x2]-predays
        probs <- dgamma(inc_days,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
        weights <- rep(weight_matrix[j,i] ,length(inc_days))
        pre_vals <- c(pre_vals,inc_days)
        pre_weight_vals <- c(pre_weight_vals,weights)
        pre_prob_vals <- c(pre_prob_vals,probs)
      }
    }
    expected_day <- sum(weight_vals*prob_vals*vals)/sum(weight_vals*prob_vals)
    expected_exposure[x] <- sum(weight_vals[vals<=expected_day])
    expected_day <- sum(pre_weight_vals*pre_prob_vals*pre_vals)/sum(pre_weight_vals*pre_prob_vals)
    expected_pre_exposure[x] <- sum(pre_weight_vals[pre_vals<=expected_day])
    ##!! if control, use all exposure, and weight 1
    ##!! if vax, store remaining exposure, and label control with remaining weight
    if(is.na(expected_exposure[x])) expected_exposure[x] <- 0
    if(is.na(expected_pre_exposure[x])) expected_pre_exposure[x] <- 0
  }
  return(list(expected_exposure,expected_pre_exposure))
}

summarise_trial <- function(netwk,ve_est_temp=0.7,eval_day=31,pre_randomisation=T){
  results <<- netwk[[1]]
  rec_day <<- results$RecruitmentDay[1]
  results$DayRemoved[is.na(results$DayRemoved)] <- results$RecruitmentDay[is.na(results$DayRemoved)] + eval_day
  potential_infectees <- netwk[[7]]
  if(pre_randomisation){
    potential_infectors <<- results # subset(results,DayRemoved>RecruitmentDay)
  }else{
    ##!! +1 for vaccination time
    potential_infectors <<- results[results$DayRemoved>results$RecruitmentDay,]
    potential_infectees <- potential_infectees[!potential_infectees%in%subset(results,DayRemoved<=RecruitmentDay)$InfectedNode]
  }
  
  trial_nodes <- NULL
  if(nrow(potential_infectors)>0){
    # get contact matrix weights
    infected_nodes <<- potential_infectors$InfectedNode
    infectious_days <<- potential_infectors$DayInfectious
    removal_days <<- potential_infectors$DayRemoved
    gwm <- get_weight_matrix(infected_nodes,potential_infectees)
    keep_participants <<- gwm[[1]]
    if(length(keep_participants)>0){
      weight_matrix <<- gwm[[2]]
      
      all_exposures <- get_exposures(potential_infectors,removal_days,infectious_days,inc_plus_vacc_shape,inc_plus_vacc_rate,rec_day,weight_matrix)
      exposures <- all_exposures[[1]]
      pretrial_exposures <- all_exposures[[2]]
      posttrial_exposures <- all_exposures[[3]]
      all_exposures <- c()
      
      # initialise trial nodes df
      lp <- length(keep_participants)
      outcome <- keep_participants%in%infected_nodes
      vaccinated <- keep_participants%in%netwk[[6]]
      trial_nodes <- as.data.frame(cbind(node=keep_participants,weight=rep(1,lp),total=exposures,pretrial=pretrial_exposures,
                                         posttrial=posttrial_exposures,time=rep(0,lp)),check.names=F,fix.empty.names=F)
      # add information
      trial_nodes$outcome <- outcome
      trial_nodes$vaccinated <- vaccinated
      # reorder to put vaccinees and infectees at the end
      #trial_nodes[order(outcome,vaccinated)]
      trial_nodes <- trial_nodes[c(which(!vaccinated&!outcome),which(vaccinated&!outcome),which(outcome)) ,]
      outcome <- trial_nodes$outcome 
      vaccinated <- trial_nodes$vaccinated
      # if using only post-randomisation time
      if(!pre_randomisation){
        trial_nodes$time[!outcome] <- trial_nodes$posttrial[!outcome]
      }else{
        # not infected, not vaccinated people have total as their exposure time
        trial_nodes$time[!outcome&!vaccinated] <- trial_nodes$total[!outcome&!vaccinated]
        # vaccinated, not infected people get included twice, once with time=pretrial, not vaccinated, then with time=posttrial, vaccinated
        vax_nodes <- trial_nodes[!outcome&vaccinated,]
        trial_nodes$time[!outcome&vaccinated] <- trial_nodes$posttrial[!outcome&vaccinated]
        if(nrow(vax_nodes)>0){
          vax_nodes$vaccinated <- F
          vax_nodes$time <- vax_nodes$pretrial
        }
      }
      ## do inf nodes separately
      inf_nodes <- trial_nodes[outcome,]
      ##
      
      ## exposures for infectees
      infectee_columns <<- which(keep_participants%in%infected_nodes)
      if(sum(infectee_columns)>0){
        expected_infectious_exposures <- get_expected_infectious_exposures()
        expected_exposure <- expected_infectious_exposures[[1]]
        expected_pre_exposure <- expected_infectious_exposures[[2]]
        
        # there is one entry for the control: the expected weight survived, which is expected_exposure+expected_pre_exposure and weight=1
        # for vax, there are three entries: 
        # weight=prob infected before vax, exposure=expected time exposed before (not vax, inf)
        # weight=prob infected after vax, exposure=expected time exposed after (vax, inf)
        # weight=prob infected after vax, exposure=time exposed before (not vax, not inf)
        
        ## weights for infectees
        infectee_weights <- get_infectee_weights(results,ve_est_temp)
        prob_after_0 <- rowSums(infectee_weights[[1]])
        ##!! weighting only for vaccinated
        inf_trial_nodes <- inf_nodes$node[inf_nodes$node%in%results$InfectedNode]
        inf_nodes$weight[inf_nodes$node%in%results$InfectedNode] <- prob_after_0[match(inf_trial_nodes,infectee_weights[[2]])]
        inf_nodes$weight[is.na(inf_nodes$weight)] <- 0
        inf_nodes$time[match(keep_participants[infectee_columns],inf_nodes$node)] <- expected_exposure
        trial_nodes <- bind_rows(trial_nodes,inf_nodes)
        
        if(pre_randomisation){
          ##!! duplicate and add in control
          inf_pre_zero <- survived_pre_zero <- inf_nodes
          inf_pre_zero$vaccinated <- survived_pre_zero$vaccinated <- F
          survived_pre_zero$outcome <- F
          inf_pre_zero$weight[inf_pre_zero$node%in%results$InfectedNode] <- 1 - inf_nodes$weight[inf_nodes$node%in%results$InfectedNode]
          inf_pre_zero$time[match(keep_participants[infectee_columns],inf_pre_zero$node)] <- expected_pre_exposure
          survived_pre_zero$time <- survived_pre_zero$pretrial
          #survived_pre_zero <- subset(survived_pre_zero,node%in%netwk[[6]])
          #inf_pre_zero$node <- -inf_pre_zero$node
          #survived_pre_zero$node <- survived_pre_zero$node/1000
          #infected_nodes <- c(infected_nodes,inf_pre_zero$node)
          trial_nodes <- bind_rows(trial_nodes,vax_nodes,inf_pre_zero,survived_pre_zero)
        }
      }
      trial_nodes$weight[is.na(trial_nodes$weight)] <- 0
      trial_nodes <- trial_nodes[trial_nodes$weight>0&trial_nodes$time>0,]
      if(nrow(trial_nodes)==0) return(NULL)
      
    }
  }
  return(trial_nodes)
}

