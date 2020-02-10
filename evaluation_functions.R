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
    prob_infectors <- sapply(1:length(rows),function(x)probability_by_lag_given_removal[[max(infector_durations[x]-1,1)]][rows[x],cols[x]])
    prob_after_0 <- sapply(1:length(rows),function(x)probability_after_day_0_given_removal[[max(infector_durations[x]-1,1)]][rows[x],cols[x]])
    normalised_prob_infectors <- prob_infectors*hh_weight/sum(prob_infectors*hh_weight)
    weight_hh_rem <- weight_hh_rem + sum(prob_after_0*normalised_prob_infectors)
  }
  return(list(c(true_positives,weight,binary,weight_hh,weight_hh_rem),how_surprising))
}

get_weighted_results_given_ve <- function(results,ve_point_est){
  weight_hh_rem <- colSums(get_infectee_weights(results,ve_point_est)[[1]])
  return(weight_hh_rem)
}

get_infectee_weights <- function(results,ve_point_est){
  # the day the cluster is recruited
  recruit_day <- results$RecruitmentDay[1]
  # the day individuals became infectious
  days_infectious <- results$DayInfectious
  # the subset of those that might have infected others
  infectors <- days_infectious[1:(nrow(results)-1)]
  # the durations for which they were infectious
  infector_durations <- results$DayRemoved[1:(nrow(results)-1)] - infectors
  infector_durations[is.na(infector_durations)|infector_durations>21] <- 21
  weight_hh_rem <- matrix(0,ncol=2,nrow=1)
  infectee_names <- c()
  # those who were infected by someone else
  infectee_index <- days_infectious>recruit_day
  if(sum(infectee_index)>0){
    weight_hh_rem <- matrix(0,ncol=2,nrow=sum(infectee_index))
    infectees <- days_infectious[infectee_index]
    infector_names <- results$InfectedNode[1:(nrow(results)-1)]
    infectee_names <- results$InfectedNode[infectee_index]
    infectee_trial <- results$inTrial[infectee_index]
    infectee_vaccinated <- results$vaccinated[infectee_index]
    for(j in 1:length(infectees)){
      if(infectee_trial[j]){
        # which infectors might have infected infectee j
        infectors_for_j <- infectors<infectees[j]
        # rows: the time lag between infector and infectee becoming infectious
        rows <- pmin(infectees[j]-infectors[infectors_for_j],nrow(probability_by_lag))
        # cols: the day the potential infector became infectious relative to recruitment day
        cols <- ref_recruit_day-recruit_day+infectors[infectors_for_j]
        # subtract the vaccine incubation period (increase the reference day from 0)
        cols <- pmax(cols - ceiling(qtruncnorm(0.5,a=0,vacc_mean,vacc_sd)),1)
        # weights for all relationships given known network
        hh_weight <- sapply(infector_names[infectors_for_j],
                            function(x)as.numeric(x%in%contact_list[[infectee_names[j]]]) +
                              as.numeric(x%in%household_list[[infectee_names[j]]])*(high_risk_scalar-1) +
                              as.numeric(x%in%high_risk_list[[infectee_names[j]]])*(high_risk_scalar-1) +
                              as.numeric(x%in%contact_of_contact_list[[infectee_names[j]]])*neighbour_scalar
        )
        #as.numeric(sapply(infector_names[infectors<infectees[j]],function(x)x%in%household_list[[infectee_names[j]]]))*(high_risk_scalar-1)+1
        ##!! using contact and removal information
        # probabilities for infectors to infect infectee j
        prob_infectors <- sapply(1:length(rows),function(x)probability_by_lag_given_removal[[max(infector_durations[x]-1,1)]][rows[x],cols[x]])
        # probabilities infected after recruitment day given infected by infector
        prob_after_0 <- sapply(1:length(rows),function(x)probability_after_day_0_given_removal[[max(infector_durations[x]-1,1)]][rows[x],cols[x]])
        # store complement
        yij <- 1 - prob_after_0
        # if vaccinated, adjust probability to be infected after day 0
        if(infectee_vaccinated[j]) prob_after_0 <- (1-ve_point_est)*prob_after_0/(yij+(1-ve_point_est)*prob_after_0)
        # recalculate probabilities for infectors, which will be the same for non-vaccinated
        prob_infectors <- prob_infectors*(yij+prob_after_0)
        # calculate normalised infector probabilities
        normalised_prob_infectors <- prob_infectors*hh_weight/sum(prob_infectors*hh_weight)
        # add to weight for vaccinated or unvaccinated
        if(infectee_vaccinated[j]){
          weight_hh_rem[j,1] <- sum(prob_after_0*normalised_prob_infectors)
        }else{
          weight_hh_rem[j,2] <- sum(prob_after_0*normalised_prob_infectors)
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

get_efficacious_probabilities <- function(results_list,vaccinees,trial_participants,max_time=10000){
  ve_estimate <- c(0.6,1)
  weight_hh_rem <- matrix(0,ncol=2,nrow=length(results_list))
  while(abs(ve_estimate[1]-ve_estimate[2])>0.005){
    v_count <- 0
    c_count <- 0
    for(iter in 1:length(results_list)){
      results <- results_list[[iter]]
      results <- results[results$DayInfected<=max_time,]
      if(nrow(results)>1){
        weights_out <- get_weighted_results_given_ve(results,ve_point_est=ve_estimate[1])
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
  }
  return(list(ve_estimate[1],pop_sizes2,weight_sums))
}

get_efficacious_probabilities2 <- function(results_list,vaccinees,trial_participants,max_time=10000){
  ve_estimate <- c(0.6,1)
  weight_hh_rem <- matrix(0,ncol=2,nrow=length(results_list))
  while(abs(ve_estimate[1]-ve_estimate[2])>0.005){
    #v_count <- 0
    #c_count <- 0
    for(iter in 1:length(results_list)){
      results <- results_list[[iter]]
      results <- results[results$DayInfected<=max_time,]
      if(nrow(results)>1){
        weights_out <- get_weighted_results_given_ve(results,ve_point_est=ve_estimate[1])
        weight_hh_rem[iter,] <- weights_out
        #v_count <- v_count + sum(results$inTrial==T&results$vaccinated==T)
        #c_count <- c_count + sum(results$inTrial==T&results$vaccinated==F)
      }
    }
    ve_estimate[2] <- ve_estimate[1]
    weight_sums <- colSums(weight_hh_rem,na.rm=T)
    pop_sizes2 <- c(sum(vaccinees), sum(trial_participants) - sum(vaccinees))
    #print(c(2,des,pop_sizes2,weight_sums[2]>0&&!any(pop_sizes2==0)))
    if(weight_sums[2]>0&&!any(pop_sizes2==0))
      ve_estimate[1] <- calculate_ve(weight_sums,pop_sizes2)
  }
  return(list(ve_estimate[1],pop_sizes2,weight_sums))
}

summarise_trial <- function(netwk,ve_est_temp=0.7){
  results <- netwk[[1]]
  rec_day <- results$RecruitmentDay[1]
  results$DayRemoved[is.na(results$DayRemoved)] <- results$RecruitmentDay[is.na(results$DayRemoved)] + eval_day
  potential_infectors <- subset(results,DayRemoved>RecruitmentDay)
  
  trial_nodes <- NULL
  if(nrow(potential_infectors)>0){
    weight_matrix <- matrix(0,nrow=nrow(potential_infectors),ncol=length(netwk[[7]]))
    for(j in 1:length(netwk[[7]])){
      j_node <- netwk[[7]][j]
      j_hr <- c(high_risk_list[[j_node]],household_list[[j_node]])
      j_contact <- contact_list[[j_node]]
      j_nb <- contact_of_contact_list[[j_node]]
      for(i in 1:nrow(potential_infectors)){
        i_node <- potential_infectors$InfectedNode[i]
        if(i_node%in%j_hr){
          weight_matrix[i,j] <- high_risk_scalar
        }else if(i_node%in%j_contact){
          weight_matrix[i,j] <- 1
        }else if(i_node%in%j_nb){
          weight_matrix[i,j] <- neighbour_scalar
        }
      }
    }
    
    sumzero <- colSums(weight_matrix)>0
    if(sum(sumzero)==0) return(NULL)
    
    keep_participants <- netwk[[7]][sumzero]
    weight_matrix <- weight_matrix[,sumzero,drop=F]
    
    potential_infectors$durations <- potential_infectors$DayRemoved-potential_infectors$DayInfectious
    potential_infectors$end_day <- potential_infectors$RecruitmentDay + 31
    potential_infectors$force_of_infection <- sapply(1:nrow(potential_infectors),function(x)
      sum(pgamma(potential_infectors$end_day[x]-potential_infectors$DayInfectious[x]:potential_infectors$DayRemoved[x],shape=incperiod_shape,rate=incperiod_rate))
    )
    
    exposures <- apply(weight_matrix,2,function(x) sum(x * potential_infectors$force_of_infection))
    trial_nodes <- data.frame(node=keep_participants,time=exposures,weight=1)
    
    ##
    
    ## exposures for infectees
    infectee_columns <- which(keep_participants%in%potential_infectors$InfectedNode)
    if(sum(infectee_columns)>0){
      expected_exposure <- c()
      for(x in 1:length(infectee_columns)){
        i = infectee_columns[x]
        x2 <- which(potential_infectors$InfectedNode==keep_participants[i])
        prob_vals <- weight_vals <- vals <- c()
        for(j  in 1:x2){
          start_day <- max(potential_infectors$DayInfectious[j],rec_day+1)
          end_day <- min(potential_infectors$DayRemoved[j],potential_infectors$DayInfectious[x2])
          days <- potential_infectors$DayInfectious[x2]-start_day:end_day
          probs <- dgamma(days,shape=incperiod_shape,rate=incperiod_rate)
          weights <- rep(weight_matrix[j,i] ,length(days))
          vals <- c(vals,days)
          weight_vals <- c(weight_vals,weights)
          prob_vals <- c(prob_vals,probs)
        }
        expected_day <- sum(weight_vals*prob_vals*vals)/sum(weight_vals*prob_vals)
        expected_exposure[x] <- sum(weight_vals[vals<=expected_day])
        if(is.na(expected_exposure[x])) expected_exposure[x] <- 0
      }
      
      ## weights for infectees
      infectee_weights <- get_infectee_weights(results,ve_est_temp)
      prob_after_0 <- rowSums(infectee_weights[[1]])
      
      inf_trial_nodes <- trial_nodes$node[trial_nodes$node%in%results$InfectedNode]
      trial_nodes$weight[trial_nodes$node%in%results$InfectedNode] <- prob_after_0[match(inf_trial_nodes,infectee_weights[[2]])]
      trial_nodes$time[match(keep_participants[infectee_columns],trial_nodes$node)] <- expected_exposure
    }
    trial_nodes$weight[is.na(trial_nodes$weight)] <- 0
    trial_nodes <- subset(trial_nodes,weight>0)
    
    trial_nodes$outcome <- trial_nodes$node%in%potential_infectors$InfectedNode
    trial_nodes$vaccinated <- trial_nodes$node%in%netwk[[6]]
    if(nrow(trial_nodes)==0) return(NULL)
  }
  return(trial_nodes)
}

iterate_ph_model <- function(netwk_list){
  ves <- c(0.6,1)
  while(abs(ves[1]-ves[2])>0.005){
    trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=ves[1])
    
    tte <- do.call(rbind,lapply(1:length(trial_summary),function(cluster)if(!is.null(trial_summary[[cluster]]))cbind(trial_summary[[cluster]],cluster)))
    if(cluster_flag==1){
      survmodel <- coxme(Surv(time,outcome)~vaccinated+(1|cluster),weights=weight,tte)
      vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(vcov(survmodel))))
      zval <- survmodel$coefficient/sqrt(vcov(survmodel))
    }else{
      survmodel <- coxph(Surv(time,outcome)~vaccinated,weights=weight,tte)
      vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*as.vector(sqrt(survmodel$var)))
      zval <- survmodel$coefficient/sqrt(survmodel$var)
    }
    pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
    #print(c(vaccEffEst,zval,pval))
    ves[2] <- ves[1]
    if(!is.na(vaccEffEst[1]))
      ves[1] <- vaccEffEst[1]
  }
  return(c(pval,ves[1]))
}
