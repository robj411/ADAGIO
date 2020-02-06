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
  # the day the cluster is recruited
  recruit_day <- results$RecruitmentDay[1]
  # the day individuals became infectious
  days_infectious <- results$DayInfectious
  # the subset of those that might have infected others
  infectors <- days_infectious[1:(nrow(results)-1)]
  # the durations for which they were infectious
  infector_durations <- results$DayRemoved[1:(nrow(results)-1)] - infectors
  infector_durations[is.na(infector_durations)|infector_durations>21] <- 21
  weight_hh_rem <- c(0,0)
  # those who were infected by someone else
  infectee_index <- days_infectious>recruit_day
  if(sum(infectee_index)>0){
    infectees <- days_infectious[infectee_index]
    infector_names <- results$InfectedNode[1:(nrow(results)-1)]
    infectee_names <- results$InfectedNode[infectee_index]
    infectee_trial <- results$inTrial[infectee_index]
    infectee_vaccinated <- results$vaccinated[infectee_index]
    for(j in 1:length(infectees)){
      if(infectee_trial[j]){
        # which infectors might have infected infectee j
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
          weight_hh_rem[1] <- weight_hh_rem[1] + sum(prob_after_0*normalised_prob_infectors)
        }else{
          weight_hh_rem[2] <- weight_hh_rem[2] + sum(prob_after_0*normalised_prob_infectors)
        }
      }
    }
  }
  return(weight_hh_rem)
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
  while(abs(ve_estimate[1]-ve_estimate[2])>0.01){
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
    if(!any(weight_sums==0))
      ve_estimate[1] <- calculate_ve(weight_sums,pop_sizes2)
  }
  return(list(ve_estimate[1],pop_sizes2,weight_sums))
}

get_efficacious_probabilities2 <- function(results_list,vaccinees,trial_participants,max_time=10000){
  ve_estimate <- c(0.6,1)
  weight_hh_rem <- matrix(0,ncol=2,nrow=length(results_list))
  while(abs(ve_estimate[1]-ve_estimate[2])>0.01){
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
    if(weight_sums[2]>0&&!any(pop_sizes2==0))
      ve_estimate[1] <- calculate_ve(weight_sums,pop_sizes2)
  }
  return(list(ve_estimate[1],pop_sizes2,weight_sums))
}
