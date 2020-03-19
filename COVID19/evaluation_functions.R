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
  weight_hh_rem <- colSums(get_infectee_weights(results,ve_point_est,tested))
  return(weight_hh_rem)
}

get_infectee_weights <- function(results,ve_point_est,tested=F){
  
  # the day the cluster is recruited
  recruit_day <- results$RecruitmentDay
  # the day individuals became infectious
  days_infectious <- results$DayInfectious
  weight_hh_rem <- c()
  # those who were infected by someone else
  infectee_index <- !is.na(recruit_day) & days_infectious>recruit_day
  if(sum(infectee_index)>0){
    infectees <- days_infectious[infectee_index] - recruit_day[infectee_index]
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
        weight_hh_rem[j] <- prob_after_0
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

response_adapt <- function(fails,pop_sizes2,max_time, adaptation='TST'){
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
    j <- max_time # t - trial_startday
    bigT <- nClusters # trial_length
    tuning_c <- ifelse(adaptation=='TS',1,(j/bigT))
    #print(tuning_c)
    p0 <- rbeta(1000,1+successes[2],1+fails[2])
    p1 <- rbeta(1000,1+successes[1],1+fails[1])
    prob1 <- sum(p1>p0)/1000
    allocation_rate <- prob1^tuning_c / (prob1^tuning_c + (1 - prob1)^tuning_c)
  }
  if(allocation_rate==0) allocation_rate <- 1e-3
  if(allocation_rate==1) allocation_rate <- 1-1e-3
  return(allocation_rate)
}

get_weights_from_all_results <- function(all_results){
  v_count <- sum(all_results$weight[all_results$infected==T&all_results$vaccinated==T])
  c_count <- sum(all_results$weight[all_results$infected==T&all_results$vaccinated==F])
  fails <- c(v_count,c_count)
  pop_sizes2 <- c(sum(all_results$weight[all_results$vaccinated==T]), 
                  sum(all_results$weight[all_results$vaccinated==F]))
  list(fails,pop_sizes2)
}

get_efficacious_probabilities <- function(results_list,vaccinees,trial_participants,
                                          max_time=10000,tested=F,randomisation_ratios=NULL,rbht_norm=0,people_per_ratio=NULL,adaptation='TST'){
  controls <- trial_participants - vaccinees
  if(is.null(randomisation_ratios)) randomisation_ratios <- rep(0.5,length(trial_participants))
  
  result_tab <- do.call(rbind,lapply(1:length(results_list),function(x){
    results <- results_list[[x]]
    y <- subset(results,!is.na(RecruitmentDay))
    ##!! could include also RecruitmentDay
    w <- subset(y,DayInfected<max_time)
    z <- subset(w,RecruitmentDay<DayInfectious)
    if(nrow(z)>0) {
      z$startDay <- x
      z$allocRatio <- randomisation_ratios[x]
      z$infected <- T
    }
    z
  }))
  
  uninf_vacc <- vaccinees - sapply(results_list,function(x)sum(x$vaccinated))
  uninf_cont <- trial_participants - vaccinees - sapply(results_list,function(x)sum(x$inTrial&!x$vaccinated))
  
  uninf <- data.frame(vaccinated=c(rep(T,sum(uninf_vacc)),rep(F,sum(uninf_cont))),
                      allocRatio=c(rep(randomisation_ratios,uninf_vacc),rep(randomisation_ratios,uninf_cont)),
                      weight=1,
                      infected=F)

  
  ve_estimate <- c(0.6,1)
  break_count <- 0
  while(abs(ve_estimate[1]-ve_estimate[2])>0.005&&break_count<5){
    result_tab$weight <- get_infectee_weights(result_tab,ve_estimate[1],tested)
    ve_estimate[2] <- ve_estimate[1]
    all_results <- rbind(result_tab[,match(colnames(uninf),colnames(result_tab))],uninf)
    if(rbht_norm==1)
      all_results$weight <- all_results$weight / (all_results$vaccinated + (-1) ^ all_results$vaccinated * all_results$allocRatio)
    if(rbht_norm<2){
      weights <- get_weights_from_all_results(all_results)
      fails <- weights[[1]]
      pop_sizes2 <- weights[[2]]
      if(fails[2]>0&&!any(pop_sizes2==0))
        ve_estimate[1] <- calculate_ve(fails,pop_sizes2)
    }else{
      excluded_people <- sapply(people_per_ratio[,2],function(p) sum(sapply(1:p,function(x){
        results <- results_list[[x]]
        y <- subset(results,!is.na(RecruitmentDay))
        w <- subset(y,DayInfected<max_time)
        sum(w$DayInfectious<=w$RecruitmentDay)})))
      people_per_ratio[,1] <- people_per_ratio[,1] - excluded_people
      M <- 1000
      new_ve <- 0
      all_results_original <- rbind(result_tab[,match(colnames(uninf),colnames(result_tab))],uninf)
      for(i in 1:M){
        #vacc_half <- round(people_per_ratio[1,1]/2)
        #first_sample <- c(sample(which(all_results_original$vaccinated==T),vacc_half,replace=F),
        #                  sample(which(all_results_original$vaccinated==F),people_per_ratio[1,1]-vacc_half,replace=F))
        first_sample <- sample(nrow(all_results_original),people_per_ratio[1,1],replace=F)
        not_sampled <- c(1:nrow(all_results_original))[-first_sample]
        all_results <- all_results_original[first_sample,]
        all_results$allocRatio <- 0.5
        #all_results$weight <- all_results$weight / (all_results$vaccinated + (-1) ^ all_results$vaccinated * all_results$allocRatio)
        weights <- get_weights_from_all_results(all_results)
        allocation_ratio <- response_adapt(weights[[1]],weights[[2]],max_time=people_per_ratio[1,2], adaptation)
        for(j in 2:(nrow(people_per_ratio)+1)){
          max_people <- nrow(all_results_original)
          max_t <- max_time
          if(j<=nrow(people_per_ratio)){
            max_people <- people_per_ratio[j,1]
            max_t <- people_per_ratio[j,2]
          }
          all_results_temp <- sample(not_sampled,max_people-people_per_ratio[j-1,1],replace=F)
          not_sampled <- not_sampled[!not_sampled%in%all_results_temp]
          temp_results <- all_results_original[all_results_temp,]
          temp_results$allocRatio <- allocation_ratio
          #temp_results$weight <- temp_results$weight / (temp_results$vaccinated + (-1) ^ temp_results$vaccinated * temp_results$allocRatio)
          all_results <- rbind(temp_results,all_results)
          weights <- get_weights_from_all_results(all_results)
          allocation_ratio <- response_adapt(weights[[1]],weights[[2]],max_time=max_t, adaptation)
        }
        all_results$weight <- all_results$weight / (all_results$vaccinated + (-1) ^ all_results$vaccinated * all_results$allocRatio)
        weights <- get_weights_from_all_results(all_results)
        fails <- weights[[1]]
        pop_sizes2 <- weights[[2]]
        if(fails[2]>0&&!any(pop_sizes2==0))
          new_ve[i] <- calculate_ve(fails,pop_sizes2)
      }
      ve_estimate[1] <- mean(new_ve,na.rm=T)
    }
    
    break_count <- break_count + 1
  }
  return(list(ve_estimate[1],pop_sizes2,fails))
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

