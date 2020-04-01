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

get_weighted_results_given_ve <- function(results,ve_point_est,contact_network=2,tested=F){
  weight_hh_rem <- colSums(get_infectee_weights(results,ve_point_est,contact_network,tested)[[1]])
  return(weight_hh_rem)
}

get_infectee_weights <- function(results,ve_point_est,contact_network=2,tested=F,correct_for_ve=T){
  
  if(contact_network==2){
    get_contact_weight <- function(x,j){
      as.numeric(x%in%contact_list[[infectee_names[j]]]) +
        as.numeric(x%in%household_list[[infectee_names[j]]])*(high_risk_scalar-1) +
        as.numeric(x%in%inv_hr_list[[infectee_names[j]]])*(high_risk_scalar-1) +
        as.numeric(x%in%contact_of_contact_list[[infectee_names[j]]])*neighbour_scalar
    }
  }else if(contact_network==1){
    get_contact_weight <- function(x,j){
      1+as.numeric(x%in%household_list[[infectee_names[j]]])
    }
  }else if(contact_network==0){
    get_contact_weight <- function(x,j){
      1
    }
  }
  
  # the day the cluster is recruited
  recruit_day <- results$RecruitmentDay
  # the day individuals became infectious
  days_infectious <- results$DayInfectious
  # the subset of those that might have infected others
  infectors <- days_infectious[1:(nrow(results)-1)]
  # the durations for which they were infectious
  infector_durations <- results$DayRemoved[1:(nrow(results)-1)] - infectors
  infector_durations[is.na(infector_durations)|infector_durations>21] <- 21
  ## initialise to return empty
  weight_hh_rem <- matrix(0,ncol=2,nrow=1)
  infectee_names <- c()
  # those who were infected by someone else
  infectee_index <- !is.na(recruit_day) & days_infectious>recruit_day
  if(sum(infectee_index)>0){
    weight_hh_rem <- matrix(0,ncol=2,nrow=sum(infectee_index))
    infectees <- days_infectious[infectee_index]
    infector_names <- results$InfectedNode[1:(nrow(results)-1)]
    infectee_names <- results$InfectedNode[infectee_index]
    infectee_trial <- results$inTrial[infectee_index]
    infectee_vaccinated <- results$vaccinated[infectee_index]
    for(j in 1:length(infectees)){
      if(infectee_trial[j]){
        if(contact_network==-1){
          prob_after_0 <- pgamma(c(days_infectious[infectee_index] - recruit_day[infectee_index])[j]-1,shape=inc_plus_vacc_shape,rate=inc_plus_vacc_rate)
        }else{
          # recruitment day for infectee j
          recj <- c(recruit_day[infectee_index])[j]
          # which infectors might have infected infectee j
          infectors_for_j <- infectors<infectees[j]
          # rows: the time lag between infector and infectee becoming infectious
          rows <- pmin(infectees[j]-infectors[infectors_for_j],nrow(probability_by_lag))
          infector_durations_for_j <- infector_durations[infectors_for_j]
          # cols: the day the potential infector became infectious relative to recruitment day
          ##!! assumes all recruited on the same day
          cols <- ref_recruit_day-recj+infectors[infectors_for_j]
          ##!! in case recruitment day exceeds 30
          cols <- pmax(cols,1)
          # weights for all relationships given known network
          hh_weight <- sapply(infector_names[infectors_for_j],get_contact_weight,j)
          ##!! using contact and removal information
          # probabilities for infectors to infect infectee j
          prob_infectors <- sapply(1:length(rows),function(x)probability_by_lag_given_removal_mat[rows[x],max(infector_durations_for_j[x]-1,1)])
          # probabilities infected after recruitment day given infected by infector
          prob_after_0 <- sapply(1:length(rows),function(x)probability_after_day_0_given_removal[[max(infector_durations_for_j[x]-1,1)]][rows[x],cols[x]])
          # adjust prob_infectors for vaccinee
          if(infectee_vaccinated[j]&correct_for_ve){
            # rows: the time lag between infector and infectee becoming infectious
            ## split into two: before day 0 and after
            before <- rep(0,sum(infectors_for_j))
            after <- rep(0,sum(infectors_for_j))
            for(k in 1:sum(infectors_for_j)){
              if(infectors[k]+infector_durations_for_j[k]<=recj+1){
                before[k] <- prob_infectors[k]
              }else if(infectors[k]>recj){
                after[k] <- prob_infectors[k]
              }else{
                before_dur <- infectors[k] + infector_durations_for_j[k] - recj
                after_dur <- infector_durations_for_j[k] - before_dur
                # rows: the time lag between infector and infectee becoming infectious
                # cols: duration infectious
                before[k] <- probability_by_lag_given_removal_mat[infectees[j]-infectors[k],max(before_dur,1)]
                after[k] <- probability_by_lag_given_removal_mat[infectees[j]-recj+1,max(after_dur-1,1)]    
              }
            }
            cf <- 1 - after
            after <- (1-ve_point_est)*after/(cf+(1-ve_point_est)*after+1e-16)
            prob_infectors <- before + after
          }
        }
        if(tested) {
          # if test positive, probability=0
          if(c(results$DayInfected[infectee_index])[j]<c(results$RecruitmentDay[infectee_index])[j]){
            prob_after_0 <- 0
          }else{
            numerator <- prob_after_0
            denom <- pgamma(c(days_infectious[infectee_index] - recruit_day[infectee_index])[j]-1,shape=incperiod_shape,rate=incperiod_rate)
            prob_after_0 <- numerator/denom
            ##!! can we set to 1 for controls?
            if(!infectee_vaccinated[j])
              prob_after_0 <- 1
          }
        }
        # store complement
        yij <- 1 - prob_after_0
        # if vaccinated, adjust probability to be infected after day 0
        if(infectee_vaccinated[j]&correct_for_ve) prob_after_0 <- (1-ve_point_est)*prob_after_0/(yij+(1-ve_point_est)*prob_after_0+1e-16)
        if(contact_network>-1){
          # recalculate probabilities for infectors, which will be the same for non-vaccinated
          prob_infectors <- prob_infectors*(yij+prob_after_0)
          # calculate normalised infector probabilities
          normalised_prob_infectors <- prob_infectors*hh_weight/sum(prob_infectors*hh_weight)
        }else{
          normalised_prob_infectors <- 1
        }
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

get_weights_from_all_results <- function(all_results){
  weight_vector <- all_results$weight
  v_index <- all_results$vaccinated==1
  i_index <- all_results$infected
  all_weight <- sum(weight_vector)
  vax_weight <- sum(weight_vector[v_index])
  v_count <- sum(weight_vector[v_index&i_index])
  c_count <- sum(weight_vector[!v_index&i_index])
  fails <- c(v_count,c_count)
  pop_sizes2 <- c(vax_weight, 
                  all_weight-vax_weight)
  list(fails,pop_sizes2)
}

response_adapt <- function(fails,pop_sizes2,days=31, adaptation='TST'){
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
    j <- days # t - trial_startday
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

trend_robust_function <- function(results_list,vaccinees,trial_participants,contact_network=0,
                                  tested=F,randomisation_ratios=NULL,adaptation='TST',people_per_ratio){
  
  
  result_lst <- lapply(1:length(results_list),function(x){
    results <- results_list[[x]]
    y <- subset(results,!is.na(RecruitmentDay))
    ##!! could include also RecruitmentDay
    #w <- subset(y,DayInfected<max_time)
    z <- subset(y,RecruitmentDay<DayInfectious)
    if(nrow(z)>0) {
      z$startDay <- x
      z$allocRatio <- randomisation_ratios[x]
      z$infected <- T
    }
    z
  })
  
  excluded <- sapply(1:0,function(x){
    sapply(results_list,function(y)
      sum(y$inTrial==T&y$RecruitmentDay>=y$DayInfectious&y$vaccinated==x)
  )})
  
  controls <- trial_participants - vaccinees
  true_trial_participants <- trial_participants - rowSums(excluded)
  true_controls <- controls - excluded[,2]
  true_vax <- vaccinees - excluded[,1]
  
  result_tab_list <- lapply(1:length(result_lst),function(x){
    y <- result_lst[[x]]#[-1,]
    if(nrow(y)>0){
      y$weight <- 0
      weightings <- get_infectee_weights(results_list[[x]],0,contact_network,tested)
      y$weight[match(weightings[[2]],y$InfectedNode)] <- rowSums(weightings[[1]])
      y <- subset(y,weight>0)
      #y <- y[,match(colnames(uninf_list[[x]]),colnames(y))]
    }
    y#rbind(y,uninf_list[[x]])
  })
  
  result_tab <- do.call(rbind,result_tab_list)
  
  M <- 500
  pval <- c()
  unique_ratios <- c(0.5,people_per_ratio[-nrow(people_per_ratio),3])
  all_results_original <- result_tab#rbind(result_tab[,match(colnames(uninf),colnames(result_tab))],uninf)
  set_indices <- lapply(1:length(unique_ratios),function(x)which(all_results_original$allocRatio==unique_ratios[x]))
  indices <- lapply(1:length(unique_ratios),function(x)which(all_results_original$allocRatio%in%unique_ratios[1:x]))
  day <- people_per_ratio[,2]
  cases_per_ratio <- c(sapply(day,function(x)sum(sapply(result_tab_list[1:x],nrow))),nrow(result_tab))
  noncases_per_ratio <- c(people_per_ratio[,1],sum(true_trial_participants)) - cases_per_ratio
  first_results <- all_results_original[indices[[1]],]#head(all_results_original,last_index[1])#
  for(i in 1:M){
    first_allocations <- rbinom(nrow(first_results),1,0.5)
    all_results_original$vaccinated[1:cases_per_ratio[1]] <- first_allocations
    first_results$vaccinated <- first_allocations
    vax <- rbinom(1,noncases_per_ratio[1],0.5)
    ve_estimate <- fast_efficacy(first_results,vax,people_per_ratio[1,1])[[1]]
    vax_weights <- first_results$weight[first_allocations==1]
    cf <- 1 - vax_weights
    first_results$weight[first_allocations==1] <- (1-ve_estimate)*vax_weights/(cf+(1-ve_estimate)*vax_weights+1e-16)
    fails <- c(sum(first_results$weight[first_results$vaccinated==1]),sum(first_results$weight[first_results$vaccinated==0]))
    popsizes <- c(vax,noncases_per_ratio[1]-vax)
    #weights <- get_weights_from_all_results(first_results)
    allocation_ratio <- response_adapt(fails,popsizes,days=day[1], adaptation=adaptation)
    for(j in 2:length(set_indices)){
      #temp_results_index <- all_results_original$allocRatio==unique_ratios[j]
      allocations <- rbinom(length(set_indices[[j]]),1,allocation_ratio)
      all_results_original$vaccinated[set_indices[[j]]] <- allocations
      if(j<length(indices)) {
        all_results <- all_results_original[indices[[j]],]#head(all_results_original,last_index[j])#
        npart <- people_per_ratio[j,1]
      }else{
        all_results <- all_results_original
        npart <- sum(true_trial_participants)
      }
      vax <- vax + rbinom(1,noncases_per_ratio[j]-noncases_per_ratio[j-1],allocation_ratio)
      ve_estimate <- fast_efficacy(all_results[1:cases_per_ratio[j],],vax,npart)[[1]]
      vax_weights <- all_results$weight[all_results$vaccinated==1]
      cf <- 1 - vax_weights
      all_results$weight[all_results$vaccinated==1] <- (1-ve_estimate)*vax_weights/(cf+(1-ve_estimate)*vax_weights+1e-16)
      fails <- c(sum(all_results$weight[all_results$vaccinated==1]),sum(all_results$weight[all_results$vaccinated==0]))
      popsizes <- c(vax,noncases_per_ratio[j]-vax)
      if(j<length(indices)) allocation_ratio <- response_adapt(fails,popsizes,days=day[j], adaptation=adaptation)
    }
    #weights <- get_weights_from_all_results(all_results)
    pval[i] <- calculate_pval(fails,popsizes)
  }
  
  return(quantile(pval,0.05))
}

fast_efficacy <- function(result_tab,vaccinees,trial_participants){
  ve_estimate <- c(0.5,1)
  break_count <- 0
  while(abs(ve_estimate[1]-ve_estimate[2])>0.005&&break_count<5&&ve_estimate[1]>0){
    
    pop_sizes2 <- c(sum(vaccinees),sum(trial_participants)-sum(vaccinees))
    if(nrow(result_tab)==0) return(list(0,pop_sizes2,c(0,0)))
    
    ve_estimate[2] <- ve_estimate[1]
    fails <- get_weighted_results_given_ve(result_tab,ve_point_est=ve_estimate[1],contact_network = -1,tested=F)
    
    if(fails[2]>0&&!any(pop_sizes2==0))
      ve_estimate[1] <- calculate_ve(fails,pop_sizes2)
    
    break_count <- break_count + 1
  }
  return(list(ve_estimate[1],pop_sizes2,fails))
      
}

get_efficacious_probabilities <- function(results_list,vaccinees,trial_participants,max_time=10000,contact_network=2,
                                          tested=F,randomisation_ratios=NULL,rbht_norm=0,people_per_ratio=NULL,adaptation='TST'){
  controls <- trial_participants - vaccinees
  if(is.null(randomisation_ratios)) randomisation_ratios <- rep(0.5,length(trial_participants))
  
  uninf_vacc <- vaccinees - sapply(results_list,function(x)sum(x$vaccinated))
  uninf_cont <- trial_participants - vaccinees - sapply(results_list,function(x)sum(x$inTrial&!x$vaccinated))
  
  uninf <- data.frame(vaccinated=c(rep(T,sum(uninf_vacc)),rep(F,sum(uninf_cont))),
                      allocRatio=c(rep(randomisation_ratios,uninf_vacc),rep(randomisation_ratios,uninf_cont)),
                      weight=1,infected=F)
  
  ve_estimate <- c(0.6,1)
  break_count <- 0
  not_nas <- lapply(1:length(results_list),function(x){
    results <- results_list[[x]]
    !is.na(results$RecruitmentDay)&results$RecruitmentDay<results$DayInfectious # subset(y,RecruitmentDay<DayInfectious)
  })
  if(contact_network==-1){
    results_tab_list <- list()
    for(x in 1:length(results_list)){
      results <- results_list[[x]]
      results$startDay <- x
      results$allocRatio <- randomisation_ratios[x]
      results_tab_list[[x]] <- results
    }
    result_tab <- do.call(rbind,results_tab_list)
    result_tab <- result_tab[unlist(not_nas),]
    result_tab$infected <- T
  }
  while(abs(ve_estimate[1]-ve_estimate[2])>0.005&&break_count<5&&ve_estimate[1]>0){
    if(contact_network>-1){
    results_tab_list <- list()
    for(x in 1:length(results_list)){
      results <- results_list[[x]]
      ##!! could include also RecruitmentDay
      w <- results#subset(results,DayInfected<max_time)
      #y <- subset(w,!is.na(RecruitmentDay))
      z <- c()
      if(sum(not_nas[[x]])>0) {
        weights <- get_infectee_weights(results=w,ve_point_est=ve_estimate[1],contact_network,tested)
        z <- w[not_nas[[x]],] # subset(y,RecruitmentDay<DayInfectious)
        z$startDay <- x
        z$allocRatio <- randomisation_ratios[x]
        z$infected <- T
        z$weight <- rowSums(weights[[1]])
      }
      results_tab_list[[x]] <- z
    }
    result_tab <- do.call(rbind,results_tab_list)
    }else{
      weights <- get_infectee_weights(results=result_tab,ve_point_est=ve_estimate[1],contact_network,tested)
      result_tab$weight <- rowSums(weights[[1]])
    }
      
    #result_tab$weight <- rowSums(get_infectee_weights(result_tab,ve_estimate[1],contact_network,tested)[[1]])
    ve_estimate[2] <- ve_estimate[1]
    
    if(nrow(result_tab)==0) return(list(0,c(sum(vaccinees),sum(trial_participants)-sum(vaccinees)),c(0,0)))
    
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
      maxJ <- nrow(people_per_ratio)
      if(nrow(all_results_original)>people_per_ratio[maxJ,1])
        maxJ <- nrow(people_per_ratio) + 1
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
        allocation_ratio <- response_adapt(weights[[1]],weights[[2]],days=people_per_ratio[1,2], adaptation)
        for(j in 2:maxJ){
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
          allocation_ratio <- response_adapt(weights[[1]],weights[[2]],days=max_t, adaptation)
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

get_efficacious_probabilities2 <- function(results_list,vaccinees,trial_participants,max_time=10000){
  ve_estimate <- c(0.6,1)
  weight_hh_rem <- matrix(0,ncol=2,nrow=length(results_list))
  break_count <- 0
  while(abs(ve_estimate[1]-ve_estimate[2])>0.005&&break_count<5){
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
    fails <- colSums(weight_hh_rem,na.rm=T)
    pop_sizes2 <- c(sum(vaccinees), sum(trial_participants) - sum(vaccinees))
    #print(c(2,des,pop_sizes2,fails[2]>0&&!any(pop_sizes2==0)))
    if(fails[2]>0&&!any(pop_sizes2==0))
      ve_estimate[1] <- calculate_ve(fails,pop_sizes2)
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
  end_day <- rec_day + 31
  # the total 'infectious force' exerted by each infectious person
  force_of_infection <- force_of_infection_after_0 <- c()
  for(x in 1:nrow(potential_infectors))
    force_of_infection[x] <- sum(pgamma_vector[end_day-infectious_days[x]:removal_days[x]])
  
  get_infectious_force_after_0 <- function(x){
    days <- infectious_days[x]:removal_days[x]
    days <- days[days>rec_day]
    sum(pgamma_vector[end_day-days])
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
  rec_day <<- max(netwk[[3]])
  results$DayRemoved[is.na(results$DayRemoved)] <- rec_day + eval_day
  potential_infectees <- netwk[[7]]
  if(pre_randomisation){
    potential_infectors <<- results # subset(results,DayRemoved>RecruitmentDay)
  }else{
    potential_infectors <<- results#[results$DayRemoved>results$RecruitmentDay,]
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
        infectee_weights <- get_infectee_weights(results,ve_est_temp,correct_for_ve=F)
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

## methods 7 (pre_randomisation=F) and 8 (pre_randomisation=T)
iterate_ph_model <- function(netwk_list,cluster_flag=0,pre_randomisation=T){
  ves <- c(0,1)
  break_count <- 0
  while(abs(ves[1]-ves[2])>0.005&&break_count<5){
    #trial_summary <- lapply(netwk_list,summarise_trial,ve_est_temp=ves[1],pre_randomisation=pre_randomisation)
    trial_summary <- list()
    for(cluster in 1:length(netwk_list)) {
      ind <- length(trial_summary) + 1
      trial_summary[[ind]] <- summarise_trial(netwk_list[[cluster]],ve_est_temp=ves[1],pre_randomisation=pre_randomisation)
      if(ind==length(trial_summary))
        trial_summary[[ind]] <- cbind(trial_summary[[ind]],cluster)
    }
    
    tte <- do.call(bind_rows,trial_summary)
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
    break_count <- break_count + 1
  }
  return(c(pval,ves[1]))
}
