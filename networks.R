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

## get parameters from data ###########################################

get_lnorm_params <- function(mean,sd){
  mu=-log(((sd/mean)^2+1)/(mean^2))/2
  sig2=2*(log(mean)-mu)
  c(mu,sqrt(sig2))
}

get_gamma_params <- function(mean,sd){
  shape=1/sd^2*mean^2
  rate = shape/mean
  c(shape,rate)
}

g1_params <- c(get_gamma_params(9.7,5.3),51)
g2_params <- c(get_gamma_params(11,4.1),47)
time_to_randomisation <- c(rgamma(g1_params[3],shape=g1_params[1],rate=g1_params[2]),rgamma(g2_params[3],shape=g2_params[1],rate=g2_params[2]))
hist(time_to_randomisation)
hist(rtruncnorm(1000,a=0,mean=mean(time_to_randomisation),sd=sd(time_to_randomisation)))
get_gamma_params(mean(time_to_randomisation),sd(time_to_randomisation))
for(i in 1:100){
  time_to_randomisation <- c(rgamma(g1_params[3]*1000,shape=g1_params[1],rate=g1_params[2]),rgamma(g2_params[3]*1000,shape=g2_params[1],rate=g2_params[2]))
  g <- get_gamma_params(mean(time_to_randomisation),sd(time_to_randomisation))
  #print(c(sum(rtruncnorm(1000,a=0,mean=mean(time_to_randomisation),sd=sd(time_to_randomisation))>20)-
  #          sum(rgamma(1000,shape=g[1],rate=g[2])>20)))
}

g1_params <- c(get_gamma_params(3.9,2.9),51)
g2_params <- c(get_gamma_params(3.8,2.6),47)
for(i in 1:100){
time_to_admission <- c(rgamma(g1_params[3]*1000,shape=g1_params[1],rate=g1_params[2]),rgamma(g2_params[3]*1000,shape=g2_params[1],rate=g2_params[2]))
g <- get_gamma_params(mean(time_to_admission),sd(time_to_admission))
#print(c(sum(rtruncnorm(1000,a=0,mean=mean(time_to_admission),sd=sd(time_to_admission))>10)-
#sum(rgamma(1000,shape=g[1],rate=g[2])>10)))
}
#hist(time_to_admission)
#hist(rtruncnorm(1000,a=0,mean=mean(time_to_admission),sd=sd(time_to_admission)))

g1_params <- c(get_lnorm_params(3.9,2.9),51)
g2_params <- c(get_lnorm_params(3.8,2.6),47)
time_to_admission <- c(rlnorm(g1_params[3],g1_params[1],g1_params[2]),rlnorm(g2_params[3],g2_params[1],g2_params[2]))
hist(time_to_admission)

## things to do ##########################################

##!! Add delay to notification/hospitalisation

##!! relationship between cluster size and time to recruit

##!! fit data. aiming for minimum 71/111833=0.00063 cases per index case. can we get data to see what happens before randomisation?

##!! what should the distribution of number infected per cluster be?

##!! add in infect from source?

##!! what is the admission time for contacts? same as for index case, or less? and should they be correlated?

##!! maybe there should be a difference between contact edges and transmission edges?
## or weight family contacts more heavily? or a behaviour change following enrollment?

##!! perhaps two types of edge, and a means to identify high-risk contacts (15%)
## Contacts are defined as individuals who, within the previous 21 days, lived in the 
## same household as the symptomatic patient, were visited by the symptomatic patient,
## or were in close physical contact with the patient's body or body fluids, linen, or clothes.
## Contacts of contacts are defined as the neighbours or extended family members to nearest 
## geographic boundary in which the local contacts of the index case reside, plus household 
## members of any high risk contacts who do not live in the same locality as the case (further 
## details are in the accompanying protocol). 

## Aiming for average number of connections 17.5
## with 70 "contacts of contacts" (i.e. neighbours, extended family and contacts of high-risk contacts)


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
plot(induced_subgraph(new_g,young_people))

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
plot(g2)

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
barplot(degreedistribution,ylab='Number of people', xlab='Number of connections',names.arg=0:(length(degreedistribution)-1),main='')
average_contacts <- sum(degreedistribution*c(1:length(degreedistribution)-1)/length(E(neighbourhood_g)))
rm(neighbour_adjacency_matrix)
pdf('three_hh.pdf'); par(mar=c(0,0,0,0))
plot(induced_subgraph(new_g,1:20),vertex.color=rep(c('navyblue','hotpink','grey'),times=c(5,11,4)), vertex.label=NA)
dev.off()
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

## set up #######################################################

# Per-time-step hazard of infection for a susceptible nodes from an infectious
# neighbour
beta_base <- 0.0065
high_risk_scalar <- 2.17
# fraction of beta_base applied to neighbours ("contacts of contacts")
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

## start ############################################################

#profvis({
for(iter in 1:nIter){
  ## select random person to start
  first_infected <- sample(g_name,1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
  hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
  inf_time <- min(inf_period,hosp_time)
  netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,inf_time,start_day=iter,from_source=per_time_step,cluster_flag=cluster_flag)
  
  results_list[[iter]] <- netwk[[1]]
  cluster_size[iter] <- netwk[[2]]
  recruit_times[iter] <- netwk[[3]]
  hosp_times[iter] <- inf_time
  
}
#})

## results #######################################################

number_infectious <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster)
}
)

number_infectious_after_randomisation <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfectious>results$RecruitmentDay)
}
)
number_infected_after_randomisation <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfected>results$RecruitmentDay)
}
)

number_infectious_before_day10 <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfectious<results$RecruitmentDay+10&results$DayInfectious>results$RecruitmentDay-2)
}
)
contacts_infectious_before_day10 <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$contact&results$DayInfectious<results$RecruitmentDay+10&results$DayInfectious>results$RecruitmentDay-2)
}
)
c(sum(contacts_infectious_before_day10),sum(number_infectious_before_day10))
number_infectious_after_randomisation_before_day10 <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10)
}
)
number_infected_after_randomisation_before_day10 <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$inCluster&results$DayInfected>results$RecruitmentDay&results$DayInfected<results$RecruitmentDay+10)
}
)

secondary_infections <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$DayInfected>hosp_times[iter])
}
)
infectious_on_recruitment <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  sum(results$DayInfectious[-1]<recruit_times[iter])
}
)
infectious_by_vaccine <- sapply(1:length(cluster_size),function(iter){
  results <- results_list[[iter]]
  c(sum(results$vaccinated&results$DayInfectious>results$RecruitmentDay+9),sum(results$inTrial&results$DayInfectious>results$RecruitmentDay+9))
}
)

c(sum(contacts_infectious_before_day10),sum(number_infectious_before_day10))
sum(number_infectious_before_day10-number_infectious_after_randomisation_before_day10)/sum(cluster_size)
sum(number_infectious_after_randomisation_before_day10)/sum(cluster_size)
sum(number_infectious_before_day10)/sum(cluster_size)

sapply(0:6,function(x)sum(number_infectious-number_infectious_before_day10==x))

quantile(cluster_size,c(0.25,0.5,0.75))

variables <- list(cs=cluster_size,ht=hosp_times,rt=recruit_times)
outcomes <- list(number_infectious_before_day10,number_infectious_after_randomisation_before_day10,number_infected_after_randomisation_before_day10,secondary_infections,infectious_on_recruitment)
zeros <- sapply(outcomes,function(x)sum(x==0)/length(x))*100
means <- sapply(outcomes,function(x)mean(x))
ranges <- apply(sapply(outcomes,function(x)quantile(x[x>0],c(0.05,0.95))),2,function(x)paste0(x,collapse='--'))
outcome_labels <- c('Number infectious, day < 10','Number infectious, 0 < day < 10',
                    'Number infected, 0 < day < 10','Number of secondary infections','Number infectious on recruitment')
df <- data.frame(cbind(zeros,means,ranges),row.names=outcome_labels,stringsAsFactors = F)
for(i in 1:2) df[,i] <- as.numeric(df[,i])
xtable(df,digits=c(1,1,2,1))

##!! doesn't include the 12 excluded
delay_before_day_ten <- c(rep(0,30+3),4,3,3,3,rep(2,5),rep(1,4))
denominator <- 3096+1461
sum(delay_before_day_ten)/denominator
hist(delay_before_day_ten)
den <- c(rnbinom(1000,3.5,0.1),rep(0,2000))
points(sort(unique(den)),sapply(sort(unique(den)),function(x)sum(den==x))/50)

probability_by_lag <- readRDS(paste0('probability_by_lag_8080.Rds'))
probability_after_day_0 <- readRDS(paste0('probability_after_day_0_8080.Rds'))
probability_by_lag_given_removal <- readRDS(paste0('probability_by_lag_given_removal_8080.Rds'))
probability_after_day_0_given_removal <- readRDS(paste0('probability_after_day_0_given_removal_8080.Rds'))

ref_recruit_day <- 30
true_vs_weight <- c(0,0,0,0,0)
how_surprising <- c()
for(i in 1:length(results_list)){
  results <- results_list[[i]]
  if(nrow(results)>1){
    weights <- get_weighted_results(results,how_surprising)
    true_vs_weight <- rbind(true_vs_weight,weights[[1]])
    how_surprising <- weights[[2]]
  }
}
hist(how_surprising)
quantile(how_surprising,c(0.95,0.975,0.99))
sapply(2:5,function(x)sum(abs(true_vs_weight[,1]-true_vs_weight[,x])))
colSums(true_vs_weight)
nrow(true_vs_weight)
pdf('count_v_estimate.pdf'); par(mfrow=c(1,1),mar=c(5,5,2,2),cex=1.5,pch=20)
plot(true_vs_weight[,1],true_vs_weight[,3],col='navyblue',ylab='Estimated count',xlab='True count',frame=F)
points(true_vs_weight[,1],true_vs_weight[,2],cex=0.8,col='hotpink')
lines(c(0,max(true_vs_weight)),c(0,max(true_vs_weight)),col='grey',lty=2,lwd=2)
legend(x=0,y=22,legend=c('Binary weights','Continuous weights'),col=c('navyblue','hotpink'),pch=20,cex=0.75,bty='n')
dev.off()

plot(sapply(results_list,nrow))

## can we infer a trend? ##################################################

direct_VE <- 0.9
nIter <- 100
adaptation <- 'TST'
pval_binary_mle2 <- ve_est2 <- c()
eval_day <- 31
latest_infector_time <- eval_day - 0
func <- get_efficacious_probabilities2
#for(per_time_step in c(0,1e-10,1e-9,1e-8,1e-7)){
for(rep in 1:1000){
  per_time_step <- 0
  #profvis({
  allocation_ratio <- 0.5
  results_list <- list()
  vaccinees <- trial_participants <- c()
  vaccinees2 <- trial_participants2 <- c()
  for(iter in 1:nIter){
    ## select random person to start
    first_infected <- sample(g_name,1)
    inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
    # Add them to e_nodes and remove from s_nodes and v_nodes
    #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
    hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
    inf_time <- min(inf_period,hosp_time)
    netwk <- simulate_contact_network(neighbour_scalar,high_risk_scalar,first_infected,inf_time,end_time=eval_day,start_day=iter,from_source=per_time_step,
                                      cluster_flag=1,allocation_ratio=allocation_ratio,direct_VE=direct_VE)
    
    results_list[[iter]] <- netwk[[1]]
    cluster_size[iter] <- netwk[[2]]
    recruit_times[iter] <- netwk[[3]]
    hosp_times[iter] <- inf_time
    rec_day <- recruit_times[iter]
    infectious_index <- results_list[[iter]]$DayInfectious<latest_infector_time+rec_day&(results_list[[iter]]$DayRemoved>rec_day|is.na(results_list[[iter]]$DayRemoved))
    infectious_names <- results_list[[iter]]$InfectedNode[infectious_index]
    infectious_ends <- pmin(results_list[[iter]]$DayRemoved[infectious_index],latest_infector_time+rec_day)
    infectious_ends[is.na(infectious_ends)] <- latest_infector_time+rec_day
    infectious_starts <- pmax(results_list[[iter]]$DayInfectious[infectious_index],rec_day)
    vaccinees[iter] <- trial_participants[iter] <- 0
    if(length(infectious_names)>0&&identical(func,get_efficacious_probabilities2)){
      popweights <- rowSums(sapply(1:length(infectious_names),function(i){
        x <- infectious_names[i]
        # prepare contacts
        contacts <- contact_list[[x]]
        c_of_c <- contact_of_contact_list[[x]]
        hr <- c(high_risk_list[[x]],household_list[[x]])
        # prepare trial participants
        vax <- netwk[[6]]
        cont <- netwk[[7]]
        # work out total risk presented by infector
        infector_weight <- sum(pgamma(eval_day-infectious_starts[i]:infectious_ends[i],shape=incperiod_shape,rate=incperiod_rate))
        # remove anyone infectious earlier
        earlier_nodes <- results_list[[iter]]$InfectedNode[results_list[[iter]]$DayInfectious<infectious_starts[i]]
        contacts <- contacts[!contacts%in%earlier_nodes]
        c_of_c <- c_of_c[!c_of_c%in%earlier_nodes]
        hr <- hr[!hr%in%earlier_nodes]
        # sum of person days times scalars
        c(sum(vax%in%contacts) + neighbour_scalar*sum(vax%in%c_of_c) + (high_risk_scalar-1)*sum(vax%in%hr),
          sum(cont%in%contacts) + neighbour_scalar*sum(cont%in%c_of_c) + (high_risk_scalar-1)*sum(cont%in%hr))*infector_weight
        }))
      if(length(netwk[[6]])>0)
        vaccinees[iter] <- popweights[1]
      trial_participants[iter] <- popweights[2]
    }
    #if(identical(func,get_efficacious_probabilities)){
      vaccinees2[iter] <- netwk[[4]]
      trial_participants2[iter] <- netwk[[5]]
    #}
    if(adaptation!=''&&iter %% eval_day == 0){
      allocation_ratio <- response_adapt(results_list,vaccinees,trial_participants,adaptation,func=func)
    }
  }
  #})
  {
  #days_infectious <- unlist(sapply(1:length(results_list),function(x) x+results_list[[x]]$DayInfectious))
  #days <- 1:nIter
  #counts <- sapply(days,function(x)sum(days_infectious==x))-1
  #dataset <- data.frame(t=days,clusters=31,count=counts)
  #dataset$clusters[1:31] <- 1:31
  #mod <- glm(count~t,data=dataset,offset=log(clusters),family=poisson(link=log))
  }
  eval_list <- func(results_list,vaccinees,trial_participants)
  pval_binary_mle2[rep]  <- calculate_pval(eval_list[[3]],eval_list[[2]])
  ve_est2[rep] <- eval_list[[1]]
  #print(c(pval_binary_mle2,ve_est2,allocation_ratio))
}
sum(pval_binary_mle2<0.05,na.rm=T)/sum(!is.na(pval_binary_mle2))
mean(ve_est2)
sd(ve_est2)
hist(rpois(1000,mean(counts-1))+1)
hist(counts)


## infections #############################################
## if we know removal day
removal_days <- c()
for(i in 1:length(results_list)){
  results <- results_list[[i]]
  removed <- subset(results,!is.na(DayRemoved))
  removal_days <- c(removal_days,removed$DayRemoved-removed$DayInfectious)
}


rm(results_list)

## mi #################################################################

evppi <- sapply(variables,function(x)
  sapply(outcomes,
         function(y)mutinformation(discretize(x),y)
  )
)

{pdf('evppi.pdf',height=10,width=4+length(outcomes)); 
#  {x11(height=10,width=4+length(outcomes)); 
    par(mar=c(12,24,3.5,10))
  labs <- rev(c('Cluster size','Index removal time','Recruitment time'))
  get.pal=colorRampPalette(brewer.pal(9,"Reds"))
  redCol=rev(get.pal(12))
  bkT <- seq(max(evppi)+1e-10, 0,length=13)
  cex.lab <- 1.5
  maxval <- round(bkT[1],digits=1)
  col.labels<- c(0,maxval/2,maxval)
  cellcolors <- vector()
  for(ii in 1:length(unlist(evppi)))
    cellcolors[ii] <- redCol[tail(which(unlist(evppi[ii])<bkT),n=1)]
  color2D.matplot(evppi,cellcolors=cellcolors,main="",xlab="",ylab="",cex.lab=2,axes=F,border='white')
  fullaxis(side=2,las=2,at=0:(length(outcomes)-1)+1/2,labels=rev(outcome_labels),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.2)
  fullaxis(side=1,las=2,at=(length(labs)-1):0+0.5,labels=labs,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.2)
  mtext(3,text='Mutual information',line=1,cex=1.3)
  color.legend(length(labs)+0.5,0,length(labs)+0.8,length(outcomes),col.labels,rev(redCol),gradient="y",cex=1.2,align="rb")
  for(i in seq(0,length(outcomes),by=1)) abline(h=i)
  for(i in seq(0,length(labs),by=1)) abline(v=i)
  #dev.off()
  }


## abc ##################################
data <- c()
# are there any cases before day 0
data[1] <- 16/11833
# % before day 0 = 16 per 11833 people = 0.0014 per person
data[2] <- 16/11833
# 0 in cluster days days 0 to 10
data[3:7] <- sapply(0:4,function(x)sum(delay_before_day_ten==x)/length(delay_before_day_ten))
# pc infectees are contacts
infected_clusters <- sum(c(11,6,10,21,6))
of_which_contacts <- sum(c(10,6,10,20,6))
data[8] <- of_which_contacts/infected_clusters
# pc infectees are high risk
of_which_hr <- sum(c(10,5,10,20,5))
data[9] <- of_which_hr/infected_clusters
cluster_sum <- sapply(0:4,function(x)sum(delay_before_day_ten==x))

zeros <- c()
sumll <- c()

distnce <- distnce_sample <- c()
distnce[1] <- 50
parameter_samples <- list(beta_hyper=list(vals=c(),star=c()),
                          nb_hyper=list(vals=c(),star=c()),
                          hr_hyper=list(vals=c(),star=c())
)
param_names=c('beta_base','neighbour_scalar','high_risk_scalar')
param_quantiles <- c(qlnorm,qbeta,function(x,...){1/qbeta(x,...)})
param_hypers <- list(c(log(0.01),0.7),c(2,2),c(5,5))
successes <- rep(0,length(parameter_samples))
step_size <- rep(0.05,length(parameter_samples))
for(i in 1:length(parameter_samples)){
  parameter_samples[[i]]$vals[1] <- parameter_samples[[i]]$star <- runif(1)
}
direct_VE <- 0
#profvis({
for(iter in 1:10000){
  
  
  for(param in 1:length(parameter_samples)){
    parameter_sample <- parameter_samples[[param]]
    star <- parameter_sample$star
    val <- parameter_sample$val
    assign(param_names[param], 
           param_quantiles[param][[1]](star,param_hypers[[param]][1],param_hypers[[param]][2]))
  
    
    est <- matrix(0,ncol=10,nrow=length(delay_before_day_ten))
    for(i in 1:length(delay_before_day_ten)){
      ## select random person to start
      first_infected <- sample(g_name,1)
      inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
      # Add them to e_nodes and remove from s_nodes and v_nodes
      #hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
      hosp_time <- rtruncnorm(length(first_infected),a=0,mean=hosp_mean_index,sd=hosp_sd_index)
      inf_time <- min(inf_period,hosp_time)
      netwk <- simulate_contact_network(beta_base,neighbour_scalar,high_risk_scalar,first_infected,inf_time,end_time=10,cluster_flag=cluster_flag)
      
      results <- netwk[[1]]
      cluster_sz <- netwk[[2]]
      #recruit_times[iter] <- netwk[[3]]
      #hosp_times[iter] <- inf_time
      
      # number before day 0
      est[i,1] <- sum(results$inCluster&results$DayInfectious<results$RecruitmentDay+1)
      
      # 0 in cluster days days 0 to 10
      est[i,2:6] <- sapply(0:4,function(x)sum(results$inCluster&results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10)==x)
      # cases in cluster days days 0 to 10
      est[i,7] <- sum(results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10&results$inCluster)
      # pc infectees are contacts
      est[i,8] <- sum(results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10&results$contact)
      # pc infectees are high risk
      est[i,9] <- sum(results$inCluster&results$DayInfectious>results$RecruitmentDay&results$DayInfectious<results$RecruitmentDay+10&results$highrisk==T)
      # cluster size
      est[i,10] <- cluster_sz
    }
    distnce_tmp <- c()
    # distance from number before day 0
    distnce_tmp[1] <- abs(data[1]-sum(est[,1])/sum(est[,10]))*100
    # distances for number of clusters with 0, 1, 2 etc cases
    distnce_tmp[2:6] <- abs(colSums(est[,2:6]) - cluster_sum)
    distnce_tmp[7] <- abs(sum(est[,8])/sum(est[,7]) - of_which_contacts/infected_clusters)*10
    distnce_tmp[8] <- abs(sum(est[,9])/sum(est[,7]) - of_which_hr/infected_clusters)*10
    distncestar <- sum(distnce_tmp[is.finite(distnce_tmp)])
    #print(c(distnce_tmp,distncestar))
    
    distnce_sample[iter] <- distncestar
    if(distncestar < 30){ 
      val[iter+1] <- star
      distnce[iter+1] <- distncestar
      successes[param] <- successes[param] + 1
    }else{
      val[iter+1] <- val[iter]
      distnce[iter+1] <- distnce[iter]
    }
    star <- val[iter+1] + rnorm(1,0,step_size[param])
    while(star<0) star <- star+ceiling(abs(star))
    while(star>1) star <- star-floor(star)
    
    if(iter%%20==0){
      acceptance_ratio <- successes[param]/20
      step_size[param] <- step_size[param]*((acceptance_ratio+0.05)/0.4)^0.25
      successes[param] <- 0
    }
    parameter_sample$star <- star
    parameter_sample$val <- val
    parameter_samples[[param]] <- parameter_sample
  }
  
  zeros[iter] <- est[2]==0
  sumll[iter] <- sum(distnce_tmp[is.finite(distnce_tmp)])
}
  #})
par(mar=c(4,4,2,2))
hist(qlnorm(runif(1000),log(0.01),0.7))
hist(qlnorm(parameter_samples[[1]]$val,log(0.01),0.7))
mean(qlnorm(parameter_samples[[1]]$val,log(0.01),0.7))
plot(parameter_samples[[1]]$val)
plot(parameter_samples[[2]]$val)
hist(qbeta(runif(1000),2,2))
hist(qbeta(parameter_samples[[2]]$val,2,2))
mean(qbeta(parameter_samples[[2]]$val,2,2))
plot(parameter_samples[[2]]$val,parameter_samples[[3]]$val)
plot(parameter_samples[[3]]$val)
hist(1/qbeta(runif(1000),5,5))
hist(1/qbeta(parameter_samples[[3]]$val,5,5))
mean(1/qbeta(parameter_samples[[3]]$val,5,5))
sum(zeros)
plot(distnce_sample)
