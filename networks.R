
setwd('~/overflow_dropbox/ADAGIO/')
source('Code/functions_network.R')
library(igraph)

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
get_gamma_params(mean(time_to_randomisation),sd(time_to_randomisation))

g1_params <- c(get_gamma_params(3.9,2.9),51)
g2_params <- c(get_gamma_params(3.8,2.6),47)
time_to_admission <- c(rgamma(g1_params[3],shape=g1_params[1],rate=g1_params[2]),rgamma(g2_params[3],shape=g2_params[1],rate=g2_params[2]))
hist(time_to_admission)
get_gamma_params(mean(time_to_admission),sd(time_to_admission))


## things to do ##########################################

##!! Add delay for development of immunisation

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
household_sizes <- sample(6:10,number_of_households,replace=T)

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

# add high-risk labels, to be used for ring vaccination, could be used to increase disease spread
high_risk_label <- rep(0,length(V(new_g)))
high_risk_label[sample(V(new_g),round(length(V(new_g))*0.15))] <- 1
new_g <- set_vertex_attr(new_g,'high_risk',value=high_risk_label)

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
# aiming for average contacts approx 60
##!! there are almost certainly duplicate edges here, so some people might get two tries to infect someone. Is that what we want?

## simulate #######################################################


spread <- function(g, s_nodes, v_nodes, e_nodes, i_nodes,c_nodes, beta, direct_VE,
                   incperiod_shape, incperiod_rate){
  # Spread will create new infected nodes from two sources: infectious nodes within the the study
  # population, and external pressure from the source population
  # Inputs:
  # g is the graph, used to find neighbours of infected nodes
  # s_nodes, e_nodes and i_nodes are susceptible, exposed and infected nodes
  # beta is the hazard of infection for one contact
  # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
  # length, currently drawn from a gamma distribution
  # connected_nodes is a list of nodes that are connected to the source population
  
  # Process: go through list of i_nodes, and choose a random number of its susceptible neighbours to
  # be infected, according to beta and choose a random number of its susceptible vaccinated neighbours to
  # be infected, according to beta and direct_VE
  infectees_susc <- c()
  infectees_vacc <- c()
  infectees_cont <- c()
  if (ncol(i_nodes)>0) {
    # Make a beta vector
    beta_v <- beta*(1-direct_VE)
    # Get a list of all neighbours of all infected nodes
    potential_contacts <- unlist(ego(g,order=1,nodes=i_nodes[1,]))
    infectees_susc <- infect_contacts(potential_contacts,node_class=s_nodes,beta_value=beta)
    ##!! this returns the infectious  node(s)
    neighbour_hh <- unlist(ego(g2,order=1,nodes=V(g)$hh[i_nodes[1,]]))
    neighbours <- unlist(ego(neighbourhood_g,order=1,nodes=i_nodes[1,]))
    infectees_susc <- unique(c(infectees_susc,infect_contacts(potential_contacts,node_class=s_nodes,beta_value=beta*neighbour_scalar)))
  } 
  
  newinfected_susc <- infectees_susc
  newinfected_vacc <- infectees_vacc
  newinfected_cont <- infectees_cont
  newinfected <- c(newinfected_susc, newinfected_vacc,newinfected_cont)
  newinfected <- unique(newinfected)
  
  if (length(newinfected)>0) {
    # Give each newly exposed node an incubation/latent period
    inc_periods <- rgamma(length(newinfected),shape=incperiod_shape,rate=incperiod_rate)
    # Add them to e_nodes and remove from s_nodes and v_nodes
    e_nodes <- cbind(e_nodes,rbind(newinfected,rep(0,length(newinfected)),inc_periods))
    s_nodes<-setdiff(s_nodes,newinfected_susc)
    v_nodes <- setdiff(v_nodes,newinfected_vacc)
    c_nodes <- setdiff(c_nodes,newinfected_cont)
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
  incubation_days <- as.vector(e_nodes[2,indices_to_remove])+1
  
  # Add one day to the length of each infected individual's time infected
  e_nodes[2,] <- e_nodes[2,]+1
  
  # Remove any progressing from e_nodes and add to i_nodes
  e_nodes <- e_nodes[,!(e_nodes[1,] %in% newinfectious),drop=FALSE]
  inf_periods <- rgamma(length(newinfectious),infperiod_shape,infperiod_rate)
  hosp_time <- rgamma(length(newinfected),shape=hosp_shape,rate=hosp_rate)
  i_nodes <- cbind(i_nodes,rbind(newinfectious,rep(0,length(newinfectious)),pmin(inf_periods,hosp_time),incubation_days))
  
  list(e_nodes, i_nodes, r_nodes, sort(newinfectious))
}

# Per-time-step hazard of infection for a susceptible nodes from an infectious
# neighbour
beta <- 0.03
# fraction of beta applied to neighbours ("contacts of contacts")
neighbour_scalar <- 0.5
# Gamma-distribution parameters of incubation and infectious period and wait times
incperiod_shape <- 3.11
incperiod_rate <- 0.32
infperiod_shape <- 1.13
infperiod_rate <- 0.226
hosp_shape_index <- 2
hosp_rate_index <- 0.5
hosp_shape <- 2
hosp_rate <- 2
recruit_shape <- 5.4
recruit_rate <- 0.47

disease_dynamics <- list(beta=beta,
                         neighbour_scalar=neighbour_scalar,
                         incperiod_shape=incperiod_shape,
                         incperiod_rate=incperiod_rate,
                         infperiod_shape=infperiod_shape,
                         infperiod_rate=infperiod_rate,
                         hosp_shape_index=hosp_shape_index,
                         hosp_rate_index=hosp_rate_index,
                         hosp_shape=hosp_shape,
                         hosp_rate=hosp_rate,
                         recruit_shape=recruit_shape,
                         recruit_shape=recruit_shape)
for(i in 1:length(disease_dynamics)) assign(names(disease_dynamics)[i],disease_dynamics[[names(disease_dynamics)[i]]])

g <<- new_g

g_name <- V(g)$name
vertices <- V(g)
number_infectious <- cluster_size <- c()
number_infected_after_randomisation <- number_infectious_after_randomisation <- c()
for(iter in 1:1000){
  trajectories <- list()
  trajectories$S <- c()
  trajectories$E <- c()
  trajectories$I <- c()
  trajectories$R <- c()
  
  e_nodes <- matrix(nrow=3,ncol=0)
  i_nodes <- matrix(nrow=4,ncol=0)
  v_nodes <- c()
  c_nodes <- c()
  s_nodes <- as.vector(V(g))
  r_nodes <- c()
  
  ## select random person to start
  first_infected <- sample(s_nodes,1)
  inf_period <- rgamma(length(first_infected),shape=infperiod_shape,rate=infperiod_rate)
  # Add them to e_nodes and remove from s_nodes and v_nodes
  hosp_time <- rgamma(length(first_infected),shape=hosp_shape_index,rate=hosp_rate_index)
  inf_time <- min(inf_period,hosp_time)
  inc_time <- rgamma(length(first_infected),shape=incperiod_shape,rate=incperiod_rate)
  i_nodes <- cbind(i_nodes,rbind(first_infected,rep(0,length(first_infected)),inf_time,inc_time))
  s_nodes<-setdiff(s_nodes,first_infected)
  
  recruitment_time <- round(rgamma(1,shape=recruit_shape,rate=recruit_rate))
  results <- data.frame("InfectedNode"=numeric(),
                        "DayInfectious"=numeric(),
                        "RecruitmentDay"=numeric())
  results <- rbind(results,data.frame("InfectedNode"=first_infected,
                                      "DayInfectious"=0,
                                      "RecruitmentDay"=recruitment_time,
                                      "DayInfected"=-inc_time))
  numinfectious <- 1
  
  sim_time <- recruitment_time + 10
  for(t in 1:sim_time){
    
    newinfectious <- c()
    if ((ncol(i_nodes)>0)||(ncol(e_nodes)>0)) 
      list[e_nodes,i_nodes,r_nodes,newinfectious] <- recover(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate)
    #sub_g <- if(ncol(i_nodes)>0) induced_subgraph(g,c(i_nodes[1,],s_nodes,v_nodes)) else NULL
    list[s_nodes,v_nodes,e_nodes,c_nodes] <- spread(g,s_nodes,v_nodes,e_nodes,i_nodes,c_nodes,
                                                    beta,direct_VE,incperiod_shape,incperiod_rate)
    numnewinfectious <- length(newinfectious)
    if (numnewinfectious>0) {
      #new_indices <- g_name %in% newinfectious
      #v_subset <- V(g)[new_indices]
      # Update results
      results <- rbind(results,data.frame("InfectedNode"=newinfectious,
                                          "DayInfectious"=t,
                                          "RecruitmentDay"=recruitment_time,
                                          "DayInfected"=i_nodes[4,match(newinfectious,i_nodes[1,])]))
      
      numinfectious <- numinfectious+numnewinfectious
    }
    
    trajectories$S <- c(trajectories$S,length(s_nodes) + length(v_nodes) + length(c_nodes))
    trajectories$E <- c(trajectories$E,ifelse(length(trajectories$E)==0,0,-diff(tail(trajectories$S,2))))
    trajectories$I <- c(trajectories$I,numnewinfectious)
    trajectories$R <- c(trajectories$R,ifelse(length(trajectories$R)==0,0,length(r_nodes)-sum(trajectories$R)))
    
    if (numinfectious>0) {
      results<-results[1:numinfectious,]
    } else {
      results <- results[1,]
    }
  }
  number_infectious[iter] <- nrow(results)-1
  number_infectious_after_randomisation[iter] <- sum(results$DayInfectious>results$RecruitmentDay)
  number_infected_after_randomisation[iter] <- sum(results$DayInfected>results$RecruitmentDay)
  cluster_size[iter] <- ego_size(g,order=1,nodes=first_infected) + ego_size(neighbourhood_g,order=1,nodes=first_infected) - 2
}

mean(number_infectious/cluster_size)*11833
sum(number_infectious==0)
mean(number_infectious_after_randomisation/cluster_size)*11833
sum(number_infectious_after_randomisation==0)
mean(number_infected_after_randomisation/cluster_size)*11833
sum(number_infected_after_randomisation==0)
quantile(cluster_size,c(0.25,0.5,0.75))
