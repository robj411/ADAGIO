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
    if(length(young_people[!young_people%in%c(first_person,ego(new_g,order=1,nodes=first_person)[[1]])])>0){
      second_person <- sample(young_people[!young_people%in%c(first_person,ego(new_g,order=1,nodes=first_person)[[1]])],1,replace=F)
      new_g <- add_edges(new_g,edges=c(first_person,second_person))
    }
  }
}

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
contact_list <<- lapply(V(new_g),function(x) {cs <- as.vector(unlist(ego(new_g,order=1,nodes=x))); cs[cs!=x]})
mean(sapply(contact_list,length))

## get neighbourhood network
# assume 5 hh per neighbourhood
n_hood_size <- 6
neighbourhood_sizes <- rep(n_hood_size,length=floor(number_of_households/n_hood_size)-1)
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
average_contacts <- sum(degreedistribution*c(1:length(degreedistribution)-1)/length(E(neighbourhood_g)))
rm(neighbour_adjacency_matrix)
# aiming for average contacts approx 60
##!! there are almost certainly duplicate edges here, so some people might get two tries to infect someone. Is that what we want?

# get list of neighbours
contact_of_contact_list <<- lapply(V(neighbourhood_g),function(x) {cofc <- as.vector(unlist(ego(neighbourhood_g,order=1,nodes=x))); cofc[cofc!=x]})

household_list <<- lapply(V(new_g),function(x){hh_members <- which(hh_labels==hh_labels[x]); as.vector(hh_members[hh_members!=x])})

# add high-risk labels, to be used for ring vaccination, could be used to increase disease spread
# assume high risk rate is constant across contacts and contacts of contacts
high_risk_rate <- sum(c(330,171,58,246,574,231))/sum(2151,1435,1104,1678,3796,2572)
high_risk_list <<- lapply(V(new_g),function(x){
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
