occupancy_data <- read_xls('ct08192011censushhtypehhsizeandageofusualresidentshouseholdsenglandandwales.xls',sheet=2)
colnames(occupancy_data) <- occupancy_data[9,]
occupancy_data <- occupancy_data[-c(1:9),]
#occupancy_data$`Household Size`[occupancy_data$`Household Size`=="11 or more people in household"] <- '11 people or more'
occupancy_data$number <- as.numeric(sapply(occupancy_data$`Household Size`,function(x)strsplit(x,' p')[[1]][1]))
occupancy_data <- subset(occupancy_data,!is.na(number))
#occupancy_data[286,5] <- as.numeric(occupancy_data[286,6])
#occupancy_data[286,7] <- as.numeric(occupancy_data[286,8])
households <- sapply(1:10,function(y)sum(as.numeric(subset(occupancy_data,number==y)$Total)))
people <- sapply(1:10,function(y)sum(as.numeric(subset(occupancy_data,number==y)$Total))*y)
barplot(households)
barplot(people)
occupancy_data$type <- 1:nrow(occupancy_data)

communal <- read_xls('communal_living.xls',sheet=1)
economic <- read_xls('communal_living_economic.xls',sheet=1)


## build network ###############################################################

number_of_households <- 500
household_types <- sample(occupancy_data$type,number_of_households,replace=T,prob=occupancy_data$Total)
colnames(occupancy_data)[1:4] <- c('description','children','adults','elderly')
household_sizes <- occupancy_data$number[household_types]

# assign individuals to households
label_start <- 0
hh <- list()
for(i in 1:number_of_households) {
  hh[[i]] <- make_full_graph(household_sizes[i]) %>%
    set_vertex_attr("name", value = label_start+1:household_sizes[i])
  label_start <- label_start + household_sizes[i]
}

# extract household data frames
attrs <- do.call(rbind,lapply(hh,function(x)igraph::as_data_frame(x,'vertices')))
# combine all
el <- do.call(rbind,lapply(hh,function(x)igraph::as_data_frame(x)))
# convert to network
new_g <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)
# save layout for plotting
pts <- 80
save_layout <- layout_nicely(induced.subgraph(new_g,1:pts))
# add household labels
hh_labels <- rep(1:number_of_households,household_sizes)
new_g <- set_vertex_attr(new_g,'hh',value=hh_labels)

household_list <<- lapply(V(new_g),function(x) {cs <- as.vector(unlist(ego(new_g,order=1,nodes=x))); cs[cs!=x]})

house_makeup <- lapply(2:nrow(occupancy_data)-1,function(x)as.numeric(c(occupancy_data$children[occupancy_data$type==x],
                                           occupancy_data$adults[occupancy_data$type==x],
                                           occupancy_data$elderly[occupancy_data$type==x])))

demographic_index <- rep(0,length(V(new_g)))

for(i in 1:number_of_households){
  residents <- which(hh_labels==i)
  occupants <- house_makeup[[household_types[i]]]
  labels <- rep(1:3,times=occupants)
  demographic_index[residents] <- labels
}
demographic_index <<- demographic_index
min_age <- c(0,10,20,30,40,50,60,70,80)
max_age <- c(9,19,29,39,49,59,69,79,150)
cfr <- c(0.002,0.006,0.03,0.08,0.15,0.6,2.2,5.1,9.3)/100
grouped_cfr <<- c(mean(cfr[1:2]),mean(cfr[3:7]),mean(cfr[8:9]))
tune <<- grouped_cfr/max(grouped_cfr)

paste0(sum(demographic_index==1),' children')
paste0(sum(demographic_index==2),' adults')
paste0(sum(demographic_index==3),' elderly')


## things ignoring:
# children, elderly, not working, communal/institutional living, social
## build workplace connections of 20 ish

worker_index <- rep(0,length(demographic_index))
worker_index[demographic_index==2] <- 1
worker_index[demographic_index==3&runif(length(worker_index))<0.2] <- 1
n_adults <- sum(worker_index)
# https://www.hse.gov.uk/contact/faqs/toilets.htm
number_workplaces <- rpois(1,n_adults/15)
workplace_index <- rep(0,length(V(new_g)))
workplace_index[worker_index==1] <- sample(1:number_workplaces,n_adults,replace=T)

## add work connections
for(i in 1:number_workplaces) {
  workers <- which(workplace_index==i)
  if(length(workers)>1)
    for(j in 2:length(workers))
      for(k in 1:(j-1)){
        new_g <- add_edges(new_g,edges=c(workers[j],workers[k]))
      }
}

#pdf('figures/contactnetwork.pdf')
#par(mar=c(1,1,1,1))
#plot(induced.subgraph(new_g,1:pts),layout=save_layout,
#     vertex.label=NA,vertex.size=10,vertex.color=c('grey','hotpink','navyblue')[demographic_index[1:pts]])
#legend('topleft',legend=c('<19','19-65','>65'),pt.cex=2,col='black',pch=21, pt.bg=c('grey','hotpink','navyblue'),bty='n')
#dev.off()

#plot.igraph(new_g,vertex.label=NA,vertex.size=1,layout=save_layout)
#cluster_sizes <- sapply(V(new_g),function(x)ego_size(new_g,order=2,nodes=x))
#hist(cluster_sizes,main='',xlab='Cluster size')
#c(mean(cluster_sizes),quantile(cluster_sizes,c(0.25,0.5,0.75)))

# plot degree distribution - aiming for mean=17.5
degreedistribution <- degree.distribution(new_g)*length(E(new_g))
#barplot(degreedistribution,ylab='Number of people', xlab='Number of connections',names.arg=0:(length(degreedistribution)-1),main='')
average_contacts <- sum(degreedistribution*c(1:length(degreedistribution)-1)/length(E(new_g)))
length(E(new_g))/length(V(new_g))*2

contact_list <<- lapply(V(new_g),function(x) {cs <- as.vector(unlist(ego(new_g,order=1,nodes=x))); cs[cs!=x]})

contact_of_contact_list <<- lapply(V(new_g),function(x) {
  cs <- as.vector(unlist(ego(new_g,order=1,nodes=x))); 
  cofcs <- as.vector(unlist(ego(new_g,order=2,nodes=x))); 
  ccs <- funique(c(cs,cofcs))
  ccs[ccs!=x]
  })

## generate random edges network for random transmission
random_g <- sample_gnp(length(V(new_g)), 10/length(V(new_g)))
random_list <<- lapply(V(random_g),function(x) {cs <- as.vector(unlist(ego(random_g,order=1,nodes=x))); cs[cs!=x]})
