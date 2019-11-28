setwd('../ADAGIO/')
library(RColorBrewer)
library(plotrix)

pval_function <- function(success0,n0,success1,n1){
  p0 <- success0/n0
  p1 <- success1/n1
  sigma0 <- p0 * ( 1 - p0 ) /n0
  sigma1 <- p1 * ( 1 - p1 ) /n1
  zval <- (p1-p0)/(sqrt(sigma0+sigma1))
  dnorm(zval)
}

## power ##################################################

VE <- 0.6
rr <- 1 - VE

n0 <- seq(100,1000,by=10)
events <- seq(1,100,by=1)
pvals_list <- list()
nsim <- 1000
for(i in 1:length(n0)){
  pvals_list[[i]] <- matrix(0,nrow=nsim,ncol=length(events))
  n1 <- n0[i]
  for(j in 1:length(events)){
    total_events <- events[j]
    for(k in 1:nsim){
      n0_events <- rbinom(1,total_events,1/(1+rr))
      n1_events <- total_events - n0_events
      pvals_list[[i]][k,j] <- pval_function(success0=n0[i]-n0_events,n0=n0[i],success1=n1-n1_events,n1=n1)
    }
  }
}
pvals <- sapply(pvals_list,function(x)colSums(x<0.05))/nsim
pvals <- pvals[nrow(pvals):1,]
get.pal=colorRampPalette(brewer.pal(9,"RdBu"))
redCol=rev(get.pal(9))
bkT <- seq(max(pvals)+1e-10, 0,length=length(redCol)+1)
cex.lab <- 1.5
maxval <- round(bkT[1],digits=1)
col.labels<- c(0,maxval/2,maxval)
cellcolors <- vector()
for(ii in 1:length(unlist(pvals)))
  cellcolors[ii] <- redCol[tail(which(unlist(pvals[ii])<bkT),n=1)]
pdf('power.pdf')
color2D.matplot(pvals,cellcolors=cellcolors,main="",xlab="Sample size",ylab="Events",cex.lab=1,axes=F,border=NA)
fullaxis(side=2,las=1,at=1:nrow(pvals),labels=events,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
fullaxis(side=1,las=2,at=1:ncol(pvals),labels=n0*2,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
color.legend(ncol(pvals)+0.5,0,ncol(pvals)+2,nrow(pvals),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
dev.off()


## T1E ##################################################

VE <- 0
rr <- 1 - VE

n0 <- seq(100,1000,by=10)
events <- seq(1,100,by=1)
pvals_list <- list()
nsim <- 1000
for(i in 1:length(n0)){
  pvals_list[[i]] <- matrix(0,nrow=nsim,ncol=length(events))
  n1 <- n0[i]
  for(j in 1:length(events)){
    total_events <- events[j]
    for(k in 1:nsim){
      n0_events <- rbinom(1,total_events,1/(1+rr))
      n1_events <- total_events - n0_events
      pvals_list[[i]][k,j] <- pval_function(success0=n0[i]-n0_events,n0=n0[i],success1=n1-n1_events,n1=n1)
    }
  }
}
pvals <- sapply(pvals_list,function(x)colSums(x<0.05))/nsim
pvals <- pvals[nrow(pvals):1,]
get.pal=colorRampPalette(brewer.pal(9,"RdBu"))
redCol=rev(get.pal(5))
bkT <- seq(max(pvals)+1e-10, 0,length=length(redCol)+1)
cex.lab <- 1.5
maxval <- round(bkT[1],digits=1)
col.labels<- c(0,maxval/2,maxval)
cellcolors <- vector()
for(ii in 1:length(unlist(pvals)))
  cellcolors[ii] <- redCol[tail(which(unlist(pvals[ii])<bkT),n=1)]
pdf('type1error.pdf')
color2D.matplot(pvals,cellcolors=cellcolors,main="",xlab="Sample size",ylab="Events",cex.lab=1,axes=F,border=NA)
fullaxis(side=2,las=1,at=1:nrow(pvals),labels=events,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
fullaxis(side=1,las=2,at=1:ncol(pvals),labels=n0*2,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
color.legend(ncol(pvals)+0.5,0,ncol(pvals)+2,nrow(pvals),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
dev.off()


## power, gs ##################################################
VE <- 0.6
rr <- 1 - VE

n0 <- seq(200,900,by=10)
events <- seq(1,100,by=1)
pvals_list <- pvals_list1 <- pvals_list2 <- list()
nsim <- 1000
for(i in 1:length(n0)){
  pvals_list[[i]] <- pvals_list1[[i]] <- pvals_list2[[i]] <- matrix(0,nrow=nsim,ncol=length(events))
  n1 <- n0[i]
  for(j in 1:length(events)){
    total_events <- events[j]
    for(k in 1:nsim){
      n0_events <- rbinom(1,total_events,1/(1+rr))
      n1_events <- total_events - n0_events
      n0_events_early <- rbinom(1,n0_events,n0[i]/1000)
      n1_events_early <- rbinom(1,n1_events,n1/1000)
      pval1 <- ifelse(n0_events_early+n1_events_early==0,0.5,pval_function(success0=n0[i]-n0_events_early,n0=n0[i],success1=n1-n1_events_early,n1=n1))
      pval2 <- ifelse(n0_events+n1_events==0,0.5,pval_function(success0=1000-n0_events,n0=1000,success1=1000-n1_events,n1=1000))
      pvals_list[[i]][k,j] <- min(pval1,pval2)
      pvals_list1[[i]][k,j] <- pval1
      pvals_list2[[i]][k,j] <- pval2
    }
  }
}
pvals <- sapply(pvals_list1,function(x)colSums(x<0.025))/nsim
pvals <- pvals[nrow(pvals):1,]
get.pal=colorRampPalette(brewer.pal(9,"RdBu"))
redCol=rev(get.pal(9))
bkT <- seq(max(pvals)+1e-10, 0,length=length(redCol)+1)
cex.lab <- 1.5
maxval <- round(bkT[1],digits=1)
col.labels<- c(0,maxval/2,maxval)
cellcolors <- vector()
for(ii in 1:length(unlist(pvals)))
  cellcolors[ii] <- redCol[tail(which(unlist(pvals[ii])<bkT),n=1)]
pdf('power2.pdf')
color2D.matplot(pvals,cellcolors=cellcolors,main="",xlab="Information time",ylab="Total events (sample size=2000)",cex.lab=1,axes=F,border=NA)
fullaxis(side=2,las=1,at=1:nrow(pvals),labels=events,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
fullaxis(side=1,las=2,at=1:ncol(pvals),labels=n0*2,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
color.legend(ncol(pvals)+0.5,0,ncol(pvals)+2,nrow(pvals),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
dev.off()


## T1E, gs ##################################################

VE <- 0
rr <- 1 - VE

n0 <- seq(200,900,by=10)
events <- seq(1,100,by=1)
pvals_list <- pvals_list1 <- pvals_list2 <- list()
nsim <- 1000
for(i in 1:length(n0)){
  pvals_list[[i]] <- pvals_list1[[i]] <- pvals_list2[[i]] <- matrix(0,nrow=nsim,ncol=length(events))
  n1 <- n0[i]
  for(j in 1:length(events)){
    total_events <- events[j]
    for(k in 1:nsim){
      n0_events <- rbinom(1,total_events,1/(1+rr))
      n1_events <- total_events - n0_events
      n0_events_early <- rbinom(1,n0_events,n0[i]/1000)
      n1_events_early <- rbinom(1,n1_events,n1/1000)
      pval1 <- ifelse(n0_events_early+n1_events_early==0,0.5,pval_function(success0=n0[i]-n0_events_early,n0=n0[i],success1=n1-n1_events_early,n1=n1))
      pval2 <- ifelse(n0_events+n1_events==0,0.5,pval_function(success0=1000-n0_events,n0=1000,success1=1000-n1_events,n1=1000))
      pvals_list[[i]][k,j] <- min(pval1,pval2)
      pvals_list1[[i]][k,j] <- pval1
      pvals_list2[[i]][k,j] <- pval2
    }
  }
}
pvals <- sapply(pvals_list2,function(x)colSums(x<0.025))/nsim
pvals <- pvals[nrow(pvals):1,]
get.pal=colorRampPalette(brewer.pal(9,"RdBu"))
redCol=rev(get.pal(5))
bkT <- seq(max(pvals)+1e-10, 0,length=length(redCol)+1)
cex.lab <- 1.5
maxval <- round(bkT[1],digits=1)
col.labels<- c(0,maxval/2,maxval)
cellcolors <- vector()
for(ii in 1:length(unlist(pvals)))
  cellcolors[ii] <- redCol[tail(which(unlist(pvals[ii])<bkT),n=1)]
pdf('type1error2.pdf')
color2D.matplot(pvals,cellcolors=cellcolors,main="",xlab="Information time",ylab="Total events (sample size=2000)",cex.lab=1,axes=F,border=NA)
fullaxis(side=2,las=1,at=1:nrow(pvals),labels=events,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
fullaxis(side=1,las=2,at=1:ncol(pvals),labels=n0*2,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
color.legend(ncol(pvals)+0.5,0,ncol(pvals)+2,nrow(pvals),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
dev.off()
