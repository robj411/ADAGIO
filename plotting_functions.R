equal_freq_disc <- function(x,bins=NULL,nbins = NROW(x)^(1/3)){
  if(is.null(bins)){
    sorted_x <- sort(x)
    bins <- quantile(x,seq(0,1,length=nbins+1))
    bins[1] <- bins[1] - 1e-6
    bins[nbins] <- bins[nbins] + 1e-6
  }
  cut_x <- cut(x,bins,labels=1:nbins)
  return(list(cut_x,bins))
}

plot_rect <- function(grid_pval,cellcolors=cellcolors,x_breaks,y_breaks,x_points,y_points){
  par(mfrow=c(1,1),mar=c(c(6,6,2,6)))
  plot(0,0,col='white',xlim=range(x_breaks),ylim=range(y_breaks),axes=F,xlab='',ylab='')
  for(i in 1:x_points){
    for(j in 1:y_points){
      rect(xleft=x_breaks[i], ybottom=y_breaks[j], xright=x_breaks[i+1], ytop=y_breaks[j+1], density = NULL, angle = 45,
           col = matrix(cellcolors,nrow=x_points,ncol=y_points,byrow = T)[i,j], border = NULL)
    }
  }
}

library(RColorBrewer)
library(infotheo)
library(plotrix)

for(type in c('t1e','power')){
  par_results <- readRDS(paste0(type,'_results.Rds'))
  for(metric in c('weight','exposure')){
    pvals <- par_results[,1]
    ycol <- which(colnames(par_results)==paste0('case_',metric))
    xcol <- which(colnames(par_results)==paste0('',metric))
    ys <- par_results[,ycol]
    xs <- par_results[,xcol]
    #par_results <- ttemat <- tte <- c()
    pvals[is.na(pvals)] <- 1
    
    binx <- equal_freq_disc(xs,nbins=15)
    x_binned <- binx[[1]]
    x_breaks <- binx[[2]]
    biny <- equal_freq_disc(ys,nbins=15)
    y_binned <- biny[[1]]
    y_breaks <- biny[[2]]
    x_points <- length(x_breaks)-1
    y_points <- length(y_breaks)-1
    grid_pval <- matrix(NA,nrow=x_points,ncol=y_points)
    for(i in 1:x_points){
      for(j in 1:y_points){
        grid_pvals <- pvals[x_binned==i&y_binned==j]
        if(length(grid_pvals)>0){
          grid_pval[i,j] <- sum(grid_pvals<0.05)/length(grid_pvals)
        }
      }
    }
    
    #grid_pval <- grid_pval[nrow(grid_pval):1,]
    get.pal=colorRampPalette(brewer.pal(9,"Spectral"))
    redCol=rev(get.pal(10))
    bkT <- seq(max(grid_pval,na.rm=T)+1e-10, min(grid_pval,na.rm=T)-1e-10,length=length(redCol)+1)
    cex.lab <- 1.5
    maxval <- round(bkT[1],digits=1)
    col.labels<- c(0,maxval/2,maxval)
    cellcolors <- vector()
    for(ii in 1:length(unlist(grid_pval)))
      if(!is.na(grid_pval[ii]))
        cellcolors[ii] <- redCol[tail(which(unlist(grid_pval[ii])<bkT),n=1)]
    pdf(paste0('ph',type,metric,'.pdf')); par(mar=c(6,6,2,6))
    #color2D.matplot(grid_pval,cellcolors=cellcolors,x_breaks,y_breaks)
    #x11()
    plot_rect(grid_pval,cellcolors=cellcolors,x_breaks,y_breaks,x_points,y_points)
    fullaxis(side=2,las=1,at=y_breaks,labels=round(y_breaks),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
    fullaxis(side=1,las=2,at=x_breaks,labels=round(x_breaks),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
    color.legend(max(x_breaks)*1.01,0,max(x_breaks)*1.03,max(y_breaks),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
    dev.off()
    
  }
}