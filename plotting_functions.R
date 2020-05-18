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
           col = matrix(cellcolors,nrow=x_points,ncol=y_points,byrow = T)[i,j], border = NA)
    }
  }
}

library(RColorBrewer)
library(infotheo)
library(plotrix)
library(latex2exp)
alphas <- c(0.01,0.02,0.03,0.04,0.05)
for(metric in c('weight','exposure')){
  for(type in c('power')){
    par_results <- readRDS(paste0('storage/',type,'_results.Rds'))
    pvals <- par_results[,1]
    ycol <- which(colnames(par_results)==paste0('case_',metric))
    xcol <- which(colnames(par_results)==paste0('',metric))
    ys <- par_results[,ycol]
    xs <- par_results[,xcol]
    #par_results <- ttemat <- tte <- c()
    pvals[is.na(pvals)] <- 1
    
    binx <- equal_freq_disc(xs,nbins=10)
    x_binned <- binx[[1]]
    x_breaks <- binx[[2]]
    biny <- equal_freq_disc(ys,nbins=10)
    y_binned <- biny[[1]]
    y_breaks <- biny[[2]]
    x_points <- length(x_breaks)-1
    y_points <- length(y_breaks)-1
    for(alpha in alphas){
    grid_pval <- matrix(NA,nrow=x_points,ncol=y_points)
    for(i in 1:x_points){
      for(j in 1:y_points){
        grid_pvals <- pvals[x_binned==i&y_binned==j]
        if(length(grid_pvals)>10){
          grid_pval[j,i] <- sum(grid_pvals<alpha)/length(grid_pvals)
        }
      }
    }
    
    #grid_pval <- grid_pval[nrow(grid_pval):1,]
    get.pal=colorRampPalette(brewer.pal(9,"Spectral"))
    redCol=rev(get.pal(10))
    bkT <- seq(1+1e-10, 0-1e-10,length=length(redCol)+1)
    cex.lab <- 1.5
    maxval <- round(bkT[1],digits=1)
    col.labels<- c(0,1/2,1)
    cellcolors <- vector()
    for(ii in 1:length(unlist(grid_pval)))
      if(!is.na(grid_pval[ii]))
        cellcolors[ii] <- redCol[tail(which(unlist(grid_pval[ii])<bkT),n=1)]
    pdf(paste0('figures/ph',type,metric,100*alpha,'.pdf')); par(mar=c(6,6,2,6),adj=0)
    #color2D.matplot(grid_pval,cellcolors=cellcolors,x_breaks,y_breaks)
    #x11()
    plot_rect(grid_pval,cellcolors=cellcolors,x_breaks,y_breaks,x_points,y_points)
    fullaxis(side=2,las=1,at=y_breaks,labels=round(y_breaks),line=-1,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.25,hadj=1)
    fullaxis(side=1,las=2,at=x_breaks,labels=round(x_breaks),line=-1,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1.25,hadj=1)
    color.legend(max(x_breaks)*1.01,0,max(x_breaks)*1.03,max(y_breaks),col.labels,rev(redCol),gradient="y",cex=1.25,align="rb")
    mtext(side=1,text='Total weighted sample size',cex=1.5,line=3)
    mtext(side=2,text='Total weighted number of confirmed cases',cex=1.5,line=3)
    string <- paste0('$\\alpha = ',alpha,'$')
    mtext(side=3,text=TeX(string),cex=1.5)
    dev.off()
    }
  }
}
