library(pracma)

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
alphas <- c(0.02,0.03,0.04,0.05)
for(metric in c('weight')){ # ,'exposure'
  for(type in c('power')){
    par_results <- readRDS(paste0('storage/',type,'_results.Rds'))
    zvals <- par_results[,1]
    ycol <- which(colnames(par_results)==paste0('case_',metric))
    xcol <- which(colnames(par_results)==paste0('',metric))
    ys <- par_results[,ycol]
    xs <- par_results[,xcol]
    #par_results <- ttemat <- tte <- c()
    zvals[is.na(zvals)] <- 0
    
    binx <- equal_freq_disc(xs,nbins=10)
    x_binned <- binx[[1]]
    x_breaks <- binx[[2]]
    biny <- equal_freq_disc(ys,nbins=10)
    y_binned <- biny[[1]]
    y_breaks <- biny[[2]]
    x_points <- length(x_breaks)-1
    y_points <- length(y_breaks)-1
    for(alpha in alphas){
      grid_pval <- grid_p <- matrix(NA,nrow=x_points,ncol=y_points)
      for(i in 1:x_points){
        for(j in 1:y_points){
          grid_zvals <- zvals[x_binned==i&y_binned==j]
          if(length(grid_zvals)>10){
            grid_pval[j,i] <- sum(grid_zvals>qnorm(1-alpha))/length(grid_zvals)
            grid_p[j,i] <- length(grid_zvals)
          }
        }
      }
      
      yvals <- sapply(2:length(y_breaks),function(x)mean(y_breaks[x:(x-1)]))
      xvals <- sapply(2:length(x_breaks),function(x)mean(x_breaks[x:(x-1)]))
      # x11(width=15); par(mfrow=c(1,2)); 
      # matplot(t(repmat(xvals,5,1)),t(grid_pval[seq(1,10,by=2),]),frame=F,xaxt='n',ylim=0:1,typ='l',xlab='',xlim=range(x_breaks))
      # fullaxis(side=1,las=2,at=x_breaks,labels=round(x_breaks),line=-1,pos=NA,outer=FALSE,font=NA,lwd=1,cex.axis=1.25,hadj=1)
      # mtext(side=1,text='Total weighted sample size',cex=1.5,line=3)
      # matplot(t(repmat(yvals,5,1)),grid_pval[,seq(1,10,by=2)],frame=F,xaxt='n',ylim=0:1,typ='l',xlab='',xlim=range(y_breaks))
      # fullaxis(side=1,las=2,at=y_breaks,labels=round(y_breaks),line=-1,pos=NA,outer=FALSE,font=NA,lwd=1,cex.axis=1.25,hadj=1)
      # mtext(side=1,text='Total weighted number of confirmed cases',cex=1.5,line=3)
      
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
      pdf(paste0('figures/colour',type,metric,100*alpha,'.pdf')); par(mar=c(6,6,2,6),adj=0)
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
      
      pdf(paste0('figures/line',type,metric,100*alpha,'.pdf'));  
      par(mar=c(5,5,2,2))
      threshold <- yvals[apply(grid_pval,2,function(x)which(x>0.8)[1])]
      y2 <- predict(lm(threshold~xvals))
      plot(xvals[!is.na(threshold)],y2,xaxt='n',yaxt='n',frame=F,lty=3,typ='l',lwd=2,ylim=range(y_breaks),xlim=range(x_breaks),xlab='',ylab='')
      threshold <- yvals[apply(grid_pval,2,function(x)which(x>0.85)[1])]
      for(i in seq(2,length(x_breaks)-1,by=2)) {
        abline(h=y_breaks[i],lwd=2,col='grey')
        abline(v=x_breaks[i],lwd=2,col='grey')
      }
      y2 <- predict(lm(threshold~xvals))
      lines(xvals[!is.na(threshold)],y2,typ='l',lwd=2,lty=2)
      threshold <- yvals[apply(grid_pval,2,function(x)which(x>0.9)[1])]
      if(alpha>0.04) {
        y2 <- predict(lm(threshold~xvals))
        lines(xvals[!is.na(threshold)],y2,typ='l',lwd=2,lty=1)
        legend(x=max(x_breaks[-length(x_breaks)]),y=max(y_breaks),legend=c(0.9,0.85,0.8),lwd=2,lty=c(1,2,3),bty='n',cex=1.25)
      }else{
        legend(x=max(x_breaks[-length(x_breaks)]),y=max(y_breaks),legend=c(0.85,0.8),lwd=2,lty=c(2,3),bty='n',cex=1.25)
      }
      fullaxis(side=2,las=1,at=y_breaks,labels=round(y_breaks),line=-1,pos=NA,outer=FALSE,font=NA,lwd=1,cex.axis=1.25,hadj=1)
      fullaxis(side=1,las=2,at=x_breaks,labels=round(x_breaks),line=-1,pos=NA,outer=FALSE,font=NA,lwd=1,cex.axis=1.25,hadj=1)
      mtext(side=1,text='Total weighted sample size',cex=1.5,line=3)
      mtext(side=2,text='Total weighted number of confirmed cases',cex=1.5,line=3)
      mtext(side=3,text=TeX(string),cex=1.5)
      dev.off()
    }
  }
}
