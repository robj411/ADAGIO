source('set_up_script.R')
registerDoParallel(cores=32)
## saves to storage/cl* and storage/res*
rm(g)
rm(new_g)
rm(random_g)
rm(hh)
## ring vaccination trial ##################################################

nTrials <- 1000
vaccine_efficacies <- c(0.7)
adaptations <- c('Ney','Ros','TST','TS','')
cluster_flags <- 0
trial_designs <- expand.grid(VE=vaccine_efficacies,cluster=cluster_flags,adapt=adaptations,stringsAsFactors = F)
trial_designs$weight <- 'continuous'
nComb <- sum(trial_designs$adapt=='')
nCombAdapt <- nComb*length(adaptations)
trial_designs$powertst <- trial_designs$VE_esttst <- trial_designs$VE_sdtst <- trial_designs$VE_estht <- trial_designs$VE_sdht <- 
  trial_designs$power <- trial_designs$VE_est <- trial_designs$VE_sd <- trial_designs$vaccinated <- trial_designs$infectious <- trial_designs$enrolled <- 0
ref_recruit_day <- 30
eval_days <- c(21,23,25,27,29)

nClusters <- 160

cls <- seq(60,nClusters,by=10)
for(eval_day in eval_days){
  
  #eval_day <- 31
  latest_infector_time <- eval_day - 0
  
  for(des in 1:5){
    set.seed(des)
    cluster_flag <<- trial_designs$cluster[des]
    direct_VE <<- trial_designs$VE[des]
    adaptation <<- trial_designs$adapt[des]
    filename <- paste0('storage/res',eval_day,des,'.Rds')
    if(file.exists(filename)){
      res <- readRDS(filename)
    }else{
      res <- foreach(tr = 1:nTrials) %dopar% {
        vaccinated_count <- infectious_count <- enrolled_count <- 0
        randomisation_ratios <- c()
        people_per_ratio <- c()
        vaccinees <- trial_participants <- c()
        infectious_by_vaccine <- excluded <- matrix(0,nrow=nClusters,ncol=2)
        results_list <- list()
        allocation_ratio <- 0.5
        netwk_list <- list()
        for(iter in 1:nClusters){
          ## select random person to start
          randomisation_ratios[iter] <- allocation_ratio
          first_infected <- sample(g_name[eligible_first_person],1)
          netwk <- simulate_contact_network(first_infected,cluster_flag=cluster_flag,end_time=eval_day,allocation_ratio=allocation_ratio,direct_VE=direct_VE,individual_recruitment_times=T,spread_wrapper=covid_spread_wrapper)
          netwk[[9]] <- allocation_ratio
          results_list[[iter]] <- netwk[[1]]
          
          vaccinees[iter] <- netwk[[4]]
          trial_participants[iter] <- netwk[[5]]
          
          ## iter corresponds to a day, so we can adapt the enrollment rate on iter=31
          if(adaptation!=''&&iter %% eval_day == 0 && sum(vaccinees)>0){
            probs <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,max_time=length(results_list),contact_network=-1,observed=observed)
            pop_sizes2 <- probs[[2]]
            fails <- probs[[3]]
            allocation_ratio <- response_adapt(fails,pop_sizes2,days=iter,adaptation)
            people_per_ratio <- rbind(people_per_ratio,c(sum(trial_participants),iter,allocation_ratio))
            #if(allocation_ratio==0) break
          }
          netwk[[10]] <- people_per_ratio
          netwk_list[[iter]] <- netwk
        }
        return(netwk_list)
      }
      print(c(des,sum(sapply(res,length)==nClusters)))
      #res_list[[des]] <- res
      saveRDS(res,filename)
    }
  }
  


  for(i in 1:length(cls)){
    cl <- cls[i]
    filename <- paste0('storage/cl',eval_day,cl,'.Rds')
    if(!file.exists(filename)){
      power <- vax <- ss <- vest <- rep(0,nCombAdapt)
      for(des in 1:5){
        set.seed(des)
        res <- readRDS(paste0('storage/res',eval_day,des,'.Rds'))
        enrolled_count <- vaccinated_count <- 0
        result_mat <- matrix(0,nrow=nTrials,ncol=4)
        adaptation <<- trial_designs$adapt[des]
        result_mat <- foreach(tr = 1:nTrials,.combine=rbind)%dopar%{
          set.seed(tr)
          netwk_list <- res[[tr]][1:cl]
          vaccinees <- sapply(netwk_list,function(netwk)netwk[[4]])
          trial_participants <- sapply(netwk_list,function(netwk)netwk[[5]])
          results_list <- lapply(netwk_list,function(x)x[[1]])
          ## regular test
          threshold <- qnorm(0.95)
          eval_list <- get_efficacious_probabilities(results_list,vaccinees,trial_participants,tested=F,contact_network=-1,observed=observed)
          zval  <- calculate_zval(eval_list[[3]],eval_list[[2]])
          vest <- eval_list[[1]]
          ## correcting for trend 
          if(adaptation!=''&eval_day<cl){
            randomisation_ratios <- sapply(netwk_list,function(netwk)netwk[[9]])
            people_per_ratio <- netwk_list[[cl]][[10]]
            threshold <- trend_robust_function(results_list,vaccinees,trial_participants,contact_network=-1,
                                               tested=F,randomisation_ratios=randomisation_ratios,adaptation=adaptation,people_per_ratio=people_per_ratio,observed=observed)
          }
          c(zval,threshold,sum(vaccinees),sum(trial_participants),vest)
          #result_mat[tr,1] <- pval
          #result_mat[tr,2] <- threshold
          #result_mat[tr,3] <- sum(vaccinees)
          #result_mat[tr,4] <- sum(trial_participants)
        }
        saveRDS(result_mat,paste0('storage/resultscl',cl,'des',des,'.Rds'))
        power[des] <- sum(result_mat[,1]>result_mat[,2],na.rm=T)/sum(!is.na(result_mat[,1])&!is.na(result_mat[,2]))
        vax[des] <- mean(result_mat[,3],na.rm=T)
        ss[des] <- mean(result_mat[,4],na.rm=T)
        vest[des] <- mean(result_mat[,5],na.rm=T)
      }
      print(c(cl,power))
      saveRDS(list(power,vax,ss,vest),filename)
    }
  }
  rm(res)

  powers <- vax <- ss <- vest <- matrix(0,nrow=5,ncol=length(cls))
  for(i in 1:length(cls)){
    cl <- cls[i]
    lst <- readRDS(paste0('storage/cl',eval_day,cl,'.Rds'))
    powers[,i] <- lst[[1]]
    vax[,i] <- lst[[2]]
    ss[,i] <- lst[[3]]
    vest[,i] <- lst[[4]]
  }
  cont <- ss - vax
  cols <- rainbow(4)
  
  pdf(paste0('figures/sspower',eval_day,'.pdf'))
  par(mar=c(5,5,2,2))
  ind <- which(!duplicated(powers[5,]))
  plot(ss[5,ind],powers[5,ind],typ='l',lwd=2,ylim=c(0,1),xlim=c(0,max(ss)),frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Sample size',ylab='Power',main=paste0('Adapt every ',eval_day,' index cases'))
  adapt_days <- floor(cl/eval_day)
  for(ad in 1:adapt_days)
    abline(v=max(ss)/cl*ad*eval_day,col='grey',lwd=2,lty=2)
  for(h in seq(0,1,by=0.1)) abline(h=h,col='grey',lwd=2,lty=1)
  for(i in 1:4) {
    ind <- which(!duplicated(powers[i,]))
    lines(ss[i,ind],powers[i,ind],col=cols[i],lwd=2)
  }
  legend(bty='n',x=0,y=1,cex=1.25,col=c('black',cols),lwd=2,lty=1,legend=c('iRCT','Ney','Ros','TST','TS'))
  dev.off()
  
  pdf(paste0('figures/vaxpower',eval_day,'.pdf'))
  par(mar=c(5,5,2,2))
  ind <- which(!duplicated(powers[5,]))
  plot(vax[5,ind],powers[5,ind],typ='l',lwd=2,ylim=c(0,1),xlim=c(0,max(vax)),frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Vaccinated',ylab='Power',main=paste0('Adapt every ',eval_day,' index cases'))
  for(h in seq(0,1,by=0.1)) abline(h=h,col='grey',lwd=2,lty=1)
  cols <- rainbow(4)
  for(i in 1:4) {
    ind <- which(!duplicated(powers[i,]))
    lines(vax[i,ind],powers[i,ind],col=cols[i],lwd=2)
  }
  legend(bty='n',x=0,y=1,cex=1.25,col=c('black',cols),lwd=2,lty=1,legend=c('iRCT','Ney','Ros','TST','TS'))
  dev.off()
  
  pdf(paste0('figures/contpower',eval_day,'.pdf'))
  par(mar=c(5,5,2,2))
  ind <- which(!duplicated(powers[5,]))
  plot(cont[5,ind],powers[5,ind],typ='l',lwd=2,ylim=c(0,1),xlim=c(0,max(cont)),frame=F,cex.axis=1.5,cex.lab=1.5,xlab='Controls',ylab='Power',main=paste0('Adapt every ',eval_day,' index cases'))
  for(h in seq(0,1,by=0.1)) abline(h=h,col='grey',lwd=2,lty=1)
  for(i in 1:4) {
    ind <- which(!duplicated(powers[i,]))
    lines(cont[i,ind],powers[i,ind],col=cols[i],lwd=2)
  }
  legend(bty='n',x=0,y=1,cex=1.25,col=c('black',cols),lwd=2,lty=1,legend=c('iRCT','Ney','Ros','TST','TS'))
  dev.off()
  
  cols <- c(cols,'grey')
  people_per_day <- mean(apply(ss,1,function(x)x/cls))
  first_over_80 <- apply(powers,1,function(x){
    interpolated <- approx(cls,x,xout=seq(min(cls),max(cls)))$y
    which(interpolated>=0.8)[1]
    })
  first_day_over_80 <- min(cls) - 1 + first_over_80
  latest_day <- max(first_day_over_80,na.rm=T)
  days_left_over <- latest_day - first_day_over_80
  vax_over_80 <- sapply(1:nrow(vax),function(x){
    interpolated <- approx(cls,vax[x,],xout=seq(min(cls),max(cls)))$y
    interpolated[first_over_80[x]]
  })
  cont_over_80 <- sapply(1:nrow(vax),function(x){
    interpolated <- approx(cls,cont[x,],xout=seq(min(cls),max(cls)))$y
    interpolated[first_over_80[x]]
  })
  vest_over_80 <- sapply(1:nrow(vax),function(x){
    interpolated <- approx(cls,vest[x,],xout=seq(min(cls),max(cls)))$y
    interpolated[first_over_80[x]]
  })
  pdf(paste0('figures/contvax',eval_day,'.pdf'))
  #x11(); 
  par(mar=c(5,5,2,2))
  plot(x=range(c(vax_over_80,cont_over_80),na.rm=T),y=range(c(cont_over_80,vax_over_80),na.rm=T),col='grey',typ='l',frame=F,
       xlab='Controls',ylab='Vaccinees',cex.lab=1.5,cex.axis=1.5,lwd=2)
  maxss <- max(c(vax_over_80,cont_over_80),na.rm=T)
  minss <- min(c(vax_over_80,cont_over_80),na.rm=T)
  text(label=paste0('Sample size = ',round(maxss+minss)),x=(minss+maxss+500)/2,y=(minss+maxss-500)/2,srt=-45,cex=1.25,pos=3)
  for(i in 1:5){
    ss_tmp <- maxss - 150*i
    lines(x=c(minss,ss_tmp),y=rev(c(minss,ss_tmp)),col='grey',lwd=2,lty=3)
    text(label=paste0(round(minss+ss_tmp)),x=(minss+ss_tmp+250)/2,y=(minss+ss_tmp-250)/2,srt=-45,cex=1.25,pos=3,col='grey')
  }
  for(i in 1:length(first_over_80))
    if(!is.na(first_over_80[i])){
      points(cont_over_80[i],vax_over_80[i],col=cols[i],pch=16,cex=2)
      if(days_left_over[i]>0)
      arrows(x0=cont_over_80[i],y0=vax_over_80[i],y1=vax_over_80[i]+days_left_over[i]*people_per_day,lwd=2,col=cols[i],lty=2)
    }
  legend(bty='n',x=max(c(vax_over_80,cont_over_80),na.rm=T)*0.9,y=mean(c(vax_over_80,cont_over_80),na.rm=T)*1.3,
         cex=1.5,col=cols,legend=c('Ney','Ros','TST','TS','iRCT'),pch=16)  
  lines(x=range(c(vax_over_80,cont_over_80),na.rm=T),y=rev(range(c(vax_over_80,cont_over_80),na.rm=T)),
        col='grey',lwd=2,lty=2)
  text(label='Fixed equal allocation',x=2100,y=2100,srt=45,cex=1.25,pos=3)
  
  dev.off()

  networks <- round((vax_over_80+cont_over_80)/people_per_day)
  netdiff <- max(networks,na.rm=T)-networks
  netdiffrange <- range(netdiff[netdiff>0],na.rm=T)
  round(netdiffrange*people_per_day)
  print(networks)
  print(vest_over_80)
  
  print(sort(sapply(ls(),function(x)object.size(get(x)))))
  
  ad <- c('Ney','Ros','TST','','FR')
  {
  for(i in c(1,2,3,5)){
    cat(paste0("\\newcommand{\\",ad[i],"networks}{",round(networks[i]),"\\xspace}\n"))
    cat(paste0("\\newcommand{\\",ad[i],"SS}{",round((vax_over_80+cont_over_80)[i],-2),"\\xspace}\n"))
    cat(paste0("\\newcommand{\\",ad[i],"vax}{",round((vax_over_80)[i],-2),"\\xspace}\n"))
    cat(paste0("\\newcommand{\\",ad[i],"cont}{",round((cont_over_80)[i],-2),"\\xspace}\n"))
  }
  cat(paste0("\\newcommand{\\ppday}{",round(people_per_day),"\\xspace}\n"))
  for(i in c(1,2)){
    cat(paste0("\\newcommand{\\",c('min','max')[i],"day}{",round(netdiffrange[i]),"\\xspace}\n"))
    cat(paste0("\\newcommand{\\",c('min','max')[i],"ppl}{",round(netdiffrange[i]*people_per_day),"\\xspace}\n"))
  }
  }
}





