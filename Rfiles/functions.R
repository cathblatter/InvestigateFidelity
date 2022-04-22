as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

bias<-function(vec, theta){
  return(mean(vec)-theta)
}

performanceMeas<-function(vecEst, vecSe, vec.p, theta, alpha=0.05){
  
  z<-qnorm(p=1-alpha/2)
  CI<-data.frame(CI_low=vecEst-z*vecSe, CI_up=vecEst+z*vecSe)
  
  res<-c(
    AvEst=mean(vecEst),
    EmpSE=sd(vecEst),
    Bias=bias(vecEst,theta),
    Coverage=mean(as.numeric((CI$CI_low<=theta)&(theta<=CI$CI_up))),
    Power=mean(as.numeric(vec.p<=alpha))
  )
  
  return(res)
  
}


Data.loss.SWD<-function(data.all, I,K, B="0", C=0, D=0){
  
  ########  data deletion depends on:  B,C, D      ##############
  
  if(C==0){data.delete<-data.all
  }else{  
    
    #number of cluster loss
    clusternr<-sample(1:I, size=C, replace=FALSE)
    
    #timepoint of cluster loss
    switch(B,
           #B0: no Cluster missing
           "0"={timepoint<-0},
           #B1: Cluster missing at random
           "1"={timepoint<-sample(2:K, size=C, replace=TRUE)  },
           #B2: Cluster is missing at beginning
           "2"={timepoint<-sample((1:round(K/4))+1, size=C, replace=TRUE)  },
           #B3: Cluster is missing at end
           "3"={timepoint<-sample((K-round(K/4)+1):K, size=C, replace=TRUE)  }           
    )
    
    dc<-NULL
    for(i in 1:C){
      
      which.tp<-timepoint[i]:K
      dc<-rbind(dc, cbind(rep(clusternr[i],length(which.tp)), timepoint[i]:K))
    }  
    colnames(dc)<-c("cluster", "measurement")
    
    ###given vector of  timepoints and corresponding vector of cluster
    #collect rows whcih has to delete
    collect.all<-NULL
    for(i in 1:dim(dc)[1]){   
      collect<-unlist(sapply(1:dim(data.all)[1], function(j){
        
        if(identical(unname(unlist(data.all[j,c("cluster", "measurement")]),force=TRUE),unname(unlist(dc[i,]),force=TRUE))){return(j)}
      }))
      collect.all<-c(collect.all,collect )
    }
    #now delete rows in data 
    data.delete<-data.all[-collect.all, ] 
    
  }
  if(D!=0){#individiual loss
    
    loss.ind<-sample(1:dim(data.delete)[1], size=D, replace=FALSE)
    data.delete<-data.delete[-loss.ind, ] 
  }
  
  return(data.delete)
}  


simulation<-function(anzSim,type, sigma.1,sigma.2=NULL,sigma.3,K,J,I,mu.0,theta,betas,
                        X, X.A, B.cond="0",C.cond=0, D.cond=0){
  
  #DEsingmatrix Daten laut Studeinedesign
  D<-completeDataDesignMatrix(J, X)
  #Desingmatrix Daten f?r reale Daten
  if(!is.null(X.A)){A<-completeDataDesignMatrix(J, X.A)
  }else{A<-NULL}
  
  
  parameters<-c(mu.0, betas, theta)
  
  #Covarianzmatrix
  if(type=="cross-sec"){
    V<-CovMat.Design(K=K, J=J, I=I, sigma.1=sigma.1, sigma.3=sigma.3)
  }
  if(type=="long"){
    V<-CovMat.Design(K=K, J=J, I=I,sigma.1, sigma.2, sigma.3)
  }
  
  
  res.all<-NULL
  for(s in 1:anzSim){ ##Repeats simulation
    
    #sample Data
    #################
    #sample I cluster from distribution of cluster, each ave the same mean vector
    sample.data<-sampleData(type, K,J,I, D, A, V, parameters )
    #Loss of Data
    data.delete<-Data.loss.SWD(sample.data, I,K,B.cond,C.cond, D.cond)   
    #print(xtabs(~cluster+measurement, data=data.delete))
    
    ####################    Estimation by linear mixed model ########################
    #timpoint 1 <- 0
    data.delete$measurement<-as.factor((as.numeric(data.delete$measurement)-1))
    #cross-sec
    if(type=="cross-sec"){
      
      #lm.res<-lmer(val~ intervention + measurement+(1|cluster)+ (1|subject), data=sample.data)
      lm.res<-lmer(val~ intervention + measurement+(1|cluster), data=data.delete) 
      lm.0<-lmer(val~ measurement+(1|cluster), data=data.delete) 
    }
    
    #longitudinal
    if(type=="long"){
      
      lm.res<-lmer(val~ intervention + measurement+(1|cluster)+ (1|subject), data=data.delete)
      #summary(lm.res)
      lm.0<-lmer(val~ measurement+(1|cluster)+ (1|subject), data=data.delete)
      #summary(lm.res)
    }
    
    
    #random effects
    res.ran<-as.data.frame(summary(lm.res)$varcor)[,5]^2
    names(res.ran)<-as.data.frame(summary(lm.res)$varcor)[,1]
    #fixed efects
    res.fix<-fixef(lm.res)[-1]
    #SE of fixed effect estimates
    SEs<-coef(summary(lm.res))[-1,"Std. Error"]
    names(SEs)<-paste("SE.",names(SEs), sep="")
    #anova for intervention
    anova.p<-anova(lm.res,lm.0)[8][2,1]
    
    res.all<-rbind(res.all, c(res.ran,res.fix, SEs, p.intervention=anova.p))
  }
  
  #Summary of all repeats
  summEst<-colMeans(res.all)
  names(summEst)<-paste(names(summEst), "Mean",".")
  perf.intervention<-performanceMeas(res.all[,"intervention"], res.all[,"SE.intervention"], res.all[,"p.intervention"], theta)
  names(perf.intervention)<-paste(names(perf.intervention),"Intervention",sep=".")
  #c(performanceMeas(est["intervention",], est["SE.intervention",], est["p.intervention",], theta))
  return(c(summEst, perf.intervention))
}


implemMatrix.parallel<-function (nC, nT, nSw, pattern) 
{
  if (length(pattern)  == nT) {
    ma <- rbind(t(replicate(nSw, rep(0, nT))), t(replicate(nC -nSw, pattern)))
    return(ma)
  }
  else {
    stop("The length of the pattern must be the same as the number of timepoints.")
  }
}

implemMatrix.crossover<-function (nC, nT, nSw,swP=NULL, pattern) 
{
  if (is.null(swP)) {
    swP <- ceiling(nT/2)
  }
  if (length(pattern)  == swP) {
    
    ma<- rbind(t(replicate(nSw, c(rep(0, swP),  pattern))), 
               t(replicate(nC - nSw, c(pattern, rep(0, nT - swP))))
               )
    return(ma)
  }
  else {
    stop("time points of intervention must be the same as the length of the pattern")
  }
}

