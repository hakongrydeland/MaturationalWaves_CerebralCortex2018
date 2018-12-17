plotvioregression <- function(tmpx,tmpy,fignum,xname,yname,xlabvar,ylabvar,col1,col2,col3) {
  print(ylabvar)
  library(outliers)
  library(Hmisc)
  #regression
  x=tmpx[!is.na(tmpx)] ;y=tmpy[!is.na(tmpy)]
  xfactor=factor(x)
  tmpcatreg = lm(as.numeric(y) ~ xfactor)
  sumtmpcatreg = summary(tmpcatreg)
  fstats = sumtmpcatreg$fstatistic
  rsquared = sumtmpcatreg$r.squared
  varofinterescovs = matrix(nrow = length(levels(xfactor))-1, ncol=4) #4 cols are: Estimate - Std. Error - t value -  Pr(>|t|) 
  for (vi in 1:(length(levels(xfactor))-1)) {
    varofinterescovs[,]=sumtmpcatreg$coefficients[2,]
  }
  
  ## funciton to get p value
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
  modelpvalue = lmp(tmpcatreg)
  
  frandp=paste("F(",fstats[2],",",fstats[3],") = ",fstats[1],", and R2 = ",rsquared," (P = ",modelpvalue,")",sep="")
  write(frandp, file = paste(outputdirsummeas,'frandp_',xname,'_',yname,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
  #write(varofinterescovs, file = paste(outputdirsummeas,'varofinterescovs_',xname,'_',yname,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
  
  pdf(paste(outputdirfigs,'fig',fignum,'_vioplot_',xname,'_',yname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  par(mar=c(4,4,2,0.5),bty="l")
  #plot(0, type='n',xlab=xlabvar, ylab=ylabvar,xlim=c(1,length(levels(xfactor))), ylim=c(min(y),max(y)))#ylim=c(min(Y),max(Y)))
  #plot(0, type='n',xlab=xlabvar, ylab=ylabvar,xlim=c(min(cbind(x,y)),max(cbind(x,y))), ylim=c(min(cbind(x,y)),max(cbind(x,y))))#ylim=c(min(Y),max(Y)))
  
  viocommand=""
  for (vi in 1:(length(levels(xfactor)))) {
    if (vi == 1) {
      viocommand=paste(viocommand,"as.numeric(tmpy[xfactor==levels(xfactor)[",vi,"]])",sep="")
    } else
      viocommand=paste(viocommand,",as.numeric(tmpy[xfactor==levels(xfactor)[",vi,"]])",sep="")
  }
  
  print(viocommand)
  par(las=1,bty="l")  ## my preferred setting
  ## set up empty plot
  tmpylim =c(min(tmpy),max(tmpy))
  plot(0, type='n', xlab=xlabvar, ylab=ylabvar, xlim=c(0,length(levels(tmpx))+1), ylim=tmpylim,ann=FALSE,xaxt="n")
  
  #vioplot(as.numeric(tmpy[xfactor==levels(xfactor)[1]]),as.numeric(tmpy[xfactor==levels(xfactor)[2]]),col=c(col2))
  if (length(levels(xfactor)) == 3) {
    eval(parse(text = paste("vioplot(",viocommand,",col=c(col3),add=TRUE)")))
  } else {
    eval(parse(text = paste("vioplot(",viocommand,",col=c(col2),add=TRUE)")))
  }
  
  for (vi in 1:(length(levels(xfactor))-1)) {
    print(vi)
    if (vi == 1) {
      vioplot(as.numeric(tmpy[xfactor==levels(xfactor)[vi]]),col=col1,add=TRUE,at=vi)
      } else {
        vioplot(as.numeric(tmpy[xfactor==levels(xfactor)[vi]]),col=col2,add=TRUE,at=vi)
      }
  }
  title(ylab=ylabvar,xlab=xlabvar)
  axis(side=1,at=c(1,2),labels=c("Wave 1","Wave 2"))
  dev.off()
  
  ### publication figure ###
  
  pdf(paste(outputdirfigs,'fig',fignum,'_pub_vioplot_',xname,'_',yname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  par(mar=c(0.6,0.6,0.5,0.5),bty="l")  #w/o labels
  #plot(0, type='n',xlab=xlabvar, ylab=ylabvar,xlim=c(1,length(levels(xfactor))), ylim=c(min(y),max(y)))#ylim=c(min(Y),max(Y)))
  #plot(0, type='n',xlab=xlabvar, ylab=ylabvar,xlim=c(min(cbind(x,y)),max(cbind(x,y))), ylim=c(min(cbind(x,y)),max(cbind(x,y))))#ylim=c(min(Y),max(Y)))
  
  viocommand=""
  for (vi in 1:(length(levels(xfactor)))) {
    if (vi == 1) {
      viocommand=paste(viocommand,"as.numeric(tmpy[xfactor==levels(xfactor)[",vi,"]])",sep="")
    } else
      viocommand=paste(viocommand,",as.numeric(tmpy[xfactor==levels(xfactor)[",vi,"]])",sep="")
  }
  
  print(viocommand)
  par(las=1,bty="l")  ## my preferred setting
  ## set up empty plot
  tmpylim =c(min(tmpy),max(tmpy))
  plot(0, type='n', xlim=c(0.6,length(levels(tmpx))+0.4), ylim=tmpylim,ann=FALSE,xlab="",ylab="",yaxt="n",xaxt="n")
  
  #vioplot(as.numeric(tmpy[xfactor==levels(xfactor)[1]]),as.numeric(tmpy[xfactor==levels(xfactor)[2]]),col=c(col2))
  if (length(levels(xfactor)) == 3) {
    eval(parse(text = paste("vioplot(",viocommand,",col=c(col3),add=TRUE)")))
  } else {
    eval(parse(text = paste("vioplot(",viocommand,",col=c(col2),add=TRUE)")))
  }
  
  for (vi in 1:(length(levels(xfactor))-1)) {
    print(vi)
    if (vi == 1) {
      vioplot(as.numeric(tmpy[xfactor==levels(xfactor)[vi]]),col=col1,add=TRUE,at=vi)
    } else {
      vioplot(as.numeric(tmpy[xfactor==levels(xfactor)[vi]]),col=col2,add=TRUE,at=vi)
    }
  }
  #add y axis w/o labels
  axis(2,labels=FALSE) 
  #add x axis w/o labels
  axis(side=1,at=c(1,2),labels=FALSE)
  dev.off()
  
}
