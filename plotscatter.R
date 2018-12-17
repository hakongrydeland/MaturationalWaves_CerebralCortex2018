plotscatter <- function(tmpx,tmpy,fignum,xname,yname,xlabvar,ylabvar) {
  
  rp.main = function(rsq,p,nsig) {
    main = bquote(r^2 ~ '=' ~ .(toString(signif(rsq,nsig))) ~ ', p =' ~ .(toString(signif(p,nsig)))) 
    return(main)
  }
  
  print(ylabvar)
  library(outliers)
  #outlierthresh=0.99
  potoulx=which(outlier(tmpx,logical = TRUE)) # potential outier x
  potouly=which(outlier(tmpy,logical = TRUE)) # potential outier y
  #which(scores(tmpx, type="z", prob=outlierthresh))
  #which(scores(tmpy, type="z", prob=outlierthresh))
  #tmpx[potoulx]=NA ; tmpy[potouly]=NA ;
  x=tmpx[!is.na(tmpx)] ;y=tmpy[!is.na(tmpy)]
  tmpcor = cor(x,y) ; print(tmpcor)
  #fit linear model
  l = lm(tmpy~x); 
  f = summary(l)$fstatistic; p = pf(f[1],f[2],f[3], lower.tail = FALSE);rsq = summary(l)$r.squared ;
  pdf(paste(outputdirfigs,'fig',fignum,'_scatter_',xname,'_',yname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  par(mar=c(4,4,2,0.5),bty="l")
  plot(0, type='n',xlab=xlabvar, ylab=ylabvar,xlim=c(min(x),max(x)), ylim=c(min(y),max(y)),main = rp.main(rsq,p,2),)#ylim=c(min(Y),max(Y)))
  #plot(0, type='n',xlab=xlabvar, ylab=ylabvar,xlim=c(min(cbind(x,y)),max(cbind(x,y))), ylim=c(min(cbind(x,y)),max(cbind(x,y))))#ylim=c(min(Y),max(Y)))
  points(x,y, col = 'black', pch = 16, cex = 0.75) # plot data#title('maturational trajectory matrix', line=+2, font.main=1)
  abline(lm(y~x), col="gray") # regression line (y~x) 
  axislims=par("usr")
  if (tmpcor >= 0) {
    #legend('topleft', legend=sprintf("r = %3.2f", tmpcor), bty='n', cex=0.7,text.col = pointcol1)
    legend(min(x),max(y),legend=sprintf("r = %3.2f", round(tmpcor, digits=2)), bty='n', cex=0.9,text.col = pointcol1)
  } else {
    #legend('topright', legend=sprintf("r = %3.2f", tmpcor), bty='n', cex=0.7,text.col = pointcol1)
    legend(max(x)-max(x)*0.3,max(y),legend=sprintf("r = %3.2f", round(tmpcor, digits=2)), bty='n', cex=0.9,text.col = pointcol1)
  }
  dev.off()
  
  #plot publication figure
  
  #get ci
  #http://rstudio-pubs-static.s3.amazonaws.com/24365_2803ab8299934e888a60e7b16113f619.html for polygon and bootstrap
  #boot <- replicate(999, mean(sample(Sample, replace=T)))
  
  #polygon and predict from frantisek (see sc.window.hakon_orig.R)
  library(seqinr) # col2alpha
  x = tmpx;
  x.pred = seq(from=min(tmpx),to=max(tmpx),by=0.0001) 
  pred = predict(l,newdata=data.frame('x'=x.pred),interval='confidence')
  #plot(tmpx,tmpy,xaxt='n',xlab=expression(paste(Delta,' CT (mm/year)',sep='')),ylab=expression(paste(Delta, ' strength')), pch = 19, col='black', xlim=c(-0.04,0.015), ylim=c(-0.04,0.03)); 
  #axis(1,at=seq(-0.04,0.01,b=0.01),labels=c(-0.04,'',-0.02,'',0,''))
  pdf(paste(outputdirfigs,'fig',fignum,'_pub_scatter_',xname,'_',yname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  par(mar=c(0.8, 0.8, 0.5,0.5) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white',bty="l")
  plot(0, type='n', xlim=c(min(tmpx),max(tmpx)), ylim=c(min(tmpy),max(tmpy)),xlab="",ylab="",yaxt="n",xaxt="n") ; axis(1,labels=FALSE) ;  axis(2,labels=FALSE) ;  #w/o labels (for figures)
  points(tmpx,tmpy, col = 'gray', pch = 16, cex = 0.75) # plot data#title('maturational trajectory matrix', line=+2, font.main=1)
  selwd=1.5 ; selty=1
  polygon(c(rev(x.pred), x.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
  lines(x.pred,pred[,1],col='black',lwd=3); 
  dev.off()
}