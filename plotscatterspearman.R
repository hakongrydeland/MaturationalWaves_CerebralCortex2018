plotscatterspearman <- function(tmpx,tmpy,fignum,xname,yname,xlabvar,ylabvar) {
  print(ylabvar)
  library(outliers)
  library(Hmisc)
  #outlierthresh=0.99
  potoulx=which(outlier(tmpx,logical = TRUE)) # potential outier x
  potouly=which(outlier(tmpy,logical = TRUE)) # potential outier y
  #which(scores(tmpx, type="z", prob=outlierthresh))
  #which(scores(tmpy, type="z", prob=outlierthresh))
  #tmpx[potoulx]=NA ; tmpy[potouly]=NA ;
  x=tmpx[!is.na(tmpx)] ;y=tmpy[!is.na(tmpy)]
  #tmpcor = cor(x,y,type="pearson")  ; print(tmpcor)
  tmpcor = rcorr(x,y,type=c("spearman"))  ; print(tmpcor)
  
  rhoandpval=paste("Spearman's r = ",tmpcor$r[1,2],", p = ",tmpcor$P[1,2],sep="")
  write(rhoandpval, file = paste(outputdirsummeas,'rhoandpval_',xname,'_',yname,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
  
  
  pdf(paste(outputdirfigs,'fig',fignum,'_scatterspearman_',xname,'_',yname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  par(mar=c(4,4,2,0.5),bty="l")
  plot(0, type='n',xlab=xlabvar, ylab=ylabvar,xlim=c(min(x),max(x)), ylim=c(min(y),max(y)))#ylim=c(min(Y),max(Y)))
  #plot(0, type='n',xlab=xlabvar, ylab=ylabvar,xlim=c(min(cbind(x,y)),max(cbind(x,y))), ylim=c(min(cbind(x,y)),max(cbind(x,y))))#ylim=c(min(Y),max(Y)))
  points(x,y, col = 'black', pch = 16, cex = 0.75) # plot data#title('maturational trajectory matrix', line=+2, font.main=1)
  abline(lm(y~x), col="gray") # regression line (y~x) 
  axislims=par("usr")
  if (tmpcor$r[1,2] >= 0) {
    #legend('topleft', legend=sprintf("r = %3.2f", tmpcor), bty='n', cex=0.7,text.col = pointcol1)
    legend(min(x),max(y),legend=sprintf("r = %3.2f", round(tmpcor$r[1,2], digits=2)), bty='n', cex=0.9,text.col = pointcol1)
  } else {
    #legend('topright', legend=sprintf("r = %3.2f", tmpcor), bty='n', cex=0.7,text.col = pointcol1)
    legend(max(x)-max(x)*0.3,max(y),legend=sprintf("r = %3.2f", round(tmpcor$r[1,2], digits=2)), bty='n', cex=0.9,text.col = pointcol1)
  }
  dev.off()
}