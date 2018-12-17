#original here: https://github.com/cran/mixtools/blob/master/R/plot.mixEM.R
#plot.mixEM <-function (x, whichplots = 1, 
#                       loglik = 1 %in% whichplots, 
#                       density = 2 %in% whichplots,
#                       xlab1="Iteration", ylab1="Log-Likelihood",
#                       main1="Observed Data Log-Likelihood", col1=1, lwd1=2,
#                       xlab2=NULL, ylab2=NULL, main2=NULL, col2=NULL, 
#                       lwd2=2,
#                       alpha = 0.05, marginal = FALSE, ...)  {
  

plotmixemcustom <- function(mixmdl,lwdnmix,ylim2,col11,col22,nbreaks,lwdn,tmpxlab,ylab,main2,ax,cutage) {
  
  mix.object=mixmdl
  
  #add a bit to maxy
  maxyaddon=0.01
  #maxyaddon=0
  
  #if (mix.object$ft == "normalmixEM") {
  k <- ncol(mix.object$posterior)
  xm <- sort(mix.object$x)
  a <- hist(xm, plot = FALSE)
  maxy <- max(max(a$density), .3989*mix.object$lambda/mix.object$sigma)
  #ylim2=maxy+maxyaddon
  main2 <- main2 ;  xlab2 <- tmpxlab  ;
  #plot histogram
  nxticks = 3; nyticks=6
  if (ax==0) {
  hist(xm, prob = TRUE, main = main2, xlab = xlab2, 
           col=col11,border = bordercol, breaks=nbreaks,axes=FALSE)
    
    xlimlow = ceiling(par("usr")[1]/10)*10 ; xlimhigh = floor(par("usr")[2]/10)*10
    ylimlow = round(par("usr")[3]) ; tmpylimhigh = (ceiling(par("usr")[4]*100)/100) ; ylimhigh = round(tmpylimhigh-(tmpylimhigh/nyticks),2)
    xseq = seq(xlimlow,xlimhigh,length=nxticks)
    yseq = seq(ylimlow,ylimhigh,length=nyticks)
    axis(side=1,at=c(xseq),labels=FALSE)
    axis(side=2,at=c(yseq),labels=FALSE)
    
    } else {
      hist(xm, prob = TRUE, main = main2, xlab = xlab2, 
           col=col11,border = bordercol, breaks=nbreaks,axes=FALSE)
      
      xlimlow = ceiling(par("usr")[1]/10)*10 ; xlimhigh = floor(par("usr")[2]/10)*10
      ylimlow = round(par("usr")[3]) ; tmpylimhigh = (ceiling(par("usr")[4]*100)/100) ; ylimhigh = round(tmpylimhigh-(tmpylimhigh/nyticks),2)
      xseq = seq(xlimlow,xlimhigh,length=nxticks)
      yseq = seq(ylimlow,ylimhigh,length=nyticks)
      axis(side=1,at=c(xseq),labels=TRUE)
      axis(side=2,at=c(yseq),labels=TRUE)
      if (cutage > 0) {
        lines(rep(cutage,nyticks),yseq)
      #lines(rep(39.5,nyticks),yseq)
      }
    }
      
      if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {
        arbmean <- TRUE
        arbvar <- TRUE
      }
      #add mix lines
      for (i in 1:k) {
        lines(xm, mix.object$lambda[i] * 
                dnorm(xm, mean = mix.object$mu[i * arbmean + (1 - arbmean)], 
                      sd = mix.object$sigma[i * arbvar + (1 - arbvar)]), 
              col = col22[i], lwd = lwdnmix)
      }
  #add density
  lines(density(xm), lty = 3, lwd = lwdn)
}