plotsplitvioplot <- function(tmpcat1,tmpcat2,tmpx,fignum,xname,yname,xlabvar,ylabvar,col1,col2) {
  
  #get libararies
  require(vioplot)
  require(devtools)
  require(digest)
  source("/Users/hakongrydeland/Dropbox/cambridge2016/projects/lifespant1wt2w/splineclusteringRcamcode/vioplot2.R")
  
  #set output file
  pdf(paste(outputdirfigs,'fig',fignum,'_splitvioplot_',xname,'_',yname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  
  #par(mar=c(4,4,2,0.5),bty="l")
  
  plot(x=NULL, y=NULL,xlim = c(0.5, 2.5), ylim=c(min(tmpx), max(tmpx)),type="n", ann=FALSE, axes=F)
  axis(1, at=c(1, 2),  labels=c(levels(factor(tmpcat1))[1], levels(factor(tmpcat1))[2]))
  axis(2)
  for (i in unique(tmpcat1)) {
      for (j in unique(tmpcat2)){
          vioplot2(tmpx[which(tmpcat1 == i & tmpcat2 == j)],
          at = ifelse(i == levels(factor(tmpcat1))[1], 1,2),
          side = ifelse(j == levels(factor(tmpcat2))[1], "left", "right"),
          col = ifelse(j == levels(factor(tmpcat2))[1], col1, col2),
          add = T)
      }
  }
  title("", xlab=xlabvar,ylab=ylabvar)
  #legend("bottomright", fill = c(col1, col2),
  #legend = c(levels(factor(tmpcat2))[1], levels(factor(tmpcat2))[2]), box.lty=0)
  
  dev.off() 
  
  #publication file without labels
  #set output file
  pdf(paste(outputdirfigs,'fig',fignum,'_pub_splitvioplot_',xname,'_',yname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  par(mar=c(0.8,0.8,0.5,0.5),bty="l")
  plot(x=NULL, y=NULL,xlim = c(0.6, 2.5), ylim=c(min(tmpx), max(tmpx)),type="n", ann=FALSE, axes=F)
  axis(1, at=c(1, 2),  labels=FALSE, levels(factor(tmpcat1))[2])
  axis(2,labels=FALSE)
  for (i in unique(tmpcat1)) {
    for (j in unique(tmpcat2)){
      vioplot2(tmpx[which(tmpcat1 == i & tmpcat2 == j)],
               at = ifelse(i == levels(factor(tmpcat1))[1], 1,2),
               side = ifelse(j == levels(factor(tmpcat2))[1], "left", "right"),
               col = ifelse(j == levels(factor(tmpcat2))[1], col1, col2),
               add = T)
    }
  }
  dev.off() 
 
}
