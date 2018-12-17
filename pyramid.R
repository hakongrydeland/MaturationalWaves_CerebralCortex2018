pyramid = function(binvector,measvector,sexvector,gapn,outputdirfigs,outputname) {

#https://stats.stackexchange.com/questions/2455/how-to-make-age-pyramid-like-plot-in-r



#this is best, because you need to include age, fiq and sex
#bin 
#percentiles
x=age
#deciles
deciles=quantile(x, probs = seq(0, 1, length = 11),na.rm=TRUE)
agelabels=vector(len=length(deciles)-1)
xx.pop=xy.pop=xy.se=xx.se=xx.sd=xy.sd=xx.l=xy.l=nwhichd12=vector(len=length(agelabels))
for (di in 1:length(agelabels)) {
  whichd12=vector(len=0)
  print(paste("di is ",di,sep=""))
  
    #get age bins based on deciles
  #get which index of participants in each bin
  if (di==1) {
    tmpdi1=floor(deciles[di]) ; print(tmpdi1)
    tmpdi2=round(deciles[di+1]) ; print(tmpdi2)
    whichd12 = which(x < tmpdi2+1) ; print(whichd12)
    nwhichd12[di]=length(whichd12)
    } else if (di==length(agelabels)) {
      tmpdi1=round(deciles[di]+1) ; print(tmpdi1)
      tmpdi2=round(deciles[di+1]) ; print(tmpdi2)
      whichd12 = which(x >= tmpdi1) ; print(whichd12)
      nwhichd12[di]=length(whichd12)
     } else {
    tmpdi1=round(deciles[di])+1 ; print(tmpdi1)
    tmpdi2=round(deciles[di+1]) ; print(tmpdi2)
    whichd12 = which(x >= tmpdi1 & x < tmpdi2+1) ; print(whichd12)
    nwhichd12[di]=length(whichd12)
  }
  agelabels[di]=c(paste(tmpdi1,"-",tmpdi2,sep=""))
  
  #get fiq values based on index
  #xx.pop[di]=round(mean(fiq[which(sex[whichd12] == "Female")])) ; print(xx.pop[di])
  #xy.pop[di]=round(mean(fiq[which(sex[whichd12] == "Male")])) ; print(xy.pop[di])
  if (length(levels(measvector))==0) {
  xx.pop[di]=mean(measvector[which(sexvector[whichd12] == "Female")]) ; print(xx.pop[di]) #rounding done below.
  xy.pop[di]=mean(measvector[which(sexvector[whichd12] == "Male")]) ; print(xy.pop[di])
  #calculate se for error bars as well
  xx.sd[di]=sd(measvector[which(sexvector[whichd12] == "Female")]) ; xx.l[di] = length(measvector[which(sexvector[whichd12] == "Female")]) ;
  xy.sd[di]=sd(measvector[which(sexvector[whichd12] == "Male")]) ; xy.l[di] = length(measvector[which(sexvector[whichd12] == "Male")]) ;
  xx.se[di]=xx.sd[di]/sqrt(xx.l[di]) ; xy.se[di]= xy.sd[di]/sqrt(xy.l[di]) ;
  } else {
    xx.pop[di]=length(measvector[which(sexvector[whichd12] == "Female")]) ; print(xx.pop[di]) #rounding done below.
    xy.pop[di]=length(measvector[which(sexvector[whichd12] == "Male")]) ; print(xy.pop[di])
    #calculate se for error bars as well
    xx.se[di]=0 ; xy.se[di]= 0
  }
}

#http://colorbrewer2.org/#type=sequential&scheme=RdPu&n=6

library(plotrix)
source("/Users/hakongrydeland/Dropbox/cambridge2016/projects/lifespant1wt2w/splineclusteringRcamcode/pyramid.plot.custom.R")
source("/Users/hakongrydeland/Dropbox/cambridge2016/projects/lifespant1wt2w/splineclusteringRcamcode/pyramid.plot.custom.pub.R")


fcolfunc <- colorRampPalette(c(rgb(201,148,199,max=255), rgb(212,185,218,max=255)))
mcolfunc <- colorRampPalette(c(rgb(152,0,67,max=255), rgb(221,28,119,max=255)))
#gapn=25
pdf(paste(outputdirfigs,outputname,'.pdf',sep=''),width=6,height=5)
#par(mar=pyramid.plot(xx.pop,xy.pop,labels=agelabels,
 #                    main="FIQ divided by age and sex",lxcol=mcolfunc(length(xy.pop)),rxcol=fcolfunc(length(xx.pop)),
#                     gap=25,show.values=TRUE,unit="Mean FIQ",ppmar=c(4,4,3.5,3.5),ndig=0))
par(mar=pyramid.plot.custom(xy.pop,xx.pop,xx.se,xy.se,labels=agelabels,lxcol=mcolfunc(length(xy.pop)),rxcol=fcolfunc(length(xx.pop)),
                     gap=gapn,show.values=TRUE,unit="N",ppmar=c(4,4,3.5,3.5),ndig=0))


dev.off()

}
