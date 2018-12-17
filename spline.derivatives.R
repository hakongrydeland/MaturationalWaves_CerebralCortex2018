rm(list = ls()) # clear all


library(R.matlab)
library(vows)
library(cluster)
library(fda)
library(clue)
library(pracma) #only needed for trimming of white space of labelnames
library(Hmisc)
library(rcompanion) ; library(corrplot) ; library(RColorBrewer) ; library(Rmisc) ; library(reshape2) 
library(Rarity) ;library(mixtools) ; library(mclust) ; library(extRemes) ; library(diptest) ; library(fpc) ; library(vioplot) ; library(gmodels)
library(lm.beta)

source("plotscatter.R")
source("plotscatterspearman.R")
source("plotsplitvioplot.R")
source("plotvioregression.R")
source("pyramid.R")
source("pyramid.pub.R")
source("plotmixemcustom.R")

#set/create input file and output directory
basedir="projects/lifespant1wt2w/"
subdatadir="results/"
measname=c("myelin_fraction70")
moi="T1w/T2w (a.u.)"
print(measname)
measdir=paste(basedir,subdatadir,measname,'/',sep='')
ifile = paste(measdir,measname,".mat",sep='') #generated in matlab, using myelinmapping.m

outputdir=paste(measdir,"spline.derivatives/",sep='')
outputdirsummeas=paste(outputdir,"summarymeasures/",sep='')
outputdirfigs=paste(outputdir,"figures/",sep='')
#create outputdir
dir.create(file.path(outputdir))
dir.create(file.path(outputdirsummeas))
dir.create(file.path(outputdirfigs))
#set directory to get SC results
struccovdir=paste(measdir,"structuralcovariance_btwm2m3/",sep="")

# load data and set up variables
mri = readMat(paste(ifile,sep=''));  #attach(mri)
origY = mri$parcmeas # dimension n.subjects * n.regions

#demog
sex = factor(mri$sex, levels = c(0,1), labels = c('Male','Female')) #female=1
age = as.vector(mri$age, mode = 'numeric')

#cog
wasimat = mri$wasi
blockdesign=wasimat[,1] #nan's
matrixreasoning=wasimat[,2] #nan's
vocabulary=wasimat[,3] #nan's
similarities=wasimat[,4] #nan's
viq=wasimat[,5] #nan's
piq=wasimat[,6] #nan's
fiq=wasimat[,7] #nan's

#get "population plot" for fiq seperately for females and males, binned by age
pyramidoutputname=paste("pyramidplot.fiq.age.sex.",measname,sep="")
gapn=25
pyramid(age,fiq,sex,gapn,outputdirfigs,pyramidoutputname)

#get "population plot" for gender, binned by age
pyramidoutputname=paste("pyramidplot.sex.age.sex.",measname,sep="")
gapn=5
pyramid(age,sex,sex,gapn,outputdirfigs,pyramidoutputname)

### motion - start ### 
# (HG: the analyses based on this, in the SI, are just re-running everything below after removing the participants with values below the 10%ile)
#tsnr from dwi (Roalf et al. 2016 NI)
dwimotion = mri$dwimotion # 9 columns - 1: mean of RMS to 1 volume; 2: mean of RMS to previous volume; 3: discards trans and rot in PE dir - mean of RMS to 1 volume; 4: discards trans and rot in PE dir - mean of RMS to previous volume;
#last four column is std instead of mean from the first 4 columns ; last column is tsnr (Roalf et al. 2016 NI)
t2dim = mri$t2dim # 
dwimotionofinterest = dwimotion[,9] #Roalf et al. 2016 NI argues that tsnr is better than rms,
#replace nan's with median
dwimotionofinterestnona = dwimotionofinterest
dwimotionofinterestnona[which(is.na(dwimotionofinterest))] = median(dwimotionofinterest[which(!is.na(dwimotionofinterest))])
dwimotionofinterest = dwimotionofinterestnona
plot(age,dwimotionofinterest)
#remove effect of dim (as in the T2w, see paper for details)
t2dimeffectdwimotion = lm(dwimotionofinterest~as.factor(t2dim)) ; summary(t2dimeffectdwimotion) ; t2dimestimate = t2dimeffectdwimotion$coefficients[[2]]
t2dim01 = t2dim ; t2dim01[which(t2dim==204)] = 1 ; t2dim01[which(t2dim==256)] = 0 ;
dwimotionofinterestt2dimcorrected = dwimotionofinterest+(t2dimestimate*t2dim01)
dwimotionofinterest=dwimotionofinterestt2dimcorrected
# plot
pdf(paste(outputdirfigs,'supfig13b.age_tsnr.pdf',sep=''),width=6,height=5)
par(mar=c(4,4,0.5,4),bty="l")
plot(age,dwimotionofinterest)
dev.off()
pdf(paste(outputdirfigs,'supfig13b.age_tsnr_pub.pdf',sep=''),width=6,height=5)
par(mar=c(0.6,0.6,0.5,0.5),bty="l")  #w/o labels
plot(age,dwimotionofinterest)
dev.off()

#threhold at a certain percentile
ptile=0.10
qthre = quantile(dwimotionofinterest,ptile) ; qthreind = which(dwimotionofinterest > qthre) ; notqthreind = which(!dwimotionofinterest > qthre) ; tsnrnotqthreind = notqthreind ;
ind=qthreind # to be used to index which participants to be used/excluded

#plot density
tmpmeas=age[notqthreind] ; tmpclust = densityMclust(tmpmeas)
pdf(paste(outputdirfigs,'supfig13a.density_excluded_tsnr.pdf',sep=''),width=6,height=5)
brr = seq(min(tmpmeas),max(tmpmeas),length=50)
plot(tmpclust, what = "density",xlab="Age", data = tmpmeas,breaks=brr,frame.plot=FALSE)
dev.off()
pdf(paste(outputdirfigs,'supfig13a.density_excluded_tsnr_pub.pdf',sep=''),width=6,height=5)
par(mar=c(0.6,0.6,0.5,0.5),bty="l")  #w/o labels
plot(tmpclust, what = "density",xlab="Age", data = tmpmeas,breaks=brr,frame.plot=FALSE)
dev.off()

#Euler number from t1w (#euler number calculated on ?h.orig.nofix in FreeSurfer) Rosen et al. NI, 2017 ("Quantitative assessment of structural image quality")
eulermotion = mri$t1weulermotion #3 columns, first is mean acros the hemispheres, which are in the next to columns.
eulermotionofinterest = eulermotion[,1] 
#plot(age,eulermotionofinterest) #shows one huge "outlier" below -600
eind=which(eulermotionofinterest > -600) ; eindout = which(eulermotionofinterest < -600) ; print(eulermotionofinterest[eindout]) ;  print(age[eindout])
# plot
pdf(paste(outputdirfigs,'supfig13d.age_euler_oneoutlier.pdf',sep=''),width=6,height=5)
par(mar=c(4,4,0.5,4),bty="l")
plot(age[eind],eulermotionofinterest[eind])
dev.off()
pdf(paste(outputdirfigs,'supfig13d.age_euler_oneoutlier_pub.pdf',sep=''),width=6,height=5)
par(mar=c(0.6,0.6,0.5,0.5),bty="l")  #w/o labels
plot(age[eind],eulermotionofinterest[eind])
dev.off()


#threhold at a certain percentile
ptile=0.10
qthre = quantile(eulermotionofinterest,ptile) ; qthreind = which(eulermotionofinterest > qthre) ; notqthreind = which(!eulermotionofinterest > qthre)
plot(age[qthreind],eulermotionofinterest[qthreind])
ind=qthre # to be used to index which participants to be used/excluded

#plot density
tmpmeas=age[notqthreind] ; tmpclust = densityMclust(tmpmeas)
pdf(paste(outputdirfigs,'supfig13c.density_excluded_euler.pdf',sep=''),width=6,height=5)
brr = seq(min(tmpmeas),max(tmpmeas),length=50)
plot(tmpclust, what = "density",xlab="Age", data = tmpmeas,breaks=brr,frame.plot=FALSE)
dev.off()
pdf(paste(outputdirfigs,'supfig13c.density_excluded_euler_pub.pdf',sep=''),width=6,height=5)
par(mar=c(0.6,0.6,0.5,0.5),bty="l")  #w/o labels
plot(tmpclust, what = "density",xlab="Age", data = tmpmeas,breaks=brr,frame.plot=FALSE)
dev.off()

#compare with dti tsnr
x=notqthreind ; y=tsnrnotqthreind
intersect(x, y)
print(paste(round((length(intersect(x,y))/length(x))*100,0),"% overlap btw dti tsnr and euler number",sep=""))

### motion - done ###

#load thickness as well
measnameCT="corrthickness"
measdirCT=paste("projects/lifespant1wt2w/",subdatadir,measnameCT,'/',sep='')
ifileCT = paste(measdirCT,measnameCT,".mat",sep='') #generated in matlab
mriCT = readMat(paste(ifileCT,sep=''));  #attach(mri)
CT = mriCT$parcmeas # dimension participants * n.regions
meanCTpart = rowMeans(CT) #mean across regions, for each participants
meanCTroi = colMeans(CT) #mean across participants


#mean across participants above tage
tage=60
#tage = gageatsecondzerocrossing ;
tagebel=15
meanCTroiabtage = colMeans(CT[which(age >= tage),]) 
meanCTroibeltage = colMeans(CT[which(age <= tage),]) 

#ROI info
hcplabelname = strTrim(mri$roiname)
labelname = as.matrix(c(hcplabelname))
hcprgbrois = mri$rgbrois
rgbrois = as.matrix(rbind(hcprgbrois))
roisize = mri$sizerois

nhcprois=length(hcplabelname)
Y = origY
dim(Y)

nsubs = dim(Y)[1] # number of subjects
nroi = dim(Y)[2] # number of regions

# Fit splines ----------------------------------------------------------
k = 8 # number of splines (3.5 gives error, 4 as well for _frac30)
bs0 = "bq" # spline bases - knots placed at quantiles of the data
norder0 = 4 # norderd - gives cubic B-splines #3
pen.order0 = 2 # pen.orderd

# fit splines global first
#set alpha level for CI for derivatives
galpha=10^-4
alpha=galpha

x=age
y=rowMeans(Y)
gmodel = semipar.mix.mp(Y = as.matrix(y), x = age, param = NULL, random = NULL, k = k, data.ran = NULL, norder = norder0, pen.order = pen.order0, knots = "quantile", store.gamm4 = TRUE)
summary(gmodel)
sum(gmodel$edf) #check df
# plot
# prepare splines for plotting using predict.gam
gtmpmod = gmodel$gamm4.list[[1]]$gam #tmp model
framelength=1000
gfake.frame <- data.frame(x=seq(min(age),max(age),length=framelength) ) # might included "Female" as above
X=gfake.frame
gdata.fit = predict(gtmpmod, newdata = X, se.fit=T)
newx=gfake.frame$x
gcurve = gdata.fit$fit
gcurvese = gdata.fit$se

# derivative (based on Deriv function http://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/)
gdata.fit.lpmat <- predict(gtmpmod,newdata = X,type="lpmatrix") #same as above, but now saved as matrix with "lpmatrix"
gdYm = diff(gdata.fit.lpmat)/diff(newx)  # the derivative of the lpmatrix
gdY = gdYm %*% coef(gtmpmod) # multiply with the coefficients to get one vector
gdX <- rowMeans(embed(newx,2)) # centers the X values for plotting
# se.derivative
gtmpmodvp = gtmpmod$Vp
gdY.sd <- rowSums(gdYm %*% gtmpmodvp * gdYm)^.5
gresidual.df <- df.residual(gtmpmod)
gtVal <- qt(1 - (galpha/2), gresidual.df)
gupr <- gdY + gtVal * gdY.sd
glwr <- gdY - gtVal * gdY.sd

# get different properties of the global curve
# max rate of change (1)
gmaxderiv = max(gdY) # save for brain plotting and histograms

# age at max rate of change (2)
gmaxderivind = which.max(gdY)
gageatmaxderiv = newx[gmaxderivind] # save for brain plotting and histograms
gfuncvalatmaxderiv = gcurve[gmaxderivind]

# age at 1st zero crossing (age at where lower deriv.se crosses zero first) - only works if inverted-U shaped?
# first, find where changes from pos to neg (given that it's an inverted-U shape)
#gmaxind = which.max(gcurve) #this does not work for 30% where the max is at the end, and there is a zero crossing at the end of the plateau for lower
gmaxind = gmaxderivind
# lwr
gnewind=as.integer(gmaxind)
setderivmeasofi=glwr # deriv measure of interest
# find were the derivmeasofi goes from pos to neg (i.e. hits/crosses zero)
tmpval=setderivmeasofi[gmaxind] #this only works if curve is inverted U-shaped, and not if U-shaped (then you would want the upperderivse)
#while (isTRUE((tmpval>0)==(setderivmeasofi[gmaxind]>0))) {gnewind=gnewind-1 ; tmpval=setderivmeasofi[gnewind]; if (!isTRUE((tmpval>0)==(setderivmeasofi[gmaxind]>0))) break; }
while (isTRUE((tmpval>0)==(setderivmeasofi[gmaxind]>0))) {gnewind=gnewind+1 ; tmpval=setderivmeasofi[gnewind]; if (!isTRUE((tmpval>0)==(setderivmeasofi[gmaxind]>0))) break; }
print(gnewind)
gageatlastpos = newx[gnewind]
gageatlastpospluss1 = newx[gnewind+1]
gageatfirstzerocrossing=mean(c(gageatlastpos,gageatlastpospluss1)) # save for brain plotting and histograms
gderivvalatageatfirstzerocrossing=mean(c(setderivmeasofi[gnewind],setderivmeasofi[gnewind+1]))
gvalatageatfirstzerocrossing=mean(c(gcurve[gnewind],gcurve[gnewind+1]))

# age at 2nd zero crossing (age at where upper deriv.se crosses zero first) - only works if inverted-U shaped?
# upr
if (gnewind==length(newx)) {
  print("rising curve - no second crossing?") 
  gageatsecondzerocrossing=newx[gnewind]
} else {
  setderivmeasofi=gupr # deriv measure of interest
  # find were the derivmeasofi goes from pos to neg (i.e. hits/crosses zero)
  tmpval=setderivmeasofi[gnewind] #this only works if curve is inverted U-shaped, and not if U-shaped (then you would want the upperderivse)
  while (isTRUE((tmpval>0)==(setderivmeasofi[gmaxind]>0))) {gnewind=gnewind+1 ; tmpval=setderivmeasofi[gnewind]; if (!isTRUE((tmpval>0)==(setderivmeasofi[gmaxind]>0))) break; }
  print(gnewind)
  if (gnewind==length(newx)) {
    print("rising curve - no second crossing?") 
    gageatsecondzerocrossing=newx[gnewind]
  } else {  
    gageatlastpos2 = newx[gnewind]
    gageatlastpos2pluss1 = newx[gnewind+1]
    gageatsecondzerocrossing=mean(c(gageatlastpos2,gageatlastpos2pluss1)) # save for brain plotting and histograms
    gderivvalatageatsecondzerocrossing=mean(c(setderivmeasofi[gnewind],setderivmeasofi[gnewind+1]))
    gvalatageatsecondzerocrossing=mean(c(gcurve[gnewind],gcurve[gnewind+1]))
  }
}

#  
#deriv 50% after 2nd zero crossing and age at this value (decline, if decline) after 2nd zero crossing 
tmpstart1=gnewind
tmpend1=length(gdY)
g50percfromsecondcrossingind=tmpstart1+((tmpend1-tmpstart1)/2)
g50pdecline=gdY[g50percfromsecondcrossingind] # save for brain plotting and histograms
gageat50pdecline = newx[g50percfromsecondcrossingind] # save for brain plotting and histograms
gvalatmaxdecline = gcurve[g50percfromsecondcrossingind]

# plot spline - fig 1B
# set colors #from: http://colorbrewer2.org/#type=sequential&scheme=PuRd&n=5 (for colors 5-7)
col1 = matrix(c(208,28,139),nrow=1,ncol=3) #spline
col2 = matrix(c(77,172,38),nrow=1,ncol=3) #deriv
col3 = matrix(c(241,182,218),nrow=1,ncol=3) #se spline lines
col4 = matrix(c(184,225,134),nrow=1,ncol=3) #se deriv lines
col5 = matrix(c(215,181,216),nrow=1,ncol=3) #max rate of change/increase (or myelination)
col6 = matrix(c(223,101,176),nrow=1,ncol=3) #1st zero crossing
col7 = matrix(c(221,28,119),nrow=1,ncol=3) #2nd zero crossing
col8 = matrix(c(152,0,67),nrow=1,ncol=3) #max rate of decline (if any)
colmat = rbind(col1,col2,col3,col4,col5,col6,col7,col8)
# set output name
pdf(paste(outputdirfigs,'fig1b_globalmean',measname,'.pdf',sep=''),width=6,height=5)
par(mar=c(4,4,0.5,0.5),bty="l") #w/ labels
tmpscatterw = 12
tmpscatterh = 10
xlabage='Age (years)'
ylabmyel='Myelination (T1w/T2w)'
plot(0, type='n', xlab=xlabage, ylab=ylabmyel, xlim=c(min(x),max(x)), ylim=c(min(y),max(y))) #w/ labels
points(x,y, col = 'black', pch = 16, cex = 0.75) # plot data#title('maturational trajectory matrix', line=+2, font.main=1) 
gtmpmodcol=rgb(colmat[1,1], colmat[1,2], colmat[1,3], max = 255) ; col1 = gtmpmodcol
gtmpmodcolse=rgb(colmat[3,1], colmat[3,2], colmat[3,3], maxColorValue=255)
selwd=1.5 ; selty=1
lines(newx,gcurve, type = 'l', lty = 1, lwd = 3, col = gtmpmodcol)
lines(newx,gcurve-2*gcurvese, type = 'l', lty = selty, lwd = selwd, col = gtmpmodcolse)
lines(newx,gcurve+2*gcurvese, type = 'l', lty = selty, lwd = selwd, col = gtmpmodcolse)
dev.off()

#plot spline + deriv in one plot - Fig. 1C
pdf(paste(outputdirfigs,'fig1c_globalmean_deriv',measname,'.pdf',sep=''),width=6,height=5)
par(mar=c(4,4,0.5,4),bty="u")
plot(0, type='n',xlab=xlabage, ylab="",xlim=c(min(x),max(x)), ylim=c(min(y),max(y)))#ylim=c(min(Y),max(Y)))
mtext(ylabmyel,side=2,line=3,col=gtmpmodcol)
lines(newx,gcurve, type = 'l', lty = 1, lwd = 3, col = gtmpmodcol)
lines(newx,gcurve-2*gcurvese, type = 'l', lty = selty, lwd = selwd, col = gtmpmodcolse)
lines(newx,gcurve+2*gcurvese, type = 'l', lty = selty, lwd = selwd, col = gtmpmodcolse)
# add derivative
#gtmpmodderivcol="red"
gtmpmodderivcol=rgb(colmat[2,1], colmat[2,2], colmat[2,3], max = 255) ; col2=gtmpmodderivcol
gtmpmodderivcolse=rgb(colmat[4,1], colmat[4,2], colmat[4,3], maxColorValue=255)
par(new=TRUE)
plot(gdX,gdY,type="l",lwd = 3,lty = 1,col=gtmpmodderivcol,xaxt="n",yaxt="n",xlab="",ylab="",axes=F,ylim=c(min(glwr),max(gupr))) # not wide enough for CI - check
lines(gdX,gupr,lty=selty,lwd=selwd,col=gtmpmodderivcolse) #check
lines(gdX,glwr,lty=selty,lwd=selwd,col=gtmpmodderivcolse) #check
abline(h=0,lty = 3,lwd = 2, col="black") # add line at zero
axis(side=4)
y2labmyel="Myelinogenesis (d(T1w/T2w)/dt)"
mtext(y2labmyel,side=4,line=3,col=gtmpmodderivcol)
dev.off()
#publication figure 1c
pdf(paste(outputdirfigs,'fig1c_pub_globalmean_deriv_pub',measname,'.pdf',sep=''),width=6,height=5)
par(mar=c(0.6,0.6,0.5,0.5),bty="u")
plot(0, type='n',xlab="", ylab="",xlim=c(min(x),max(x)), ylim=c(min(y),max(y)))
lines(newx,gcurve, type = 'l', lty = 1, lwd = 3, col = gtmpmodcol)
lines(newx,gcurve-2*gcurvese, type = 'l', lty = selty, lwd = selwd, col = gtmpmodcolse)
lines(newx,gcurve+2*gcurvese, type = 'l', lty = selty, lwd = selwd, col = gtmpmodcolse)
# add derivative
gtmpmodderivcol=rgb(colmat[2,1], colmat[2,2], colmat[2,3], max = 255) ; col2=gtmpmodderivcol
gtmpmodderivcolse=rgb(colmat[4,1], colmat[4,2], colmat[4,3], maxColorValue=255)
par(new=TRUE)
plot(gdX,gdY,type="l",lwd = 3,lty = 1,col=gtmpmodderivcol,xaxt="n",yaxt="n",xlab="",ylab="",axes=F,ylim=c(min(glwr),max(gupr))) # not wide enough for CI - check
lines(gdX,gupr,lty=selty,lwd=selwd,col=gtmpmodderivcolse) #check
lines(gdX,glwr,lty=selty,lwd=selwd,col=gtmpmodderivcolse) #check
abline(h=0,lty = 3,lwd = 2, col="black") # add line at zero
axis(side=4,labels=FALSE)
dev.off()

## plot derivative with milestone points - Fig. 1d
# derivative
#set to 1 if you want to plot point
plotp1 = 1 ; plotp2 = 1 ;  plotp3 = 1 ; plotp4 = 0 ;
#set to 1 if you want to plot text and lines
mplotp1x = 1 ; mplotp1y = 0; mplotp2 = 0 ;  mplotp3 = 1 ; mplotp4 = 1 ; mplotp5 = 1 ;
pdf(paste(outputdirfigs,'fig1d_deriv_milestonepoints_m234',measname,'.pdf',sep=''),width=6,height=5)
par(mar=c(4,4,0.5,0.5),bty="l")
tmpylab=y2labmyel
plot(0, type='n', xlab=xlabage, ylab=tmpylab, xlim=c(min(gdX),max(gdX)), ylim=c(min(glwr),max(gupr)))#ylim=c(min(Y),max(Y)))
lines(gdX,gdY, type = 'l', lty = 1, lwd = 3, col = gtmpmodderivcol)
lines(gdX,gupr,lty=selty,lwd=selwd,col=gtmpmodderivcolse) #check
lines(gdX,glwr,lty=selty,lwd=selwd,col=gtmpmodderivcolse) #check
abline(h=0,lty = 1,lwd = 1, col="black") # add line at zero
# add milestone markers
# first point
if (mplotp1x == 1) {
  if (plotp1 == 1) {
    if (plotp4 == 1) {
      pointcol1=rgb(colmat[5,1], colmat[5,2], colmat[5,3], max = 255)
    } else {
      pointcol1=rgb(colmat[6,1], colmat[6,2], colmat[6,3], max = 255)
    }
  }
  # add line down to x axis
  pointlinecol="gray"
  poinlinelw=2
  axislims=par("usr")
  tmpsty=mean(c(axislims[3],min(glwr)))
  tmpypos=tmpsty/2 #hight of text 
  tmpyposdist=abs(tmpypos*0.21) #dist from text above (white space above)
  tmpseqy1=seq(tmpsty,tmpypos-tmpyposdist,0.001)
  tmpseqy2=seq(tmpypos,gmaxderiv,0.001)
  tmpseqx1=rep(gageatmaxderiv,1,length(tmpseqy1))
  tmpseqx2=rep(gageatmaxderiv,1,length(tmpseqy2))
  lines(tmpseqx1,tmpseqy1, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  lines(tmpseqx2,tmpseqy2, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 1, at = gageatmaxderiv,labels = "",col=pointcol1)
  if (mplotp1y == 1) {
    text(x=gageatmaxderiv, y = tmpypos, pos=1, labels = "M2", col = pointcol1)
  } else {
    text(x=gageatmaxderiv, y = tmpypos, pos=1, labels = "M1", col = pointcol1)
  }
}
if (mplotp1y == 1) {
  # add line across to y axis
  tmpstx=mean(c(axislims[1],min(gdX)))
  tmpseqx=seq(tmpstx-0.5,gageatmaxderiv,0.1)
  tmpseqy=rep(gmaxderiv,1,length(tmpseqx))
  lines(tmpseqx,tmpseqy, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 2, at = gmaxderiv,labels = "",col=pointcol1)
  text(x=min(gdX), y = gmaxderiv, pos=1, labels = "M1", col = pointcol1,offset=-0.15)
}
# add first point
if (plotp1 == 1) {
  points(gageatmaxderiv,gmaxderiv,col=pointcol1,lw=5) # add point after to hide some of the line at the point
}
# second point
if (mplotp2 == 1) {
  if (plotp4 == 1) {
    pointcol2=rgb(colmat[6,1], colmat[6,2], colmat[6,3], max = 255)
  } else {
    pointcol2=rgb(colmat[7,1], colmat[7,2], colmat[7,3], max = 255)
  }
  # add line down to x axis
  tmpypos=tmpsty/2 #hight of text 
  tmpyposdist=abs(tmpypos*0.21) #dist from text above (white space above)
  tmpsty=mean(c(axislims[3],min(glwr)))
  
  #tmpseqy=seq(axislims[3],gderivvalatageatfirstzerocrossing,0.001)
  tmpseqy1=seq(tmpsty,tmpypos-tmpyposdist,0.001)
  tmpseqx1=rep(gageatfirstzerocrossing,1,length(tmpseqy1))
  tmpseqy2=seq(tmpypos,gderivvalatageatfirstzerocrossing,0.001)
  tmpseqx2=rep(gageatfirstzerocrossing,1,length(tmpseqy2))
  lines(tmpseqx1,tmpseqy1, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  lines(tmpseqx2,tmpseqy2, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 1, at = gageatfirstzerocrossing,labels = "",col=pointcol2)
  if (mplotp1y == 1) {
    text(x=gageatfirstzerocrossing, y = tmpsty/2, pos=1, labels = "M3", col = pointcol2)
  } else {
    text(x=gageatfirstzerocrossing, y = tmpsty/2, pos=1, labels = "M2", col = pointcol2)
  }
}
# add second point
if (plotp2 == 1) {
  if (plotp4 == 1) {
    pointcol2=rgb(colmat[6,1], colmat[6,2], colmat[6,3], max = 255)
  } else {
    pointcol2=rgb(colmat[7,1], colmat[7,2], colmat[7,3], max = 255)
  }
  points(gageatfirstzerocrossing,gderivvalatageatfirstzerocrossing,col=pointcol2,lw=5)
}
# third point
if (mplotp3 == 1) {
  if (plotp4 == 1) {
    pointcol3=rgb(colmat[7,1], colmat[7,2], colmat[7,3], max = 255)
  } else {
    pointcol3=rgb(colmat[8,1], colmat[8,2], colmat[8,3], max = 255)
  }
  # add line down to x axis
  tmpsty=mean(c(axislims[3],min(glwr)))
  #tmpseqy=seq(axislims[3],gderivvalatageatsecondzerocrossing[1],0.001)
  tmpypos=tmpsty/2 #hight of text 
  tmpyposdist=abs(tmpypos*0.21) #dist from text above (white space above)
  tmpseqy1=seq(tmpsty,tmpypos-tmpyposdist,0.001)
  tmpseqx1=rep(gageatsecondzerocrossing,1,length(tmpseqy1))
  tmpseqy2=seq(tmpypos,gderivvalatageatsecondzerocrossing,0.001)
  tmpseqx2=rep(gageatsecondzerocrossing,1,length(tmpseqy2))
  lines(tmpseqx1,tmpseqy1, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  lines(tmpseqx2,tmpseqy2, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 1, at = gageatsecondzerocrossing,labels = "",col=pointcol2)
  if (mplotp1y == 1) {
    text(x=gageatsecondzerocrossing, y = tmpsty/2, pos=1, labels = "M4", col = pointcol3)
  } else {
    text(x=gageatsecondzerocrossing, y = tmpsty/2, pos=1, labels = "M3", col = pointcol3)
  }
}
# add third point
if (plotp3 == 1) {
  points(gageatsecondzerocrossing,gderivvalatageatsecondzerocrossing,col=pointcol3,lw=5)
}
# fourth point
pointcol4=rgb(colmat[8,1], colmat[8,2], colmat[8,3], max = 255)
if (mplotp4 == 1) {
  # add line across to y axis
  tmpstx=mean(c(axislims[1],min(gdX)))
  tmpseqx=seq(tmpstx-0.5,gageat50pdecline,0.1)
  tmpseqy=rep(g50pdecline,1,length(tmpseqx))
  lines(tmpseqx,tmpseqy, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 2, at = g50pdecline,labels = "",col=pointcol4)
  text(x=min(gdX), y = g50pdecline, pos=1, labels = "M5", col = pointcol4,offset=-0.15)
}
# add fourth point
if (plotp4 == 1) {
  points(gageat50pdecline,g50pdecline,col=pointcol4,lw=5) # add point after to hide some of the line at the point
}
dev.off()

## plot derivative with milestone points - publication style (w/o labels)
# derivative
#set to 1 if you want to plot point
plotp1 = 1 ; plotp2 = 0 ;  plotp3 = 1 ; plotp4 = 1 ;
#set to 1 if you want to plot text and lines
mplotp1x = 0 ; mplotp1y = 1; mplotp2 = 0 ;  mplotp3 = 0 ; mplotp4 = 1 ; mplotp5 = 1 ;
pdf(paste(outputdirfigs,'fig1d_pub_deriv_milestonepoints_m234',measname,'.pdf',sep=''),width=6,height=5)
par(mar=c(0.6,0.6,0.5,0.5),bty="l")
tmpylab=y2labmyel
plot(0, type='n', xlim=c(min(gdX),max(gdX)), ylim=c(min(glwr),max(gupr)),xlab="",ylab="",yaxt="n",xaxt="n") ; axis(1,labels=FALSE) ;  axis(2,labels=FALSE) ;  #w/o labels (for figures)
lines(gdX,gdY, type = 'l', lty = 1, lwd = 3, col = gtmpmodderivcol)
lines(gdX,gupr,lty=selty,lwd=selwd,col=gtmpmodderivcolse) #check
lines(gdX,glwr,lty=selty,lwd=selwd,col=gtmpmodderivcolse) #check
abline(h=0,lty = 1,lwd = 1, col="black") # add line at zero
# add milestone markers
# first point
if (mplotp1x == 1) {
  if (plotp1 == 1) {
    if (plotp4 == 1) {
      pointcol1=rgb(colmat[5,1], colmat[5,2], colmat[5,3], max = 255)
    } else {
      pointcol1=rgb(colmat[6,1], colmat[6,2], colmat[6,3], max = 255)
    }
  }
  # add line down to x axis
  pointlinecol="gray"
  poinlinelw=2
  axislims=par("usr")
  tmpsty=mean(c(axislims[3],min(glwr)))
  tmpypos=tmpsty/2 #hight of text 
  tmpyposdist=abs(tmpypos*0.21) #dist from text above (white space above)
  tmpseqy1=seq(tmpsty,tmpypos-tmpyposdist,0.001)
  tmpseqy2=seq(tmpypos,gmaxderiv,0.001)
  tmpseqx1=rep(gageatmaxderiv,1,length(tmpseqy1))
  tmpseqx2=rep(gageatmaxderiv,1,length(tmpseqy2))
  lines(tmpseqx1,tmpseqy1, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  lines(tmpseqx2,tmpseqy2, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 1, at = gageatmaxderiv,labels = "",col=pointcol1)
  if (mplotp1y == 1) {
    text(x=gageatmaxderiv, y = tmpypos, pos=1, labels = "", col = pointcol1)
  } else {
    text(x=gageatmaxderiv, y = tmpypos, pos=1, labels = "", col = pointcol1)
  }
}
if (mplotp1y == 1) {
  # add line across to y axis
  tmpstx=mean(c(axislims[1],min(gdX)))
  #tmpseqx=seq(axislims[1],gageatmaxderiv,0.1)
  tmpseqx=seq(tmpstx-0.5,gageatmaxderiv,0.1)
  tmpseqy=rep(gmaxderiv,1,length(tmpseqx))
  lines(tmpseqx,tmpseqy, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 2, at = gmaxderiv,labels = "",col=pointcol1)
  text(x=min(gdX), y = gmaxderiv, pos=1, labels = "M1", col = pointcol1,offset=-0.15)
}
# add first point
if (plotp1 == 1) {
  points(gageatmaxderiv,gmaxderiv,col=pointcol1,lw=5) # add point after to hide some of the line at the point
}
# second point
if (mplotp2 == 1) {
  if (plotp4 == 1) {
    pointcol2=rgb(colmat[6,1], colmat[6,2], colmat[6,3], max = 255)
  } else {
    pointcol2=rgb(colmat[7,1], colmat[7,2], colmat[7,3], max = 255)
  }
  # add line down to x axis
  tmpypos=tmpsty/2 #hight of text 
  tmpyposdist=abs(tmpypos*0.21) #dist from text above (white space above)
  tmpsty=mean(c(axislims[3],min(glwr)))
  
  #tmpseqy=seq(axislims[3],gderivvalatageatfirstzerocrossing,0.001)
  tmpseqy1=seq(tmpsty,tmpypos-tmpyposdist,0.001)
  tmpseqx1=rep(gageatfirstzerocrossing,1,length(tmpseqy1))
  tmpseqy2=seq(tmpypos,gderivvalatageatfirstzerocrossing,0.001)
  tmpseqx2=rep(gageatfirstzerocrossing,1,length(tmpseqy2))
  lines(tmpseqx1,tmpseqy1, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  lines(tmpseqx2,tmpseqy2, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 1, at = gageatfirstzerocrossing,labels = "",col=pointcol2)
  if (mplotp1y == 1) {
    text(x=gageatfirstzerocrossing, y = tmpsty/2, pos=1, labels = "", col = pointcol2)
  } else {
    text(x=gageatfirstzerocrossing, y = tmpsty/2, pos=1, labels = "", col = pointcol2)
  }
}
# add second point
if (plotp2 == 1) {
  points(gageatfirstzerocrossing,gderivvalatageatfirstzerocrossing,col=pointcol2,lw=5)
}
# third point
if (mplotp3 == 1) {
  if (plotp4 == 1) {
    pointcol3=rgb(colmat[7,1], colmat[7,2], colmat[7,3], max = 255)
  } else {
    pointcol3=rgb(colmat[8,1], colmat[8,2], colmat[8,3], max = 255)
  }
  # add line down to x axis
  tmpsty=mean(c(axislims[3],min(glwr)))
  tmpypos=tmpsty/2 #hight of text 
  tmpyposdist=abs(tmpypos*0.21) #dist from text above (white space above)
  tmpseqy1=seq(tmpsty,tmpypos-tmpyposdist,0.001)
  tmpseqx1=rep(gageatsecondzerocrossing,1,length(tmpseqy1))
  tmpseqy2=seq(tmpypos,gderivvalatageatsecondzerocrossing,0.001)
  tmpseqx2=rep(gageatsecondzerocrossing,1,length(tmpseqy2))
  lines(tmpseqx1,tmpseqy1, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  lines(tmpseqx2,tmpseqy2, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 1, at = gageatsecondzerocrossing,labels = "",col=pointcol2)
  if (mplotp1y == 1) {
    #text(x=gageatmaxderiv, y = tmpsty, pos=1, labels = "M2", col = pointcol1) # does not currently work in the saved version
    text(x=gageatsecondzerocrossing, y = tmpsty/2, pos=1, labels = "", col = pointcol3)
  } else {
    text(x=gageatsecondzerocrossing, y = tmpsty/2, pos=1, labels = "", col = pointcol3)
  }
}
# add third point
if (plotp3 == 1) {
  points(gageatsecondzerocrossing,gderivvalatageatsecondzerocrossing,col=pointcol3,lw=5)
}
# fourth point
pointcol4=rgb(colmat[8,1], colmat[8,2], colmat[8,3], max = 255)
if (mplotp4 == 1) {
  # add line across to y axis
  tmpstx=mean(c(axislims[1],min(gdX)))
  tmpseqx=seq(tmpstx-0.5,gageat50pdecline,0.1)
  tmpseqy=rep(g50pdecline,1,length(tmpseqx))
  lines(tmpseqx,tmpseqy, type = 'l', lty = 3, lwd = poinlinelw, col = pointlinecol)
  # add tick mark and text
  axis(side = 2, at = g50pdecline,labels = "",col=pointcol4)
  text(x=min(gdX), y = g50pdecline, pos=1, labels = "", col = pointcol4,offset=-0.15)
}
# add fourth point
if (plotp4 == 1) {
  points(gageat50pdecline,g50pdecline,col=pointcol4,lw=5) # add point after to hide some of the line at the point
}
dev.off()


## save milestones measures # please note that the M1, M2, and M3 mentioned in the paper, are called m2, m3, and m4 here, while the rate of maturation and decline, are called m1 and m5, respectively
gmilestones=paste("the milestones were ",gmaxderiv," aribtrary units, ",gageatmaxderiv," years, ",gageatfirstzerocrossing, " years, and ",gageatsecondzerocrossing," years, respectively, and the rate of decline at 50% between age at start of decline and max age (at ",gageat50pdecline," years), was ",g50pdecline,sep="")
write(gmilestones, file = paste(outputdirsummeas,'global_milestones','.galpha',galpha,'.',measname,'.txt', sep = "")) 
#write m2 and m3 for SC in matlab
write(gageatmaxderiv, file = paste(outputdirsummeas,'global_milestone2','.galpha',galpha,'.',measname,'.txt', sep = "")) 
write(gageatfirstzerocrossing, file = paste(outputdirsummeas,'global_milestone3','.galpha',galpha,'.',measname,'.txt', sep = "")) 
write(gageatsecondzerocrossing, file = paste(outputdirsummeas,'global_milestone4','.galpha',galpha,'.',measname,'.txt', sep = "")) 
write(g50pdecline, file = paste(outputdirsummeas,'global_milestone5','.galpha',galpha,'.',measname,'.txt', sep = "")) 
write(gageat50pdecline, file = paste(outputdirsummeas,'global_milestone6','.galpha',galpha,'.',measname,'.txt', sep = ""))


## fit splines per ROI ##
semipar.output = semipar.mix.mp(Y = Y, x = age, param = NULL, random = NULL, k = k, data.ran = NULL, norder = norder0, pen.order = pen.order0, knots = "quantile", store.gamm4 = TRUE)

# extract statistics, and curves for plotting
p.spl = pwdf = pwlsp = vector(length = nroi)
curve.all = curve.se.all = matrix(nrow = nroi, ncol = framelength)
deriv.all = deriv.se.all = matrix(nrow = nroi, ncol = framelength)

# set variables for max age, etc
maxestval  = tmpmaxderivindabs = tmpmaxderivabs = tmpmaxderiv = tmpageatmaxderiv = tmpageatfirstzerocrossing = tmpageatsecondzerocrossing = tmp50pdecline = tmpageat50pdecline = vector(length = nroi)
xvalatyindclosestoprctofi = yvatatyindclosestoprctofi = maxd1 = ageatmaxd1 = rateofchangeatprcofi = xvalatyindclosestoprctofiaging = rateofchangeatprcofiaging = vector(length = nroi)

#set name of output file with splines
pdf(paste(outputdirsummeas,'splines_wsummeas',measname,'alpha',toString(alpha),'.pdf',sep=''),width=6,height=5)

for (n in 1:nroi) {
  
  # p-values
  p.spl[n] = summary(semipar.output$gamm4.list[[n]]$gam)$s.pv
  pwdf[n] = sum(semipar.output$gamm4.list[[n]]$gam$edf) # pointwise effective degrees of freedom 
  pwlsp[n] = log(semipar.output$gamm4.list[[n]]$gam$sp) # pointwise log smoothing parameters
  
  # prepare splines for plotting using predict.gam
  tmpmod = semipar.output$gamm4.list[[n]]$gam #tmp model
  fake.frame <- data.frame(x=seq(min(age),max(age),length=framelength) ) # might included "Female" as above
  data.fit = predict(tmpmod, newdata = fake.frame, se.fit=T)
  newx=gfake.frame$x
  tmpcurve = data.fit$fit
  tmpcurvese = data.fit$se
  
  # store curves and knot locations for later plotting
  curve.all[n,] = tmpcurve
  curve.se.all[n,] = tmpcurvese
  
  # derivative (based on Deriv function http://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/) 
  tmpdata.fit.lpmat <- predict(tmpmod,newdata = X,type="lpmatrix") #same as above, but now saved as matrix with "lpmatrix"
  tmpdYm = diff(tmpdata.fit.lpmat)/diff(newx)  # the derivative of the lpmatrix
  tmpdY = tmpdYm %*% coef(tmpmod) # multiply the mat from above with the coefficients to get one vector
  tmpdX <- rowMeans(embed(newx,2)) # centers the X values for plotting
  # se.derivative
  tmpmodvp = tmpmod$Vp
  tmpdY.sd <- rowSums(tmpdYm %*% tmpmodvp * tmpdYm)^.5
  tmpresidual.df <- df.residual(tmpmod)
  tmptVal <- qt(1 - (alpha/2), tmpresidual.df)
  tmpupr <- tmpdY + tmptVal * tmpdY.sd
  tmplwr <- tmpdY - tmptVal * tmpdY.sd
  
  # get different properties of the curve
  # max rate of change (1)
  tmpmaxderiv[n] = max(tmpdY) # save for brain plotting and histograms
  
  # age at max rate of change (2)
  tmpmaxderivind = which.max(tmpdY)
  tmpmaxderivindabs[n] = which.max(abs(tmpdY))
  tmpmaxderivabs[n]=abs(tmpdY[tmpmaxderivindabs[n]])
  tmpageatmaxderiv[n] = newx[tmpmaxderivind] # save for brain plotting and histograms
  tmpfuncvalatmaxderiv = tmpcurve[tmpmaxderivind]
  # age at 1st zero crossing (age at where lower deriv.se crosses zero first) - only works if inverted-U shaped?
  # first, find first point where change is sig (i.e. CI not containing zero)
  tmpnewind=1
  while (cbind(tmplwr[tmpnewind] <= 0 & tmpupr[tmpnewind] >= 0)) {tmpnewind=tmpnewind+1 ; if (!isTRUE(tmpnewind <= length(tmplwr))) break; if (!cbind(tmplwr[tmpnewind] <= 0 & tmpupr[tmpnewind] >= 0)) break; }
  print(tmpnewind)
  #second, find where changes from pos to neg (given that it's an inverted-U shape)
  # lwr
  setderivmeasofi=tmplwr # deriv measure of interest
  # find were the derivmeasofi goes from pos to neg (i.e. hits/crosses zero)
  tmpval=setderivmeasofi[tmpnewind] #this only works if curve is inverted U-shaped, and not if U-shaped (then you would want the upperderivse)
  while (isTRUE(tmpval>=0)) {tmpnewind=tmpnewind+1 ; tmpval=setderivmeasofi[tmpnewind]; if (!isTRUE(tmpval>=0)) break; }
  print(tmpnewind)
  tmpageatlastpos = newx[tmpnewind]
  tmpageatlastpospluss1 = newx[tmpnewind+1]
  tmpageatfirstzerocrossing[n]=mean(c(tmpageatlastpos,tmpageatlastpospluss1)) # save for brain plotting and histograms
  tmpderivvalatageatfirstzerocrossing=mean(c(setderivmeasofi[tmpnewind],setderivmeasofi[tmpnewind+1]))
  tmpvalatageatfirstzerocrossing=mean(c(tmpcurve[tmpnewind],tmpcurve[tmpnewind+1]))
  
  # age at 2nd zero crossing (age at where upper deriv.se crosses zero first) - only works if inverted-U shaped?
  # upr
  setderivmeasofi=tmpupr # deriv measure of interest
  # find were the derivmeasofi goes from pos to neg (i.e. hits/crosses zero)
  tmpval=setderivmeasofi[tmpnewind] #this only works if curve is inverted U-shaped, and not if U-shaped (then you would want the upperderivse)
  while (isTRUE(tmpval>=0)) {tmpnewind=tmpnewind+1 ; tmpval=setderivmeasofi[tmpnewind]; if (!isTRUE(tmpval>=0)) break; }
  print(tmpnewind)
  if (tmpnewind==length(newx)) {
    print("no second crossing?") 
    tmpageatsecondzerocrossing[n]=newx[tmpnewind] # save for brain plotting and histograms
    tmpderivvalatageatsecondzerocrossing=mean(c(setderivmeasofi[tmpnewind],setderivmeasofi[tmpnewind+1]))
    tmpvalatageatsecondzerocrossing=mean(tmpcurve[tmpnewind]) #mean is a fix to get a non-matrix vluea
  } else {
    tmpageatlastpos2 = newx[tmpnewind]
    tmpageatlastpos2pluss1 = newx[tmpnewind+1]
    tmpageatsecondzerocrossing[n]=mean(c(tmpageatlastpos2,tmpageatlastpos2pluss1)) # save for brain plotting and histograms
    tmpderivvalatageatsecondzerocrossing=mean(c(setderivmeasofi[tmpnewind],setderivmeasofi[tmpnewind+1]))
    tmpvalatageatsecondzerocrossing=mean(c(tmpcurve[tmpnewind],tmpcurve[tmpnewind+1]))
  }
  
  #deriv 50% after 2nd zero crossing and age at this value (decline, if decline) after 2nd zero crossing 
  if (tmpnewind==length(newx)) {
    print("same as 2nd zero crossing") 
    tmpageat50pdecline[n]=newx[tmpnewind] # save for brain plotting and histograms
    tmp50pdecline[n]=tmpdY[length(tmpdY)] # save for brain plotting and histograms
  } else {
    tmpstart1=tmpnewind
    tmpend1=length(tmpdY)
    tmp50percfromsecondcrossingind=tmpstart1+((tmpend1-tmpstart1)/2)
    tmp50pdecline[n]=tmpdY[tmp50percfromsecondcrossingind] # save for brain plotting and histograms
    tmpageat50pdecline[n] = newx[tmp50percfromsecondcrossingind] # save for brain plotting and histograms
    tmpvalat50pdecline = tmpcurve[tmp50percfromsecondcrossingind]
  }
  
  ## plot
  # plot spline - no milestone points
  # set name and colors
  tmplabelname = labelname[n,]
  print(tmplabelname)
  tmprgb=rgb(rgbrois[n,1], rgbrois[n,2], rgbrois[n,3], maxColorValue=1)
  par(mar=c(4,4,2,0.5),bty="l")
  y=Y[,n]
  x=age
  plot(0, type='n', xlab='age (years)', ylab="T1w/T2w (a.u.)", xlim=c(min(x),max(x)), ylim=c(min(y),max(y)))#ylim=c(min(Y),max(Y)))
  points(x,y, col = 'black', pch = 16, cex = 0.75) # plot data#title('maturational trajectory matrix', line=+2, font.main=1)
  tmpmodcol=tmprgb
  tmpmodcolse=rgb(rgbrois[n,1], rgbrois[n,2], rgbrois[n,3], 0.6,names = NULL,maxColorValue=1)
  selwd=1.5
  selty=1
  lines(newx,tmpcurve, type = 'l', lty = 1, lwd = 3, col = tmpmodcol)
  lines(newx,tmpcurve-2*tmpcurvese, type = 'l', lty = selty, lwd = selwd, col = tmpmodcolse)
  lines(newx,tmpcurve+2*tmpcurvese, type = 'l', lty = selty, lwd = selwd, col = tmpmodcolse)
  title(tmplabelname)
  
  #plot spline + deriv in one plot
  par(mar=c(4,4,2,4),bty="u")
  plot(0, type='n',xlab='age (years)', ylab="",xlim=c(min(x),max(x)), ylim=c(min(y),max(y)))#ylim=c(min(Y),max(Y)))
  mtext("T1w/T2w (a.u.)",side=2,line=3,col=tmpmodcol)
  points(x,y, col = 'black', pch = 16, cex = 0.75) # plot data#title('maturational trajectory matrix', line=+2, font.main=1)
  lines(newx,tmpcurve, type = 'l', lty = 1, lwd = 3, col = tmpmodcol)
  lines(newx,tmpcurve-2*tmpcurvese, type = 'l', lty = selty, lwd = selwd, col = tmpmodcolse)
  lines(newx,tmpcurve+2*tmpcurvese, type = 'l', lty = selty, lwd = selwd, col = tmpmodcolse)
  title(tmplabelname)
  
  # add milestone points (spline)
  # first point
  pointcol1=rgb(colmat[5,1], colmat[5,2], colmat[5,3], max = 255) # same as for global
  points(tmpageatmaxderiv[n],tmpfuncvalatmaxderiv,col=pointcol1,lw=5)
  # second point
  pointcol2=rgb(colmat[6,1], colmat[6,2], colmat[6,3], max = 255)
  points(tmpageatfirstzerocrossing[n],tmpvalatageatfirstzerocrossing,col=pointcol2,lw=5)
  # third point
  pointcol3=rgb(colmat[7,1], colmat[7,2], colmat[7,3], max = 255)
  points(tmpageatsecondzerocrossing[n],tmpvalatageatsecondzerocrossing,col=pointcol3,lw=5)
  # fourth point - only if diff from third
  if (tmpageat50pdecline[n] != tmpageatsecondzerocrossing[n]) {
    pointcol4=rgb(colmat[8,1], colmat[8,2], colmat[8,3], max = 255)
    points(tmpageat50pdecline[n],tmpvalat50pdecline,col=pointcol4,lw=5)
  }
  
  # add derivative
  #gtmpmodderivcol="red"
  #gtmpmodderivcol=rgb(colmat[2,1], colmat[2,2], colmat[2,3], max = 255)
  #gtmpmodderivcolse=rgb(colmat[4,1], colmat[4,2], colmat[4,3], maxColorValue=255)
  tmpmodderivcol=gtmpmodderivcol
  tmpmodderivcolse=gtmpmodderivcolse
  par(new=TRUE)
  plot(tmpdX,tmpdY,type="l",lwd = 3,lty = 1,col=tmpmodderivcol,xaxt="n",yaxt="n",xlab="",ylab="",axes=F,ylim=c(min(glwr),max(gupr))) # not wide enough for CI - check
  lines(tmpdX,tmpupr,lty=selty,lwd=selwd,col=tmpmodderivcolse) #check
  lines(tmpdX,tmplwr,lty=selty,lwd=selwd,col=tmpmodderivcolse) #check
  abline(h=0,lty = 1,lwd = 1, col="black") # add line at zero
  axis(side=4)
  mtext("Rate of change",side=4,line=3,col=gtmpmodderivcol)
}
dev.off()


#save milestones
write(tmpmaxderiv, file = paste(outputdirsummeas,'m1','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpageatfirstzerocrossing)) #as row vector
write(tmpageatmaxderiv, file = paste(outputdirsummeas,'m2','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpageatfirstzerocrossing)) #as row vector
write(tmpageatfirstzerocrossing, file = paste(outputdirsummeas,'m3','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpageatfirstzerocrossing)) #as row vector
write(tmpageatsecondzerocrossing, file = paste(outputdirsummeas,'m4','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpageatfirstzerocrossing)) #as row vector
write(tmp50pdecline, file = paste(outputdirsummeas,'m5','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpageatfirstzerocrossing)) #as row vector
write(tmpageat50pdecline, file = paste(outputdirsummeas,'m6','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpageatfirstzerocrossing)) #as row vector


## write to file for brain plotting
#have header and save all to one file?
#- max rate of change
write(tmpmaxderiv, file = paste(outputdirsummeas,'maxrateofchange','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector
#- ageatmaxestimatedvalue
write(tmpageatmaxderiv, file = paste(outputdirsummeas,'ageatmaxrateofchange','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector
#ageatfirstzerocrossing
write(tmpageatfirstzerocrossing, file = paste(outputdirsummeas,'ageatfirstzerocrossing','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(maxestval)) #as row vector
#ageatsecondzerocrossing
write(tmpageatsecondzerocrossing, file = paste(outputdirsummeas,'ageatsecondzerocrossing','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(maxestval)) #as row vector
#50pdecline
write(tmp50pdecline, file = paste(outputdirsummeas,'rateofchangedecline','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(maxestval)) #as row vector

#write matrices to file for correlation
write(t(curve.all), file = paste(outputdir,'curvesall.',measname,'.txt', sep = ""), ncolumns=dim(curve.all)[2])
write(t(deriv.all), file = paste(outputdirsummeas,'derivcurvesall',measname,'.txt', sep = ""), ncolumns=dim(deriv.all)[2])
#write mean across subjects as well
write(colMeans(Y), file = paste(outputdirsummeas,'meanxsubjects.',measname,'.txt', sep = ""), ncolumns=length(colMeans(Y)))


#all milestones matrix - 4 (all)
tmpscatterlabels=c("Peak growth rate","Peak growth age","Onset stability age","Onset decline age","Decline rate")
nmiles=5
milestones.all = matrix(nrow = nmiles , ncol = nroi)
milestones.all[1,] = tmpmaxderiv
milestones.all[2,] = tmpageatmaxderiv #NB! NOT Age first to get age at x-axis
milestones.all[3,] =tmpageatfirstzerocrossing
milestones.all[4,] =tmpageatsecondzerocrossing
milestones.all[5,]=tmp50pdecline

#first, reorder, putting M1 first, and transpose
milestones1to5 = milestones.all[1:5,]

M=t(milestones1to5)
#second, add colum and row names,
corrlabels=c("Peak growth","Peak growth age","Onset stability age","Onset decline age","Decline rate")
colnames(M) <- corrlabels
rownames(M) <- 1:dim(M)[1]


##correlate derivative measures (M1 and M5) with mean T1w/T2w levels, and with each other
#m1-avg
yname="maxderiv";xname="avgmyelin"
tmpcor=rcorr(tmpmaxderiv,colMeans(Y),type="spearman")
rhoandpval=paste("Spearman's r = ",tmpcor$r[1,2],", p = ",tmpcor$P[1,2],sep="")
write(rhoandpval, file = paste(outputdirsummeas,'rhoandpval_',xname,'_',yname,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#m5-avg
yname="50pdecline";xname="avgmyelin"
tmpcor=rcorr(tmp50pdecline,colMeans(Y),type="spearman")
rhoandpval=paste("Spearman's r = ",tmpcor$r[1,2],", p = ",tmpcor$P[1,2],sep="")
write(rhoandpval, file = paste(outputdirsummeas,'rhoandpval_',xname,'_',yname,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#m1-m5
yname="maxderiv";xname="50pdecline"
tmpcor=rcorr(tmp50pdecline,colMeans(Y),type="spearman")
rhoandpval=paste("Spearman's r = ",tmpcor$r[1,2],", p = ",tmpcor$P[1,2],sep="") ; print(rhoandpval)
write(rhoandpval, file = paste(outputdirsummeas,'rhoandpval_',xname,'_',yname,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

# mixture distribution
rseed <- .Random.seed
tmpmeass=c("M1","M2","M2b","M3","M4","M4b","M5")
tmpclustclassification=matrix(nrow=length(tmpmeass),ncol=length(tmpmaxderiv))
for (ti in 1:length(tmpmeass)) {
  tmpmeasn=tmpmeass[ti]
  if (tmpmeasn == "M1") {
    tmpmeas=tmpmaxderiv
    m1mclust = densityMclust(tmpmeas)
    tmpclust=m1mclust
    tmpxlab="Rate of change - myelination"
    ylim2=120
    cutage = 0
  } else if (tmpmeasn == "M5") {
    tmpmeas=tmp50pdecline
    m5mclust = densityMclust(tmpmeas)
    tmpclust=m5mclust
    tmpxlab="Rate of change - decline"
    ylim2=80
    cutage = 0
  } else if (tmpmeasn == "M2") {
    tmpmeas=tmpageatmaxderiv
    m2mclust = densityMclust(tmpmeas,G=2)
    tmpclust=m2mclust
    tmpxlab="Age (years)"
    ylim2=0.3
    cutage = 14
  } else if (tmpmeasn == "M2b") {
    tmpmeas=tmpageatmaxderiv[tmpageatmaxderiv<25]
    m2bmclust = densityMclust(tmpmeas,G=2)
    tmpclust=m2bmclust
    tmpxlab="Age (years)"
    ylim2=0.3
    cutage = 14
  } else if (tmpmeasn == "M3") {
    tmpmeas=tmpageatfirstzerocrossing
    m3mclust = densityMclust(tmpmeas)
    tmpclust=m3mclust
    tmpxlab="Age (years)"
    ylim2=0.3
    #cutage = 38 #this should probably be done in a none eyeball manner by looking at where the distributions cross above (mixmdl)
    cutage = 39.5 #this should probably be done in a none eyeball manner by looking at where the distributions cross above (mixmdl)
  } else if (tmpmeasn == "M4") {
    tmpmeas=tmpageatsecondzerocrossing
    m4mclust = densityMclust(tmpmeas)
    tmpclust=m4mclust
    tmpxlab="Age (years)"
    #ylim2=1.2
    ylim2=0.3
  } else if (tmpmeasn == "M4b") {
    tmpm4nolowoutlier = tmpageatsecondzerocrossing[tmpageatsecondzerocrossing > 50]
    tmpmeas = tmpm4nolowoutlier
    m4mclustnolowoutlier = densityMclust(tmpmeas,G = 2)
    tmpclust=m4mclustnolowoutlier
    tmpxlab="Age (years)"
    #ylim2=1.2
    ylim2=0.3
    cutage = 80
  }
  print(tmpmeasn)
  
  #write results to files
  write(tmpclust$classification, file = paste(outputdirsummeas,tmpmeasn,'_mcclustclassification','.alpha',alpha,'.',measname,'.txt', sep = ""),ncolumns = length(tmpclust$classification)) #as row vector
  write(tmpclust$G, file = paste(outputdirsummeas,tmpmeasn,'_mcclustNcomponents','.alpha',alpha,'.',measname,'.txt', sep = ""),ncolumns = length(tmpclust$G)) #as row vector
  
  #plot qqplot
  lwdn=5
  pdf(paste(outputdirfigs,'supfig9.qqplot',tmpmeasn,'.pdf',sep=''),width=6,height=5)
  par(mar=c(4,4,2,0.5),bty="l")
  #qqplot
  qqnorm(tmpmeas,lwd=lwdn)
  qqline(tmpmeas,lwd=lwdn,col="grey")
  dev.off()
  
  #run mixmdl
  maxrestartsn=200 ;
  mixmdl = normalmixEM(tmpmeas,maxrestarts=maxrestartsn) ;
  
  ax=1
  pdf(paste(outputdirfigs,'fig2b.mixmdlplot',tmpmeasn,'.pdf',sep=''),width=tmpscatterw,height=tmpscatterh)
  lwdnmix=5 ; ylab="Density" ; col11 = "lightgrey"; bordercol = "white" ; col22=c(col1, col2) ; nbreaks = 20 ; main2="Density Curves"
  par(mar=c(4,4,0.5,0.5),bty="l");
  plotmixemcustom(mixmdl,lwdnmix,ylim2,col11,col22,nbreaks,lwdn,tmpxlab,ylab,main2,ax,cutage)
  dev.off()
  
}


# M5 vs M1
fignum="3b" ; xname="maxderiv" ; yname="50pdecline"; 
xlabvar="Developmental rate of change (a.u)" ; ylabvar="Aging rate of change (a.u)"
plotscatter(tmpmaxderiv,tmp50pdecline,fignum,xname,yname,xlabvar,ylabvar) 
plotscatterspearman(tmpmaxderiv,tmp50pdecline,fignum,xname,yname,xlabvar,ylabvar) #don't really need both, scatter is the same, and spearman corr is calculated above


#Compute Hartigans' Dip Test Statistic for Unimodality
diptestp=vector(length=dim(milestones.all[1:5,])[1])
dipp.tantrump=vector(length=dim(milestones.all[1:5,])[1])
booststraplrts=vector(length=dim(milestones.all[1:5,])[1])
booststrapp=vector(length=dim(milestones.all[1:5,])[1])

for (mi in 1:dim(milestones.all[1:5,])[1]) {
  vals=milestones.all[mi,]
  diptest = dip.test(vals, simulate.p.value = TRUE, B = 5000)
  #diptest = dip.test(vals)
  #summary(diptest)
  diptestp[mi]=diptest$p.value #nb! tmpageatmaxderiv is first
  tmpdipp.tantrum = dipp.tantrum(vals,diptest$statistic,M=10000)
  dipp.tantrump[mi]=tmpdipp.tantrum$p.value
  
  #another test - log likelihood - from here: https://stats.stackexchange.com/questions/138223/how-to-test-if-my-distribution-is-multimodal
  #paramtric
  vals.gmm = densityMclust(vals)
  if (!vals.gmm$G==1) {
    #vs 1
    vals.gmm.1 = densityMclust(vals, G=1)
    logLik(vals.gmm.1)
    # 'log Lik.' -1114.274 (df=2) for mi=3 for 70%
    tmploglik = logLik(vals.gmm)-logLik(vals.gmm.1)
    # 'log Lik.' 118.3881 (df=8) for mi=3 for 70%
    tmppvalvs1=1-pchisq(tmploglik[1], df=vals.gmm$df-vals.gmm.1$df)  # 
    print(tmppvalvs1)
  }
  if (!vals.gmm$G==2) {
    #vs 2
    vals.gmm.2 = densityMclust(vals, G=2)
    logLik(vals.gmm.2)
    # 'log Lik.' -1114.274 (df=2) for mi=3 for 70%
    tmploglik = logLik(vals.gmm)-logLik(vals.gmm.2)
    # 'log Lik.' 118.3881 (df=8) for mi=3 for 70%
    tmppvalvs2=1-pchisq(tmploglik[1], df=vals.gmm$df-vals.gmm.2$df)  
    print(tmppvalvs2)
  }
  #Bootstrap Likelihood Ratio Test For The Number Of Mixture Components (using mclust)
  B = 10000;    
  vals.boot = mclustBootstrapLRT(vals, model = "E",nboot = B,maxG = 2) #"V" is univariate mixture, unequal variance, "E" is uni mix, equal var, see mclustModelNames for details
  #maxG is the maximum number of components tested, can have more, but if only want to test for 1 vs 2, set maxG = 2
  booststraplrts[mi] = vals.boot$obs[1] #the likelihood ratio test statistic (LRTS).
  booststrapp[mi] = vals.boot$p.value[1]
}

#write boot lrts and p to file
write(booststraplrts, file = paste(outputdirsummeas,'booststraplrts','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(booststraplrts)) #as row vector
write(booststrapp, file = paste(outputdirsummeas,'booststrapp','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(booststrapp)) #as row vector

#write p values to file
write(diptestp, file = paste(outputdirsummeas,'diptestp','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(diptestp)) #as row vector
write(dipp.tantrump, file = paste(outputdirsummeas,'dipp.tantrump','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(diptestp)) #as row vector

#test m2b and m4b as well
diptestp2 = dipp.tantrump2 = vector(length=2)
for (mii in 1:2) {
  if (mii == 1) {
    vals=tmpageatmaxderiv[tmpageatmaxderiv<25]
  } else {
    vals = tmpageatsecondzerocrossing[tmpageatsecondzerocrossing < 80 & tmpageatsecondzerocrossing > 50]
  }
  diptest = dip.test(vals, simulate.p.value = TRUE, B = 5000)
  #diptest = dip.test(vals)
  #summary(diptest)
  diptestp2[mii]=diptest$p.value #nb! tmpageatmaxderiv is first
  tmpdipp.tantrum = dipp.tantrum(vals,diptest$statistic,M=10000)
  dipp.tantrump2[mii]=tmpdipp.tantrum$p.value
}

#write p values to file
write(diptestp2, file = paste(outputdirsummeas,'diptestp2','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(diptestp)) #as row vector
write(dipp.tantrump2, file = paste(outputdirsummeas,'dipp.tantrump2','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(diptestp)) #as row vector


### Two-by-two plots contingency tables ###
catmilestones.all = catmilestones.R = catmilestones.L = matrix(nrow = nmiles , ncol = nroi) #"categorization" of milestones

# 2-Way Frequency Table - M2 & M3
tmpageatmaxderivval = 14
tmpageatfirstzerocrossingval = 39.5 #based on approximate identifiaction of where the distributions cross above (mixmdl)
tmpageatsecondzerocrossingval = 85

A=tmpageatmaxderiv > tmpageatmaxderivval
B=tmpageatfirstzerocrossing > tmpageatfirstzerocrossingval

## a mean and ci vals
aearlyvals = tmpageatmaxderiv[which(A==FALSE)]
alatevals = tmpageatmaxderiv[which(A==TRUE)]
cilev = 0.95 ;
ciearly = CI(aearlyvals, ci = cilev) ;ciearlyup = ciearly[[1]]; ciearlymean = ciearly[[2]]; ciearlylow = ciearly[[3]];
cilate= CI(alatevals, ci = cilev) ; cilateup = cilate[[1]]; cilatemean = cilate[[2]]; cilatelow = cilate[[3]];
ndec=2; write(paste("early wave (n= ",length(aearlyvals),") mean age of ",round(ciearlymean,ndec)," (95% CI = ",round(ciearlylow,ndec),"-",round(ciearlyup,ndec),")","\nlate wave (n= ",length(alatevals),") mean age of ",round(cilatemean,ndec)," (95% CI = ",round(cilatelow,ndec),"-",round(cilateup,ndec),")",sep=""), file = paste(outputdirsummeas,'m2_earlylatewavemeanci','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
## b mean and ci vals
bearlyvals = tmpageatfirstzerocrossing[which(B==FALSE)]
blatevals = tmpageatfirstzerocrossing[which(B==TRUE)]
cilev = 0.95 ;
ciearly = CI(bearlyvals, ci = cilev) ;ciearlyup = ciearly[[1]]; ciearlymean = ciearly[[2]]; ciearlylow = ciearly[[3]];
cilate= CI(blatevals, ci = cilev) ; cilateup = cilate[[1]]; cilatemean = cilate[[2]]; cilatelow = cilate[[3]];
ndec=2; write(paste("early wave (n= ",length(bearlyvals),") mean age of ",round(ciearlymean,ndec)," (95% CI = ",round(ciearlylow,ndec),"-",round(ciearlyup,ndec),")","\nlate wave (n= ",length(blatevals),") mean age of ",round(cilatemean,ndec)," (95% CI = ",round(cilatelow,ndec),"-",round(cilateup,ndec),")",sep=""), file = paste(outputdirsummeas,'m3_earlylatewavemeanci','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

A[A==TRUE]="M2-w2"
A[A==FALSE]="M2-w1"
B[B==TRUE]="M3-w2"
B[B==FALSE]="M3-w1"

catmilestones.all[1,]=A
catmilestones.all[2,]=B

#table 1 data
##get mean age and CI M1 and M5 vals of early and late M2&M3 wave (and put in data frame)
##m1
tmpmeas=tmpmaxderiv
#m1 vals for m2
tmpwavefac = A
tmpdframe = data.frame(tmpmeas,tmpwavefac)
#trad
m1meanci_m2w1w2_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
write.table(m1meanci_m2w1w2_trad, file = paste(outputdirsummeas,'table1_m1meanci_m2w1w2_trad','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#boot
m1meanci_m2w1w2_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
write.table(m1meanci_m2w1w2_boot, file = paste(outputdirsummeas,'table1_m1meanci_m2w1w2_boot','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector


#m1 vals for m3
tmpwavefac = B
tmpdframe = data.frame(tmpmeas,tmpwavefac)
#trad
m1meanci_m3w1w2_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
write.table(m1meanci_m3w1w2_trad, file = paste(outputdirsummeas,'table1_m1meanci_m3w1w2_trad','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#boot
m1meanci_m3w1w2_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
write.table(m1meanci_m3w1w2_boot, file = paste(outputdirsummeas,'table1_m1meanci_m3w1w2_boot','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

##m5
#m5 vals for m2
tmpmeas=tmp50pdecline
tmpwavefac = A
tmpdframe = data.frame(tmpmeas,tmpwavefac)
#trad
m5meanci_m2w1w2_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
write.table(m5meanci_m2w1w2_trad, file = paste(outputdirsummeas,'table1_m5meanci_m2w1w2_trad','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#boot
m5meanci_m2w1w2_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
write.table(m5meanci_m2w1w2_boot, file = paste(outputdirsummeas,'table1_m5meanci_m2w1w2_boot','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#m5 vals for m3
tmpwavefac = B
tmpdframe = data.frame(tmpmeas,tmpwavefac)
#trad
m5meanci_m3w1w2_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
write.table(m5meanci_m3w1w2_trad, file = paste(outputdirsummeas,'table1_m5meanci_m3w1w2_trad','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#boot
m5meanci_m3w1w2_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
write.table(m5meanci_m3w1w2_boot, file = paste(outputdirsummeas,'table1_m5meanci_m3w1w2_boot','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

#regression - results not reported (reported as mean and CI instead), but p-values used for multiple comparisons correction below
##regression with M1 and M5
m2factor =  factor(A)
#m2-m1
glmm1m2cat = lm(tmpmaxderiv ~ m2factor)
sumglmm1m2cat =summary(glmm1m2cat)
anovaglmm1m2cat = anova(glmm1m2cat)
#m2-m5
glmm5m2cat = lm(tmp50pdecline ~ m2factor)
sumglmm5m2cat =summary(glmm5m2cat)
anovaglmm5m2cat = anova(glmm5m2cat)
#m3-m1
Bfactor = factor(B)
m3factor=Bfactor #for use for von Economo later
glmm1m3cat = lm(tmpmaxderiv ~ m3factor)
sumglmm1m3cat = summary(glmm1m3cat)
anovaglmm1m3cat = anova(glmm1m3cat)
#m3-m5
glmm5m3cat = lm(tmp50pdecline ~ m3factor)
sumglmm5m3cat = summary(glmm5m3cat)
anovaglmm5m3cat = anova(glmm5m3cat)
#write results to file
#m1
m1m2reg=paste("Reg coeff = ",round(glmm1m2cat$coefficients[2],digits=4),", F(",anovaglmm1m2cat$Df[1],",",anovaglmm1m2cat$Df[2],") = ",round(anovaglmm1m2cat$`F value`[1],digits=1),", p = ",round(anovaglmm1m2cat$`Pr(>F)`[1],digits=4),", Rsquared = ",round(sumglmm1m2cat$r.squared,digits=3),sep="")
m1m3reg=paste("Reg coeff = ",round(glmm1m3cat$coefficients[2],digits=4),", F(",anovaglmm1m3cat$Df[1],",",anovaglmm1m3cat$Df[2],") = ",round(anovaglmm1m3cat$`F value`[1],digits=1),", p = ",round(anovaglmm1m3cat$`Pr(>F)`[1],digits=4),", Rsquared = ",round(sumglmm1m3cat$r.squared,digits=3),sep="")
#write(m1m2reg, file = paste(outputdirsummeas,'m1m2catreg','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#write(m1m3reg, file = paste(outputdirsummeas,'m1m3catreg','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#m5
m5m2reg=paste("Reg coeff = ",glmm5m2cat$coefficients[2],", F(",anovaglmm5m2cat$Df[1],",",anovaglmm5m2cat$Df[2],") = ",anovaglmm5m2cat$`F value`[1]," p = ",anovaglmm5m2cat$`Pr(>F)`[1],", Rsquared = ",sumglmm5m2cat$r.squared,sep="")
m5m3reg=paste("Reg coeff = ",glmm5m3cat$coefficients[2],", F(",anovaglmm5m3cat$Df[1],",",anovaglmm5m3cat$Df[2],") = ",anovaglmm5m3cat$`F value`[1]," p = ",anovaglmm5m3cat$`Pr(>F)`[1],", Rsquared = ",sumglmm5m3cat$r.squared,sep="")
#write(m5m2reg, file = paste(outputdirsummeas,'m5m2catreg','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#write(m5m3reg, file = paste(outputdirsummeas,'m5m3catreg','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

#vioplots
#m1-m2
tmpboxplotlabs = c(tmpscatterlabels,tmpscatterlabels,tmpscatterlabels)
fignum="3cm2i" ; xname = xlabvar=tmpboxplotlabs[2] ; yname = ylabvar = tmpscatterlabels[1]
plotvioregression(m2factor,tmpmaxderiv,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol)
#m5-m2
fignum="3cm2ii" ; yname = ylabvar = tmpscatterlabels[5]
plotvioregression(m2factor,tmp50pdecline,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol)

#m1-m3
fignum="3cm3i" ; xname = xlabvar=tmpboxplotlabs[3] ; yname =  ylabvar = tmpscatterlabels[1]
plotvioregression(m3factor,tmpmaxderiv,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol)
#m5-m3
fignum="3cm3ii" ; yname =  ylabvar = tmpscatterlabels[5]
plotvioregression(m3factor,tmp50pdecline,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol)

#get index per ROI if ROI belong to which of the 4 groups
m2w1w2=vector(length=nroi) #m2
m3w1w2=vector(length=nroi) #m3
m2w1w2m3w1w2=vector(length=nroi) #m2 and m3

m2w1w2[which(A=="M2-w1")]=1 
m3w1w2[which(B=="M3-w1")]=1 
m2w1w2[which(A=="M2-w2")]=2
m3w1w2[which(B=="M3-w2")]=2
m2w1w2m3w1w2[which(A=="M2-w1" & B=="M3-w1")]=1
m2w1w2m3w1w2[which(A=="M2-w1" & B=="M3-w2")]=2
m2w1w2m3w1w2[which(A=="M2-w2" & B=="M3-w1")]=3
m2w1w2m3w1w2[which(A=="M2-w2" & B=="M3-w2")]=4

#save for brain plotting
write(m2w1w2, file = paste(outputdirsummeas,'m2w1w2','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector
write(m3w1w2, file = paste(outputdirsummeas,'m3w1w2','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector
write(m2w1w2m3w1w2, file = paste(outputdirsummeas,'m2w1w2m3w1w2','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector

# 2-Way Frequency Table - M2 & M4
A=tmpageatmaxderiv > tmpageatmaxderivval
B=tmpageatsecondzerocrossing > tmpageatsecondzerocrossingval

## b mean and ci vals
bearlyvals = tmpageatsecondzerocrossing[which(B==FALSE)]
blatevals = tmpageatsecondzerocrossing[which(B==TRUE)]
cilev = 0.95 ;
ciearly = CI(bearlyvals, ci = cilev) ;ciearlyup = ciearly[[1]]; ciearlymean = ciearly[[2]]; ciearlylow = ciearly[[3]];
cilate= CI(blatevals, ci = cilev) ; cilateup = cilate[[1]]; cilatemean = cilate[[2]]; cilatelow = cilate[[3]];
ndec=2; write(paste("early wave (n= ",length(bearlyvals),") mean age of ",round(ciearlymean,ndec)," (95% CI = ",round(ciearlylow,ndec),"-",round(ciearlyup,ndec),")","\nlate wave (n= ",length(blatevals),") mean age of ",round(cilatemean,ndec)," (95% CI = ",round(cilatelow,ndec),"-",round(cilateup,ndec),")",sep=""), file = paste(outputdirsummeas,'m4_earlylatewavemeanci','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

A[A==TRUE]="M2-w2"
A[A==FALSE]="M2-w1"
B[B==TRUE]="M4-w2"
B[B==FALSE]="M4-w1"

catmilestones.all[3,]=B

##get mean age and CI M1 and M5 vals of early and late M4 wave
##m1
#m1 vals for m4
tmpmeas=tmpmaxderiv
tmpwavefac = B
tmpdframe = data.frame(tmpmeas,tmpwavefac)
#trad
m1meanci_m4w1w2_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
#boot
m1meanci_m4w1w2_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
write.table(m1meanci_m4w1w2_boot, file = paste(outputdirsummeas,'table1_m1meanci_m4w1w2_boot','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

##m5
#m5 vals for m4
tmpmeas=tmp50pdecline
tmpdframe = data.frame(tmpmeas,tmpwavefac)
#trad
m5meanci_m4w1w2_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
#boot
m5meanci_m4w1w2_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
write.table(m5meanci_m4w1w2_boot, file = paste(outputdirsummeas,'table1_m5meanci_m4w1w2_boot','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

##regression with M1 and M5
#m4-m1
m4factor=factor(B)
glmm1m4cat = lm(tmpmaxderiv ~ Bfactor)
sumglmm1m4cat = summary(glmm1m4cat)
anovaglmm1m4cat = anova(glmm1m4cat)
#m4-m5
glmm5m4cat = lm(tmp50pdecline ~ Bfactor)
sumglmm5m4cat = summary(glmm5m4cat)
anovaglmm5m4cat = anova(glmm5m4cat)
#write results to file
#m1
m1m4reg=paste("Reg coeff = ",glmm1m4cat$coefficients[2],", Rsquared = ",sumglmm1m4cat$r.squared,", F(",anovaglmm1m4cat$Df[1],",",anovaglmm1m4cat$Df[2],") = ",anovaglmm1m4cat$`F value`[1]," p = ",anovaglmm1m4cat$`Pr(>F)`[1],sep="")
#write(m1m4reg, file = paste(outputdirsummeas,'m1m4reg','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
#m5
m5m4reg=paste("Reg coeff = ",glmm5m4cat$coefficients[2],", Rsquared = ",sumglmm5m4cat$r.squared,", F(",anovaglmm5m4cat$Df[1],",",anovaglmm5m4cat$Df[2],") = ",anovaglmm5m4cat$`F value`[1]," p = ",anovaglmm5m4cat$`Pr(>F)`[1],sep="")
#write(m5m4reg, file = paste(outputdirsummeas,'m5m4reg','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector

#vioplots as well
#tmpboxplotlabs = c(tmpscatterlabels,tmpscatterlabels,tmpscatterlabels)
#fignum="3ci" ; xname = xlabvar=tmpboxplotlabs[4] ; yname =  ylabvar = tmpscatterlabels[1]
#plotvioregression(m4factor,tmpmaxderiv,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol)
#fignum="3di" ; yname =  ylabvar = tmpscatterlabels[5]
#plotvioregression(m4factor,tmp50pdecline,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol)

#put all (m2,m3 and m4 corr with m1 and m5) p-values in one vector for multiple comparisons correction
tmpbvals=c(glmm1m2cat$coefficients[2],glmm1m3cat$coefficients[2],glmm5m2cat$coefficients[2],glmm5m3cat$coefficients[2],glmm1m4cat$coefficients[2],glmm5m4cat$coefficients[2])
tmppvals=c(anovaglmm1m2cat$`Pr(>F)`[1],anovaglmm1m3cat$`Pr(>F)`[1],anovaglmm5m2cat$`Pr(>F)`[1],anovaglmm5m3cat$`Pr(>F)`[1],anovaglmm1m4cat$`Pr(>F)`[1],anovaglmm5m4cat$`Pr(>F)`[1])
p = tmppvals ; print(max(p))
fdrtmppvals_bh=p.adjust(p, method = "BH", n = length(p)) ; print(max(fdrtmppvals_bh))
fdrtmppvals_by=p.adjust(p, method = "BY", n = length(p)) ;print(max(fdrtmppvals_by))
fdrtmppvals_bonf=p.adjust(p, method = "bonferroni", n = length(p)) ;print(max(fdrtmppvals_by))
ptable = rbind(tmpbvals,p,fdrtmppvals_by,fdrtmppvals_bonf)
write.table(ptable,file = paste(outputdirsummeas,'corrptable','.alpha',alpha,'.',measname,'.txt', sep = ""),col.names=TRUE,row.names=TRUE)

#get index per ROI if ROI belong to which of the 4 groups
m4w1w2=vector(length=nroi) #m4 (m2saved above)
m4w1w2[which(B=="M4-w1")]=1
m4w1w2[which(B=="M4-w2")]=2
#save for brain plotting
write(m4w1w2, file = paste(outputdirsummeas,'m4w1w2','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector


### CT: Test if steeper CT decline in wave 2 regions (which would be in line with the "last in, first out" hypothesis
#CT-M4
B=tmpageatsecondzerocrossing > tmpageatsecondzerocrossingval
B[B==TRUE]="M4-w2"
B[B==FALSE]="M4-w1"
Bfactor = factor(B)
wave1n="M4-w1" ; wave2n = "M4-w2" ; tmpmeaswCT = tmpageatsecondzerocrossing ; tmpfactorwCT = B
#mean across participants above tage
#set age to use as "cut-off"
tage=round(gageatfirstzerocrossing,1)
#tagebel=15
meanCTroiabtage = colMeans(CT[which(age >= tage),]) 
meanCTroibeltage = colMeans(CT[which(age <= tage),]) 

##reg CT-age with wave
#first get CT-age slopes
CTageslopesabtage = vector(length = dim(CT)[2])
CTageslopesabtageF=CTageslopesabtaget=CTageslopesabtagedir=CTageslopesabtage
myelinageslopesabtage = vector(length = dim(CT)[2])
myelinageslopespabtage = vector(length = dim(CT)[2])
CTageslopesabtagepercofCT = vector(length = dim(CT)[2])

for (cii in 1:dim(CT)[2]) {
  #get participants in a given age range
  ind = which(age>=tage)
  #ct-age slopes above
  tmpregab = lm(CT[ind,cii]~age[ind]) ; sumtmpregab = summary(tmpregab) ; 
  CTageslopesabtage[cii] = tmpregab$coefficients[[2]] #regression slope
  CTageslopesabtageF[cii] = sumtmpregab$fstatistic[[1]] ; #F value
  CTageslopesabtaget[cii] = sumtmpregab$coefficients[2,3] ; #t-value
  CTageslopesabtagepercofCT[cii] = ((CTageslopesabtage[cii]-median(CT[,cii]))/median(CT[,cii]))*100 ;
}

#boxplot CT-age slopes
boxplot(CTageslopesabtage~tmpfactorwCT,notch=TRUE)
#get sign of slope (1 for pos, and -1 for neg)
CTageslopesabtagedir[which(CTageslopesabtage >= 0)] = 1 ; CTageslopesabtagedir[which(CTageslopesabtage < 0)] = -1 ;
#multiply with -1 or 1 for neg and pos slopes respectively, 
CTageslopesabtageFtimesDirSlope=CTageslopesabtageF*CTageslopesabtagedir 
boxplot(CTageslopesabtageFtimesDirSlope~tmpfactorwCT,notch=TRUE) 
boxplot(CTageslopesabtaget~tmpfactorwCT,notch=TRUE) 
#wilcox.test
wtestslope = wilcox.test(CTageslopesabtage~tmpfactorwCT,alternative = c("two.sided")) ; wtestsloperes = paste("W=",wtestslope$statistic[[1]], ", P=",wtestslope$p.value,sep="")
print(wtestsloperes)
wtestswF = wilcox.test(CTageslopesabtageFtimesDirSlope~tmpfactorwCT,alternative = c("two.sided")) ; wtestswFres = paste("W=",wtestswF$statistic[[1]], ", P=",wtestswF$p.value,sep="") ;

#regression - no T1w/T2w as covariate 
tmpregCTagewaveab = lm(CTageslopesabtage~as.factor(tmpfactorwCT))
summary(tmpregCTagewaveab) #sig
##regression -  T1w/T2w as covariate
avgmyelinabtage = colMeans(Y[which(age > tage),])
#F value
tmpregCTagewaveabavgmyelin = lm(CTageslopesabtageFtimesDirSlope~as.factor(tmpfactorwCT)+avgmyelinabtage)
sumtmpregCTagewaveabavgmyelin= summary(tmpregCTagewaveabavgmyelin) #sig
#unstand beta
#tmpregCTagewaveabavgmyelinres = paste("bwave=",sumtmpregCTagewaveabavgmyelin$coefficients[2,1],", P=",sumtmpregCTagewaveabavgmyelin$coefficients[2,4],", adjR=",sumtmpregCTagewaveabavgmyelin$adj.r.squared,sep="")
#stand beta
lmbetatmpregCTagewaveabavgmyelin = lm.beta(tmpregCTagewaveabavgmyelin)
tmpregCTagewaveabavgmyelinres = paste("Bwave=",lmbetatmpregCTagewaveabavgmyelin$standardized.coefficients[[2]],", P=",sumtmpregCTagewaveabavgmyelin$coefficients[2,4],", adjR=",sumtmpregCTagewaveabavgmyelin$adj.r.squared,sep="")

#write result to file
write(wtestsloperes,file = paste(outputdirsummeas,'wilcoxtestslope','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector
write(wtestswFres,file = paste(outputdirsummeas,'wilcoxtestswF','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector
write(tmpregCTagewaveabavgmyelinres,file = paste(outputdirsummeas,'regressionageCTwavemyelin','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector


#Two-by-two plots contingency table
# 2-Way Frequency Table - M3 & M4

A=tmpageatfirstzerocrossing > tmpageatfirstzerocrossingval
B=tmpageatsecondzerocrossing > tmpageatsecondzerocrossingval
A[A==TRUE]="M3-w2"
A[A==FALSE]="M3-w1"
B[B==TRUE]="M4-w2"
B[B==FALSE]="M4-w1"

#get index per ROI if ROI belong to which of the 4 groups
m3w1w2=vector(length=nroi) #m3 (m2 and m4 saved above)

m3w1w2[which(A=="M3-w1")]=1
m3w1w2[which(A=="M3-w2")]=2

#save for brain plotting
write(m3w1w2, file = paste(outputdirsummeas,'m3w1w2','.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) #as row vector


### Relate milestones to other measures ###

##Von economo
vonecodir="voneconomo/"
vonecofile=paste(vonecodir,"hcpvoneconomoclasses.txt",sep='')
vonecovec = read.delim(vonecofile,header=FALSE,sep=",")
origvonecovec = as.numeric(vonecovec)

#set von economo regions not from 1-7 to NaN
voneconomoclasesnames=c("PM", "AC1", "AC2", "PSS", "PS","LB","IC")
univonecovec = sort(unique(origvonecovec))

#get vector of classes names for anova
#both hemi
factorvoneconomovec=c(rep("cx",length(origvonecovec))) 
sortedunivonecovec=sort(univonecovec)
tmppallvoneco = tmpkruskalpallvoneco = tmpfallvoneco = tmppallvonecochisq = tmpchiallvonecochisq = vector()
for (ni in 1:length(sortedunivonecovec)) {
  tmpvecind = origvonecovec == sortedunivonecovec[ni]
  factorvoneconomovec[tmpvecind]=voneconomoclasesnames[ni]
}
#right
factorvoneconomovec_right=c(rep("cx",length(origvonecovec)/2)) 
for (ni in 1:length(sortedunivonecovec)) {
  tmpvecind = origvonecovec[1:180] == sortedunivonecovec[ni]
  factorvoneconomovec_right[tmpvecind]=voneconomoclasesnames[ni]
}
#left
factorvoneconomovec_left=c(rep("cx",length(origvonecovec)/2)) 
for (ni in 1:length(sortedunivonecovec)) {
  tmpvecind = origvonecovec[181:360] == sortedunivonecovec[ni]
  factorvoneconomovec_left[tmpvecind]=voneconomoclasesnames[ni]
}

factorvoneconomo=factor(factorvoneconomovec)
vonecoinds = matrix(nrow = length(univonecovec), ncol = length(t(origvonecovec)))
for (ni in 1:length(univonecovec)) {
  vonecoinds[ni,] = origvonecovec == univonecovec[ni]
}

#for average myelin, and then for each milestone (all, rh, lh), loop over and plot a violin plot
tmpboxplotlabs = c(tmpscatterlabels,tmpscatterlabels,tmpscatterlabels)

hcpind=c(seq(1,nhcprois))
hcpindright=c(seq(1,nhcprois/2))
hcpindleft=c(seq((nhcprois/2)+1,nhcprois))
#bothhemi, right and left
for (mi in 1:(dim(milestones.all)[1]*3)) {
  print(paste("milestone # is ",mi,sep=""))  
  if (mi == 1) {
    #control - myelin level
    for (hi in 1:3) {
      print(hi)
      if (hi == 1) {
        heminame="bothhemi"
        tmpind=hcpind
        tmpfactorvoneconomovec=factorvoneconomovec
      }  else if (hi == 3) {
        heminame="leftthemi"
        tmpind=hcpindleft
        tmpfactorvoneconomovec=factorvoneconomovec_left     
      } else {
        heminame="righthemi"
        tmpind=hcpindright
        tmpfactorvoneconomovec=factorvoneconomovec_right
      }
      meanxrois = colMeans(Y[,tmpind])
      x1 <- meanxrois[vonecoinds[1,tmpind]]
      x2 <- meanxrois[vonecoinds[2,tmpind]]
      x3 <- meanxrois[vonecoinds[3,tmpind]]
      x4 <- meanxrois[vonecoinds[4,tmpind]]
      x5 <- meanxrois[vonecoinds[5,tmpind]]
      x6 <- meanxrois[vonecoinds[6,tmpind]]
      x7 <- meanxrois[vonecoinds[7,tmpind]]
      print(heminame)
      
      #kruskal-wallis
      tmpkruskalwallis = kruskal.test(meanxrois~factor(tmpfactorvoneconomovec))
      tmpoutput2=paste("Chi sq=",round(tmpkruskalwallis$statistic,digits=1),",df=",tmpkruskalwallis$parameter,", P = ",tmpkruskalwallis$p.value,sep="")
      print(tmpoutput2)
      write(tmpoutput2, file = paste(outputdirsummeas,'KruskalWallis_ranksumtest_myelinlevelvoneconomovoneconomo_',heminame,'.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) 
      
    }
  }
}
if (mi < 6) {
  heminame="bothhemi"
  #milestone of interst
  miind=mi
  tmpind=hcpind
  tmpfactorvoneconomovec=factorvoneconomovec
}  else if (mi > 10) {
  miind=mi-10
  heminame="leftthemi"
  #milestone of interst
  tmpind=hcpindleft
  tmpfactorvoneconomovec=factorvoneconomovec_left
} else {
  heminame="righthemi"
  #milestone of interst
  miind=mi-5
  tmpind=hcpindright
  tmpfactorvoneconomovec=factorvoneconomovec_right
}
print(heminame)
#get data based on info above
mmoi = milestones.all[miind,tmpind] #age at max rate, max rate, 1st zero cross, 2nd zero cross, then same for rh and lh
mmoiname = tmpboxplotlabs[miind]
x1 <- mmoi[vonecoinds[1,tmpind]]
x2 <- mmoi[vonecoinds[2,tmpind]]
x3 <- mmoi[vonecoinds[3,tmpind]]
x4 <- mmoi[vonecoinds[4,tmpind]]
x5 <- mmoi[vonecoinds[5,tmpind]]
x6 <- mmoi[vonecoinds[6,tmpind]]
x7 <- mmoi[vonecoinds[7,tmpind]]


#kruskal-wallis
tmpkruskalwallis = kruskal.test(mmoi~factor(tmpfactorvoneconomovec))
tmpoutput2=paste("Chi sq=",round(tmpkruskalwallis$statistic,digits=1),",df=",tmpkruskalwallis$parameter,", P = ",tmpkruskalwallis$p.value,sep="")
tmpkruskalpallvoneco = c(tmpkruskalpallvoneco,tmpkruskalwallis$p.value) #put all kruskal-wallis p vals in onw vector for multiple comparisons correction
print(tmpoutput2)
write(tmpoutput2, file = paste(outputdirsummeas,'KruskalWallis_ranksumtest_voneconomo_',mmoiname,heminame,'.alpha',alpha,'.',measname,'.txt', sep = ""), ncolumns = length(tmpmaxderiv)) 


#multiple comparisons correction
#kruskal
pkruskal = tmpkruskalpallvoneco ; print(max(pkruskal)) ; print(length(which(pkruskal < 0.05)))
fdrtmppvalskruskal_by=p.adjust(pkruskal, method = "BY", n = length(pkruskal)) ;print(max(fdrtmppvalskruskal_by)) ; print(length(which(fdrtmppvalskruskal_by < 0.05)))
kruskalpvalscorr=rbind(pkruskal,fdrtmppvalskruskal_by)
write.table(kruskalpvalscorr,file = paste(outputdirsummeas,'corrptablekruskalvoneconomo','.alpha',alpha,'.',measname,'.txt', sep = ""),col.names=TRUE,row.names=TRUE)


#von economo bar charts for m2, m3, and m4
xlabel="Class - Hemisphere" ; ylabel="# regions"
tmppallvonecochisq = tmpchiallvonecochisq = vector()
for (mi in 1:(dim(milestones.all)[1])) {
  print(mi)
  if (mi == 1 | mi==5) next
  print("mi not age milestone - skipping")
  if (mi == 2) {
    tmpfactor=m2factor
    themearg = theme(legend.position="none") #change, as for m4, if want for all figures
    xlabelarg="" #change to xlabel if want for all figure
    ylabelarg=ylabel
    #print(paste("mi is ",mi,", xlab is ",xlabelarg,", and ylab is ",ylabelarg,sep=""))
  } else if (mi == 3) {
    tmpfactor=m3factor
    themearg = theme(legend.position="none") #change, as for m4, if want for all figures
    xlabelarg=xlabel
    ylabelarg="" #change to ylabel if want for all figure
    #print(paste("mi is ",mi,"xlab is ",xlabelarg,", and ylab is ",ylabelarg,sep=""))
  } else if (mi == 4) {
    tmpfactor=m4factor
    themearg=theme(legend.title=element_blank())
    xlabelarg="" #change to xlabel if want for all figure
    ylabelarg="" #change to ylabel if want for all figure
  }
  mmoiname=tmpboxplotlabs[mi]
  print(mmoiname)
  hemiindvec =rep(c("R","L"),each=length(origvonecovec)/2)
  df = data.frame('class' = factorvoneconomovec,'hemi'=hemiindvec,'wave'=tmpfactor)
  
  dfcast = acast(df,class~hemi~wave,fun.aggregate = length)
  melteddfcast = melt(dfcast)
  colnames(melteddfcast) = c("class","hemi","wave","value")
  #plot
  pdf(paste(outputdirfigs,'supfig11_barchart_voneconomo_',mmoiname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  #need print first to get it to save in script/function
  tmpgp = ggplot() +
    geom_bar(data=melteddfcast, aes(y = value, x = hemi, fill = wave), stat="identity",position='stack') +
    #theme_bw() + 
    theme_minimal() +
    facet_grid( ~ class,switch = "y") +  theme(strip.background = element_blank()) +
    scale_fill_manual(values=c(col1, col2),labels=c("Wave 1", "Wave 2"))
  gp = tmpgp + labs(x = xlabelarg,y=ylabelarg) + themearg
  #gp                  
  print(gp)
  dev.off()
  #publication plot (w/o labels)
  pdf(paste(outputdirfigs,'fig11_barchart_voneconomo_pub_',mmoiname,'.alpha',alpha,'.',measname,'.pdf',sep=''),width=6,height=5)
  #need print first to get it to save in script/function
  tmpgp = ggplot() +
    geom_bar(data=melteddfcast, aes(y = value, x = hemi, fill = wave), stat="identity",position='stack') +
    #theme_bw() + 
    #theme_minimal() +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y =element_blank(),
          axis.text.y=element_blank(),
          strip.text.x=element_blank()) +
    facet_grid( ~ class,switch = "y") +  theme(strip.background = element_blank()) +
    scale_fill_manual(values=c(col1, col2),labels=c("", ""))
  #tmpgp
  #gp = tmpgp + labs(x = xlabel,y=ylabel) + theme(legend.title=element_blank())
  gp = tmpgp +  themearg
  #gp                  
  print(gp)
  dev.off()
  
  #per hemi
  for (hi in 1:2) {
    if (hi==1) {
      tmpind=1:180
      heminame="right"
    } else {
      tmpind=181:360
      heminame="left"
    }
    
    print(heminame)
    tmpcrosstable=CrossTable(factorvoneconomo[tmpind],tmpfactor[tmpind])
    tmpchisqtest = chisq.test(tmpcrosstable$t)
    tmppallvonecochisq = c(tmppallvonecochisq,tmpchisqtest$p.value)
    tmpchiallvonecochisq = c(tmpchiallvonecochisq,tmpchisqtest$statistic[[1]])
    #write results to file
    tmpchisqtestres=paste("chisquare= ",round(tmpchisqtest$statistic,digits = 1),", df=",tmpchisqtest$parameter,", p=",tmpchisqtest$p.value,sep="")
    write(tmpchisqtestres, file = paste(outputdirsummeas,'chisquare_voneconomo_',heminame,mmoiname,'.alpha',alpha,'.',measname,'.txt', sep = "")) 
    
  }
  
}

#############################
### Connectivity - degree ###
#############################

tmppallgraph = tmpballgraph = vector()
tmpboxplotlabs = c(tmpscatterlabels,tmpscatterlabels,tmpscatterlabels)
#meancidataframe for total,intra,inter deg
meancidataframe_totdeg = meancidataframe_partcoefrfmri = matrix(nrow = 6,ncol=3) ; cn = c("mean","lower","upper") ; colnames(meancidataframe_totdeg) = cn ; rn = c("m2w1","m2w2","m3w1","m3w2","m4w1","m4w2") ; rownames(meancidataframe_totdeg) = rn ;
meancidataframe_intradeg = meancidataframe_interdeg = meancidataframe_totdegrfmri = meancidataframe_totdeg ; colnames(meancidataframe_partcoefrfmri) = cn ; rownames(meancidataframe_partcoefrfmri) = rn ; 
#bothhemi, right and left
for (mi in 1:(dim(milestones.all)[1])) {
  
  tmpind=hcpind
  #milestone of interst
  mmoi = milestones.all[mi,] #age at max rate, max rate, 1st zero cross, 2nd zero cross, then same for rh and lh
  mmoiname = tmpboxplotlabs[mi]
  print(mmoiname)
  
  #get graph meas
  #unthresholded
  tmpgmoiname="meanweigtheddegree"
  ifile = paste(struccovdir,tmpgmoiname,sep='') #generated in matlab, using structuralcovariance_X.m
  tmpgmoi = read.delim(ifile,header=FALSE,sep=",")
  
  #tresholded
  cost=0.1
  gamma=1.12 #SC from age range btw M2 and M3
  print(paste("gamma is ",gamma,sep=""))
  
  #degree
  tmpgmoiname="degree"
  ifile = paste(struccovdir,tmpgmoiname,"_cost",cost,sep='') #generated in matlab, using structuralcovariance_btwm2m3.m
  tmpgmoiall = read.delim(ifile,header=FALSE,sep=",") ; vcolnames = names(tmpgmoiall) ; #get col names to use below for rfmri
  tmpgmoi = as.numeric(tmpgmoiall[tmpind])
  
  #add regression with mmoi as dichotomous and then divided by hemisphere
  if (mi == 2 | mi== 3 | mi==4)  {
    print("running splitvioreg 1")
    catmmoi=catmilestones.all[mi-1,tmpind]
    #split vioplot
    xlabvar=mmoiname ; ylabvar=yname
    fignum="4b" ; xname=mmoiname; yname=tmpgmoiname ; 
    plotsplitvioplot(hemiindvec,catmmoi,tmpgmoi,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol) #
    
    #reg - wave + hemi + wave*hemi
    glmtmpgmoimXcat = lm(as.numeric(tmpgmoi) ~ factor(catmmoi) + factor(hemiindvec) + factor(catmmoi)*factor(hemiindvec))
    sumglmtmpgmoimXcat = summary(glmtmpgmoimXcat) ; anovaglmtmpgmoimXcat= anova(glmtmpgmoimXcat)
    #write results to file
    modelpvalue=pf(sumglmtmpgmoimXcat$fstatistic[1], sumglmtmpgmoimXcat$fstatistic[2], sumglmtmpgmoimXcat$fstatistic[3],lower.tail = FALSE)
    #
    wavet_t = sumglmtmpgmoimXcat$coefficients[2,3] ; wave_t_df=sumglmtmpgmoimXcat$fstatistic[[3]] ;
    tmpgmoimXcatreg=paste("Rsquared = ",round(sumglmtmpgmoimXcat$r.squared,digits=2),", model p = ",modelpvalue,", B1 (wave) = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),", b1 (wave) t(df) = ",round(wavet_t,digits=2),"(",round(wave_t_df,digits=3),")",", B1 p = ",sumglmtmpgmoimXcat$coefficients[2,4],", B2 (hemi) = ",round(glmtmpgmoimXcat$coefficients[3],digits=3),", B2 p = ",sumglmtmpgmoimXcat$coefficients[3,4],", B3 (wave*hemi) = ",round(glmtmpgmoimXcat$coefficients[4],digits=3),", B3 p = ",sumglmtmpgmoimXcat$coefficients[4,4],sep="")
    write(tmpgmoimXcatreg, file = paste(outputdirsummeas,"catregwavehemi_",gsub(" ","",xlabvar),ylabvar,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
    tmppallgraph = c(tmppallgraph,sumglmtmpgmoimXcat$coefficients[2,4],sumglmtmpgmoimXcat$coefficients[4,4])
    tmpballgraph = c(tmpballgraph,glmtmpgmoimXcat$coefficients[[2]],glmtmpgmoimXcat$coefficients[[4]])
    
    #reg - both and each hemi separately
    for (hi in 1:3) {
      if (hi == 1) {
        heminame="both"
        hemiind=1:360
      } else if (hi == 2) {
        heminame="right"
        hemiind=1:180
      } else {
        heminame="left"
        hemiind=181:360
      }
      
      #print(paste("hemi is ",heminame,", and hemiind is ",hemiind,sep=""))
      glmtmpgmoimXcat = lm(as.numeric(tmpgmoi[hemiind]) ~ factor(catmmoi[hemiind]))
      sumglmtmpgmoimXcat = summary(glmtmpgmoimXcat) ; anovaglmtmpgmoimXcat= anova(glmtmpgmoimXcat)
      #write results to file
      tmpgmoimXcatreg=paste("Rsquared = ",round(sumglmtmpgmoimXcat$r.squared,digits=2),",  p = ",anovaglmtmpgmoimXcat$`Pr(>F)`[1],", B = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),sep="")
      write(tmpgmoimXcatreg, file = paste(outputdirsummeas,"catreg_",xlabvar,ylabvar,heminame,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
    }
    
    ## get mean and ci gmeas vals for each wave
    tmpmeas=tmpgmoi
    #degree vals for m2,m3,m4
    tmpwavefac = catmmoi
    tmpdframe = data.frame(tmpmeas,tmpwavefac)
    if (mi == 4) {
      #doing man-whitney of stat sig (as boot CI overlap sligthly)
      degdiffm4w1w2_mannwhit=wilcox.test(tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[1])], y = tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[2])],alternative = c("two.sided"))
    }
    #trad
    tmp_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
    #boot
    tmp_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
    if (mi == 2) {
      tind = 1 
    } else if (mi == 3) {
      tind = 3
    } else {
      tind = 5
    }
    
    meancidataframe_totdeg[tind,1] = tmp_boot[1,3] ; meancidataframe_totdeg[tind,2] = tmp_boot[1,6] ; meancidataframe_totdeg[tind,3] = tmp_boot[1,7]
    meancidataframe_totdeg[tind+1,1] = tmp_boot[2,3] ; meancidataframe_totdeg[tind+1,2] = tmp_boot[2,6] ; meancidataframe_totdeg[tind+1,3] = tmp_boot[2,7]
  }
  
  write.table(meancidataframe_totdeg, file = paste(outputdirsummeas,'table1_totaldegree_allwaves','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
  
  #intramodulardegree
  tmpgmoiname="intramodulardegree"
  ifile = paste(struccovdir,tmpgmoiname,"_cost",cost,'_gamma',gamma,sep='') #generated in matlab, using structuralcovariance_btwm2m3.m
  tmpgmoiall = read.delim(ifile,header=FALSE,sep=",")
  tmpgmoi = as.numeric(tmpgmoiall[tmpind])
  
  #add regression with mmoi as dichotomous
  if (mi == 2 | mi==3 | mi==4)  {
    print("running splitvioreg 1")
    catmmoi=catmilestones.all[mi-1,tmpind]
    #split vioplot
    xlabvar=mmoiname ; ylabvar=yname
    fignum="4b" ; xname=mmoiname; yname=tmpgmoiname ; 
    plotsplitvioplot(hemiindvec,catmmoi,tmpgmoi,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol) #
    
    #reg - wave + hemi + wave*hemi
    glmtmpgmoimXcat = lm(as.numeric(tmpgmoi) ~ factor(catmmoi) + factor(hemiindvec) + factor(catmmoi)*factor(hemiindvec))
    sumglmtmpgmoimXcat = summary(glmtmpgmoimXcat) ; anovaglmtmpgmoimXcat= anova(glmtmpgmoimXcat)
    #write results to file
    modelpvalue=pf(sumglmtmpgmoimXcat$fstatistic[1], sumglmtmpgmoimXcat$fstatistic[2], sumglmtmpgmoimXcat$fstatistic[3],lower.tail = FALSE)
    wavet_t = sumglmtmpgmoimXcat$coefficients[2,3] ; wave_t_df=sumglmtmpgmoimXcat$fstatistic[[3]] ;
    tmpgmoimXcatreg=paste("Rsquared = ",round(sumglmtmpgmoimXcat$r.squared,digits=2),", model p = ",modelpvalue,", B1 (wave) = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),", b1 (wave) t(df) = ",round(wavet_t,digits=2),"(",round(wave_t_df,digits=3),")",", B1 p = ",sumglmtmpgmoimXcat$coefficients[2,4],", B2 (hemi) = ",round(glmtmpgmoimXcat$coefficients[3],digits=3),", B2 p = ",sumglmtmpgmoimXcat$coefficients[3,4],", B3 (wave*hemi) = ",round(glmtmpgmoimXcat$coefficients[4],digits=3),", B3 p = ",sumglmtmpgmoimXcat$coefficients[4,4],sep="")
    #tmpgmoimXcatreg=paste("Rsquared = ",round(sumglmtmpgmoimXcat$r.squared,digits=2),",  p = ",anovaglmtmpgmoimXcat$`Pr(>F)`[1],", B1 (wave) = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),", B2 (hemi) = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),", B3 (wave*hemi) = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),sep="")
    write(tmpgmoimXcatreg, file = paste(outputdirsummeas,"catregwavehemi_",gsub(" ","",xlabvar),ylabvar,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
    tmppallgraph = c(tmppallgraph,sumglmtmpgmoimXcat$coefficients[2,4],sumglmtmpgmoimXcat$coefficients[4,4])
    tmpballgraph = c(tmpballgraph,glmtmpgmoimXcat$coefficients[[2]],glmtmpgmoimXcat$coefficients[[4]])
    
    #reg - both and each hemi separately
    for (hi in 1:3) {
      if (hi == 1) {
        heminame="both"
        hemiind=1:360
      } else if (hi == 2) {
        heminame="right"
        hemiind=1:180
      } else {
        heminame="left"
        hemiind=181:360
      }
      
      #print(paste("hemi is ",heminame,", and hemiind is ",hemiind,sep=""))
      glmtmpgmoimXcat = lm(as.numeric(tmpgmoi[hemiind]) ~ factor(catmmoi[hemiind]))
      sumglmtmpgmoimXcat = summary(glmtmpgmoimXcat) ; anovaglmtmpgmoimXcat= anova(glmtmpgmoimXcat)
      #write results to file
      tmpgmoimXcatreg=paste("Rsquared = ",round(sumglmtmpgmoimXcat$r.squared,digits=2),",  p = ",anovaglmtmpgmoimXcat$`Pr(>F)`[1],", B = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),sep="")
      write(tmpgmoimXcatreg, file = paste(outputdirsummeas,"catreg_",xlabvar,ylabvar,heminame,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
    }
    ## get mean and ci gmeas vals for each wave
    tmpmeas=tmpgmoi
    #degree vals for m2,m3,m4
    tmpwavefac = catmmoi
    tmpdframe = data.frame(tmpmeas,tmpwavefac)
    if (mi == 4) {
      #doing man-whitney of stat sig (as boot CI overlap sligthly)
      intradegdiffm4w1w2_mannwhit=wilcox.test(tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[1])], y = tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[2])],alternative = c("two.sided"))
    }
    #trad
    tmp_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
    #boot
    tmp_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
    if (mi == 2) {
      tind = 1 
    } else if (mi == 3) {
      tind = 3
    } else {
      tind = 5
    }
    
    meancidataframe_intradeg[tind,1] = tmp_boot[1,3] ; meancidataframe_intradeg[tind,2] = tmp_boot[1,6] ; meancidataframe_intradeg[tind,3] = tmp_boot[1,7]
    meancidataframe_intradeg[tind+1,1] = tmp_boot[2,3] ; meancidataframe_intradeg[tind+1,2] = tmp_boot[2,6] ; meancidataframe_intradeg[tind+1,3] = tmp_boot[2,7]
    
  }
  
  write.table(meancidataframe_intradeg, file = paste(outputdirsummeas,'table1_intradegree_allwaves','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
  
  ### intermodulardegree ###
  
  tmpgmoiname="intermodulardegree"
  ifile = paste(struccovdir,tmpgmoiname,"_cost",cost,'_gamma',gamma,sep='') #generated in matlab, using structuralcovariance_btwm2m3.m
  tmpgmoiall = read.delim(ifile,header=FALSE,sep=",")
  tmpgmoi = as.numeric(tmpgmoiall[tmpind])

  #add regression with mmoi as dichotomous
  if (mi == 2 | mi==3 | mi==4)  {
    print("running splitvioreg 1")
    catmmoi=catmilestones.all[mi-1,tmpind]
    #split vioplot
    xlabvar=mmoiname ; ylabvar=yname
    fignum="4b" ; xname=mmoiname; yname=tmpgmoiname ; 
    #plotsplitvioplot(catmmoi,hemiindvec,tmpgmoi,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol)
    plotsplitvioplot(hemiindvec,catmmoi,tmpgmoi,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol) #
    
    #reg - wave + hemi + wave*hemi
    glmtmpgmoimXcat = lm(as.numeric(tmpgmoi) ~ factor(catmmoi) + factor(hemiindvec) + factor(catmmoi)*factor(hemiindvec))
    sumglmtmpgmoimXcat = summary(glmtmpgmoimXcat) ; anovaglmtmpgmoimXcat= anova(glmtmpgmoimXcat)
    #write results to file
    modelpvalue=pf(sumglmtmpgmoimXcat$fstatistic[1], sumglmtmpgmoimXcat$fstatistic[2], sumglmtmpgmoimXcat$fstatistic[3],lower.tail = FALSE)
    wavet_t = sumglmtmpgmoimXcat$coefficients[2,3] ; wave_t_df=sumglmtmpgmoimXcat$fstatistic[[3]] ;
    tmpgmoimXcatreg=paste("Rsquared = ",round(sumglmtmpgmoimXcat$r.squared,digits=2),", model p = ",modelpvalue,", B1 (wave) = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),", b1 (wave) t(df) = ",round(wavet_t,digits=2),"(",round(wave_t_df,digits=3),")",", B1 p = ",sumglmtmpgmoimXcat$coefficients[2,4],", B2 (hemi) = ",round(glmtmpgmoimXcat$coefficients[3],digits=3),", B2 p = ",sumglmtmpgmoimXcat$coefficients[3,4],", B3 (wave*hemi) = ",round(glmtmpgmoimXcat$coefficients[4],digits=3),", B3 p = ",sumglmtmpgmoimXcat$coefficients[4,4],sep="")
    write(tmpgmoimXcatreg, file = paste(outputdirsummeas,"catregwavehemi_",gsub(" ","",xlabvar),ylabvar,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
    tmppallgraph = c(tmppallgraph,sumglmtmpgmoimXcat$coefficients[2,4],sumglmtmpgmoimXcat$coefficients[4,4])
    tmpballgraph = c(tmpballgraph,glmtmpgmoimXcat$coefficients[[2]],glmtmpgmoimXcat$coefficients[[4]])
    tmppallgraph_wave = c(tmppallgraph,sumglmtmpgmoimXcat$coefficients[2,4],sumglmtmpgmoimXcat$coefficients[4,4])
    tmpballgraph_wave = c(tmpballgraph,glmtmpgmoimXcat$coefficients[[2]],glmtmpgmoimXcat$coefficients[[4]])
    
    
    #reg - both and each hemi separately
    for (hi in 1:3) {
      if (hi == 1) {
        heminame="both"
        hemiind=1:360
      } else if (hi == 2) {
        heminame="right"
        hemiind=1:180
      } else {
        heminame="left"
        hemiind=181:360
      }
      
      #print(paste("hemi is ",heminame,", and hemiind is ",hemiind,sep=""))
      glmtmpgmoimXcat = lm(as.numeric(tmpgmoi[hemiind]) ~ factor(catmmoi[hemiind]))
      sumglmtmpgmoimXcat = summary(glmtmpgmoimXcat) ; anovaglmtmpgmoimXcat= anova(glmtmpgmoimXcat)
      #write results to file
      tmpgmoimXcatreg=paste("Rsquared = ",round(sumglmtmpgmoimXcat$r.squared,digits=2),",  p = ",anovaglmtmpgmoimXcat$`Pr(>F)`[1],", B = ",round(glmtmpgmoimXcat$coefficients[2],digits=3),sep="")
      write(tmpgmoimXcatreg, file = paste(outputdirsummeas,"catreg_",xlabvar,ylabvar,heminame,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
    }
    ## get mean and ci gmeas vals for each wave
    tmpmeas=tmpgmoi
    #degree vals for m2,m3,m4
    tmpwavefac = catmmoi
    tmpdframe = data.frame(tmpmeas,tmpwavefac)
    if (mi == 4) {
      #doing man-whitney of stat sig (as boot CI overlap sligthly)
      interdegdiffm4w1w2_mannwhit=wilcox.test(tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[1])], y = tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[2])],alternative = c("two.sided"))
    }
    #trad
    tmp_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
    #boot
    tmp_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
    if (mi == 2) {
      tind = 1 
    } else if (mi == 3) {
      tind = 3
    } else {
      tind = 5
    }
    meancidataframe_interdeg[tind,1] = tmp_boot[1,3] ; meancidataframe_interdeg[tind,2] = tmp_boot[1,6] ; meancidataframe_interdeg[tind,3] = tmp_boot[1,7]
    meancidataframe_interdeg[tind+1,1] = tmp_boot[2,3] ; meancidataframe_interdeg[tind+1,2] = tmp_boot[2,6] ; meancidataframe_interdeg[tind+1,3] = tmp_boot[2,7]
  }
  
  write.table(meancidataframe_interdeg, file = paste(outputdirsummeas,'table1_interdegree_allwaves','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
  
  ##########################################
  ### connectivity measures from rs-fMRI ###
  ##########################################
  
  ### degree - rfmri
  #set wavelet wave
  w=2 ;
  tmpgmoiname="degree"
  subdatadir="results/validt1wt2w_204_256_dimregressed/"
  finalrfmri="denoi_filtered_func_data_wds_12mot_csf_fraction_0.5_nearest"
  rfmrimeasdir=paste("projects/lifespant1wt2w/",subdatadir,'rfmri/',finalrfmri,'/avgmats/',sep='')
  ifile = paste(rfmrimeasdir,tmpgmoiname,"_wavelets",w,"_cost",cost,sep='')
  tmpgmoi = read.delim(ifile,header=FALSE,sep=",") 
  tmpgmoi = as.numeric(tmpgmoi[tmpind])
  
  #add regression with mmoi as dichotomous and then divided by hemisphere
  if (mi == 2 | mi== 3 | mi==4)  {
    print("running splitvioreg 1")
    catmmoi=catmilestones.all[mi-1,tmpind]
    #split vioplot
    fignum="5b" ; xname=mmoiname; yname=tmpgmoiname ; 
    xlabvar=mmoiname ; ylabvar=paste(yname,"_rfmri",sep="")
    plotsplitvioplot(hemiindvec,catmmoi,tmpgmoi,fignum,xname,yname,xlabvar,ylabvar,gtmpmodcol,gtmpmodderivcol) #
    
    ## get mean and ci gmeas vals for each wave
    tmpmeas=tmpgmoi
    #degree vals for m2,m3,m4
    tmpwavefac = catmmoi
    tmpdframe = data.frame(tmpmeas,tmpwavefac)
    
    #trad
    tmp_trad=group.CI(tmpmeas ~ tmpwavefac,data=tmpdframe,ci = 0.95)
    #boot
    tmp_boot=groupwiseMean(tmpmeas ~ tmpwavefac,data=tmpdframe,conf= 0.95,digits = 3,R=10000,boot = TRUE,traditional = FALSE,normal= FALSE,basic= FALSE,percentile  = FALSE,bca= TRUE)
    if (mi == 2) {
      tind = 1 
    } else if (mi == 3) {
      tind = 3
    } else {
      tind = 5
    }
    
    meancidataframe_totdegrfmri[tind,1] = tmp_boot[1,3] ; meancidataframe_totdegrfmri[tind,2] = tmp_boot[1,6] ; meancidataframe_totdegrfmri[tind,3] = tmp_boot[1,7]
    meancidataframe_totdegrfmri[tind+1,1] = tmp_boot[2,3] ; meancidataframe_totdegrfmri[tind+1,2] = tmp_boot[2,6] ; meancidataframe_totdegrfmri[tind+1,3] = tmp_boot[2,7]
    
    #do  man-whitney of stat sig (necessary to check difference if boot CI overlap sligthly)
    degdiff_mannwhit=wilcox.test(tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[1])], y = tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[2])],alternative = c("two.sided"))
    write.table(paste("Wilcoxon rank sum test between wave 1 and 2: W = ",degdiff_mannwhit$statistic[[1]],", P = ",degdiff_mannwhit$p.value," (two-sided)",sep=""), file = paste(outputdirsummeas,'degdiff_mannwhit_',xlabvar,'_',ylabvar,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
    
  }
  
  write.table(meancidataframe_totdegrfmri, file = paste(outputdirsummeas,'table1_totaldegree_rfmri_wavelets',w,'_allwaves','.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
  
  #do  man-whitney of stat sig (necessary to check difference if boot CI overlap sligthly)
  degdiff_mannwhit=wilcox.test(tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[1])], y = tmpmeas[which(tmpwavefac == levels(factor(tmpwavefac))[2])],alternative = c("two.sided"))
  write.table(paste("Wilcoxon rank sum test between wave 1 and 2: W = ",degdiff_mannwhit$statistic[[1]],", P = ",degdiff_mannwhit$p.value," (two-sided)",sep=""), file = paste(outputdirsummeas,'degdiff_mannwhit_',xlabvar,'_',ylabvar,'.alpha',alpha,'.',measname,'.txt', sep = "")) #as row vector
  
}

#multiple comparisons correction
p = tmppallgraph ; print(max(p)) ; print(length(which(p < 0.05)))
fdrtmppvals_bh=p.adjust(p, method = "BH", n = length(p)) ; print(max(fdrtmppvals_bh)) ; print(length(which(fdrtmppvals_bh < 0.05)))
fdrtmppvals_by=p.adjust(p, method = "BY", n = length(p)) ;print(max(fdrtmppvals_by)) ; print(length(which(fdrtmppvals_by < 0.05)))
degpvalscorr=rbind(tmpballgraph,p,fdrtmppvals_by)
write.table(degpvalscorr,file = paste(outputdirsummeas,'corrptabledegree','.alpha',alpha,'.',measname,'.txt', sep = ""),col.names=TRUE,row.names=TRUE)


