# Program to make figures that appear in figure 2 of manuscript.  The top panels are Timeseries of biomass for collapsed and non-collapse stocks.  The bottom panels show population productivity and exploitation rate for collapsed and non-collapsed stocks



workdir<- # set workdir to path that contains allforagedata.csv
setwd(workdir)

# load data
datafilename<-"allforagedata.csv"

setwd(workdir)
data<-read.csv(file=datafilename,header=TRUE)
###### Initialize Data ###################################
stockid<-unique(data$ASSESSID)
n.stocks<-length(stockid)
min.year<-min(data$Year)
max.year<-max(data$Year)
n.years<-max.year-min.year+1
year.list<-seq(min.year,max.year,by=1)
outputB<-matrix(NA,nrow=n.years,ncol=n.stocks)
# go through each stock, get TB/mean(TB) by year and place in matrix
biomass<-tapply(data$TB,data$ASSESSID,mean,na.rm=TRUE)
min.biomass<-rep(NA,n.stocks)
max.biomass<-rep(NA,n.stocks)
first.year.listB<-rep(NA,n.stocks)
last.year.listB<-rep(NA,n.stocks)
outputSP<-matrix(NA,nrow=n.years,ncol=n.stocks)
outputC<-matrix(NA,nrow=n.years,ncol=n.stocks)
# function to calculate surplus production
calc.SP<-function(X) X[-1,1]-X[-nrow(X),1]+X[-nrow(X),2]

for (i in 1:n.stocks){
  stockindex<-which(data$ASSESSID==stockid[i])
  TB.tmp<-data$TB[stockindex] 
  Year.tmp<-data$Year[stockindex]
  TC.tmp<-data$TC[stockindex] # total catch
  TL.tmp<-data$TL[stockindex] # total landings
  # use TC if available, otherwise use TL
  TC<-TL.tmp
  if(any(!is.na(TC.tmp))) {
    TC<-TC.tmp
  }
  n.TB<-length(which(!is.na(TB.tmp)))
  first.year.TB<-Year.tmp[which(!is.na(TB.tmp))[1]]
  last.year<-max(Year.tmp)
  last.year.listB[i]<-last.year
  first.year<-min(Year.tmp)
  ly.index<-which(year.list==last.year)
  fy.index<-which(year.list==first.year)
  if (n.TB==0) {
    B.list<-rep(NA,last.year-first.year+1)
    C.list<-rep(NA,last.year-first.year+1)
  } else {
    B.list<-TB.tmp/mean(TB.tmp,na.rm=TRUE)
    C.list<-TC/mean(TB.tmp,na.rm=TRUE) # catch standardized by mean biomass
    first.year.listB[i]<-first.year.TB
  }
  outputB[fy.index:ly.index,i]<-B.list
  outputC[fy.index:ly.index,i]<-C.list
  outputSP[fy.index:(ly.index-1),i]<-calc.SP(cbind(B.list,C.list))
} 
colnames(outputB)<-stockid
rownames(outputB)<-seq(min.year,max.year)
## Repeat for Exploitation Rate
outputF<-matrix(NA,nrow=n.years,ncol=n.stocks)
first.year.listF<-rep(NA,n.stocks)
last.year.listF<-rep(NA,n.stocks)
# turn F into an approximate Exploitation Fraction 1-exp(-F)
data$F<-rep(1,nrow(data))-exp(-data$F)
for (i in 1:n.stocks){
  stockindex<-which(data$ASSESSID==stockid[i])
  F.tmp<-data$F[stockindex] # Fishing mortality rate
  ER.tmp<-data$ER[stockindex] # exploitation rate catch / biomass
  
  Year.tmp<-data$Year[stockindex]
  # figure out whether to use F.tmp or ER.tmp (use ER unless unavailable)
  n.F<-length(which(!is.na(F.tmp)))
  n.ER<-length(which(!is.na(ER.tmp)))
  
  first.year.F<-Year.tmp[which(!is.na(F.tmp))[1]]
  first.year.ER<-Year.tmp[which(!is.na(ER.tmp))[1]]
  first.year<-min(Year.tmp)
  last.year<-max(Year.tmp)
  if(n.ER>1){
    # only use F if ER is not available
    F.list<-ER.tmp
    first.year.listF[i]<-first.year.ER
  } else {
    F.list<-F.tmp
    first.year.listF[i]<-first.year.F
  }
  ly.index<-which(year.list==last.year)
  fy.index<-which(year.list==first.year)
  outputF[fy.index:ly.index,i]<-F.list
  last.year.listF[i]<-last.year
} 
colnames(outputF)<-stockid
rownames(outputF)<-seq(min.year,max.year)


################## done initializing data ###########################

library(RColorBrewer)
library(KernSmooth)

# Plotting Functions Here

# function to interpolate F along a segment to make the plot more smooth
interp.F.col<-function(first.F,last.F,n.interps,F.index.list,F.colors,first.B,last.B){
  interp.F<-approx(c(0,1),c(first.F,last.F),xout=seq(0,1,length.out=n.interps))$y
  max.F<-max(F.index.list)
  interp.F<-replace(interp.F,which(interp.F>=max.F),max.F)
  # round the interp.F to nearest 0.01
  interp.F<-round(100*interp.F)/100
  colors.2.use<-F.colors[match(interp.F,F.index.list)]
  interp.B<-approx(c(0,1),c(first.B,last.B),xout=seq(0,1,length.out=n.interps))$y
  time.list<-seq(0,1,length.out=n.interps)
  return(list(interp.B,colors.2.use,time.list))
}

## Function to generate percentiles of kernel density in each year
plot.kernel.ts<-function(fitted.states,percentiles,kernel.range){
  short.year.list<--10:10
  n.stocks=nrow(fitted.states)
  n.years<-length(short.year.list)
  
  kernel.output<-matrix(NA,nrow=n.years,ncol=length(percentiles))
  
  for (i in 1:n.years){
    year.data<-fitted.states[!is.na(fitted.states[,i]),i]
    
    # use bkde kernel density smoother, which automatically fits bandwidth based on variance in x
    smoothed<-bkde(year.data,kernel="normal",range.x=kernel.range)
    deltax<-smoothed$x[2]-smoothed$x[1]
    prob<-smoothed$y*deltax # turn density into probabilities
    # adjust so all density is betweeen 0 and 1
    prob<-prob/sum(prob)
    cumprob<-rep(NA,length(prob))
    for (j in 1:length(prob)){
      cumprob[j]<-sum(prob[1:j])
    }
    
    kernel.output[i,]=approx(cumprob,smoothed$x,percentiles)$y
  }
  return(kernel.output)
}

# add transparency to a list of colors
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


###########################
std.inner.mar<-c(.7,0.35,0.25,0.3)
par(las=1,mai=std.inner.mar,omi=c(1.5,2.0,1.5,2.0))
nf<-layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=TRUE),widths=c(4,1.5,4,1.5))
txt.mult<-1.25

collapse.threshold<-0.25 
collapse.count<-0
# stocks that have at least 25 years of both B and F data
count.years.stocks<-function(x) length(which(!is.na(x)))

n.years.B<-apply(outputB,2,count.years.stocks)
n.years.F<-apply(outputF,2,count.years.stocks)

#extract stocks with at least 25 years of data
B.index<-which(n.years.B>=25)
F.index<-which(n.years.F>=25)
stock.index<-intersect(B.index,F.index)
outputB<-outputB[,stock.index]
outputF<-outputF[,stock.index]
outputSP<-outputSP[,stock.index]
years<-rownames(outputB)
stock.collapse<-c()
year.collapse<-c()
prepost.data<-vector("list",ncol(outputB))

for (i in 1:length(stock.index)){
  B.tmp<-floor(100*outputB[,i])/100 # round down to nearest 0.01
  F.tmp<-outputF[,i]
  min.B<-min(B.tmp,na.rm=TRUE) # round down to nearest 0.01
  
  if (min.B<=collapse.threshold){
    collapse.count=collapse.count+1
    year.index<-which(B.tmp<=collapse.threshold)[1]
    first.notNA.index<-which(!is.na(B.tmp))[1]
    
    year.collapse<-c(year.collapse,year.index)
    # only save output if there is at least three years of pre collapse data
    if (year.index>=(first.notNA.index+3)){
      stock.collapse<-c(stock.collapse,i)
      min.year.index<-max(1,(year.index-10)) # flags so that index is never less than 1
      max.year.index<-min(year.index+10,length(years)) # and index does not exceed length of data
      time.list<-seq(min.year.index,max.year.index)-rep(year.index,length(seq(min.year.index,max.year.index)))
      if (max.year.index>nrow(outputSP)){
        # extract surplus production for these years
        tmp.SP<-c(outputSP[min.year.index:(max.year.index-1),i],NA)
        
      } else {
        tmp.SP<-c(outputSP[min.year.index:(max.year.index),i])
      }
      # save all information on this stock into prepost.data
      data.2.use<-cbind(time.list,outputB[min.year.index:max.year.index,i], outputF[min.year.index:max.year.index,i],tmp.SP)
      prepost.data[[i]]<-data.2.use
      
    }
  }
}

# rest of commands are for plotting

#  plot B time series as a function of years since collapse, with color proportional to F
color.palette<-colorRampPalette(rev(brewer.pal(10,"RdYlBu")))
F.index.list<-seq(0,0.91,by=0.01)
F.index.list=round(100*F.index.list)/100
F.colors<-color.palette(length(F.index.list))
# setup first plot 
# plot the collapses
plot(0,0,xlim=c(-10,10),ylim=c(0,3),type="n",main="",ylab="",xlab="Year Since Collapse",cex.axis=txt.mult,cex.lab=txt.mult,xaxs="i")
par(las=0)
mtext(text="Standardized",side=2,outer=FALSE,line=4.5)
mtext(text="Biomass",side=2,outer=FALSE,line=3)
par(las=1)
# loop through collapsed stocks and plot each
for (i in 1:length(stock.collapse)){
  stock.index.plot<-stock.collapse[i]
  stock.data<-prepost.data[[stock.index.plot]]
  #Loop through years
  for (t in 2:nrow(stock.data)){
    if (any(is.na(stock.data[(t-1):t,3]))) {
      col.2.use="black"
      lty="dotted"
      # do this if there is no F data in year interval
      segments(stock.data[t-1,1],stock.data[t-1,2],stock.data[t,1],stock.data[t,2],col="black",lty="dotted",lwd=2)
    } else {  
      # otherwise code line color based on F
      interp.plot<-interp.F.col(stock.data[t-1,3],stock.data[t,3],100,F.index.list,F.colors,stock.data[t-1,2],stock.data[t,2])
      B.plot<-interp.plot[[1]]
      col.plot<-interp.plot[[2]]
      time.plot<-rep(stock.data[t-1,1],length(B.plot))+interp.plot[[3]]
      for (j in 2:100){
        #turn each year into 200 segments where F is interpolated
        segments(time.plot[j-1],B.plot[j-1],time.plot[j],B.plot[j],col=col.plot[j],lwd=2,lty="solid")
      }  
    }
    
  }  
}

# Make colormap of F levels
par(mai=c(0.6,0,.5,.3),xpd=TRUE)
plot(0,0,type="n",xlim=c(0,5),ylim=c(0,.9),axes=FALSE,ylab="",xlab="")
# loop through colors, make squares
xmin<-0.5
xmax<-1.5
for (i in 2:length(F.index.list)){
  F.color<-F.colors[i-1]
  x<-c(xmin,xmax,xmax,xmin)
  y<-c(F.index.list[i-1],F.index.list[i-1],F.index.list[i],F.index.list[i])
  polygon(x,y,col=F.color,border=F.color)
}
F.plot.list<-c(0,0.2,0.4,0.6,0.8)
par(xpd=TRUE)
text(x=rep(1.5,length(F.plot.list)),y=F.plot.list,labels=F.plot.list,pos=4,cex=txt.mult)
text(x=0.0,y=1.1,labels="Fishing",pos=4,cex=txt.mult)
text(x=0.0,y=1.0,labels="Rate",pos=4,cex=txt.mult)



# Plot the F distributions as a function of year since collapse
# setup matrix to hold results
year.f.data<-matrix(NA,nrow=length(stock.collapse),ncol=21)
for (i in 1:length(stock.collapse)){
  stock.f.data<-prepost.data[[stock.collapse[[i]]]][,3]
  for (t in 1:length(stock.f.data)){
    year.f.data[i,t]<-stock.f.data[t]
  } 
}

col.2.plot<-add.alpha(F.colors[length(F.colors)],0.5)
# set up plot 
span=0.2
degree=2
n.points=41
family="gaussian"
ylims<-c(-.15,0.55)
xlims<-c(-10,10)
short.year.list<--10:10
# just set up plot location, nothing is in here
par(las=1,mai=std.inner.mar)
f.collapse.mean<-apply(year.f.data,2,mean,na.rm=T)
#apply a loess smoother to the mean
f.smooth=loess.smooth(-10:10, f.collapse.mean, span = span, degree = degree,
                      family = family, evaluation = n.points)
x<-f.smooth$x
y<-f.smooth$y
plot(x,y,type="n",lwd=2,main="",xlab="Year Since Collapse",ylab="Exploitation Rate",xlim=xlims,ylim=ylims,xaxs="i",yaxs="i",axes=FALSE,cex.axis=txt.mult,cex.lab=txt.mult)
box()
axis(1,at=seq(-10,10,by=5),cex.axis=txt.mult)
y.at.lab<-seq(-.0,.5,by=0.25)
y.at<-(y.at.lab)
axis(side=2,at=y.at,labels=y.at.lab,cex.axis=txt.mult)
# plot a +/- standard error polygon
f.collapse.SE<-apply(year.f.data,2,sd,na.rm=T)/sqrt(apply(year.f.data,2,function(x) length(which(!is.na(x)))))

bottom.poly<-f.collapse.mean-f.collapse.SE
top.poly<-f.collapse.mean+f.collapse.SE

bottom.poly.smooth<-loess.smooth(-10:10, bottom.poly, span = span, degree = degree,
                                 family = family, evaluation = n.points*2)
top.poly.smooth<-loess.smooth(-10:10, top.poly, span = span, degree = degree,
                              family = family, evaluation = n.points*2)

apply(year.f.data,2,function(x) length(which(!is.na(x))))
x.poly<-c(bottom.poly.smooth$x,rev(top.poly.smooth$x))
y.poly<-c(bottom.poly.smooth$y,rev(top.poly.smooth$y))
polygon(x.poly,y.poly,col=col.2.plot,border=col.2.plot)
lines(x,y,lwd=2,col=F.colors[length(F.colors)])


### Make plot of Surplus Production
year.SP.data<-matrix(NA,nrow=length(stock.collapse),ncol=21)
for (i in 1:length(stock.collapse)){
  stock.SP.data<-prepost.data[[stock.collapse[[i]]]][,4]
  for (t in 1:length(stock.SP.data)){
    year.SP.data[i,t]<-stock.SP.data[t]
  } 
}

# to improve plotting clarity, remove outliers (highest and lowest observation in each year)
meansd.SP.minmax<-function(x){
  min.x.index<-which(x==min(x,na.rm=T))
  max.x.index<-which(x==max(x,na.rm=T))
  new.x<-x[-c(min.x.index,max.x.index)]
  mean.x<-mean(new.x,na.rm=T)
  mean.se<-sd(new.x,na.rm=T)/sqrt(length(which(!is.na(new.x))))
  return(c(mean.x,mean.se))
  
  return(c(mean.x,mean.se))
  
}

SP.collapse.total<-apply(year.SP.data,2,meansd.SP.minmax)
SP.collapse.mean<-SP.collapse.total[1,]
SP.collapse.SE<-SP.collapse.total[2,]

col.2.plot<-add.alpha(F.colors[1],0.5)

# set up plot 
par(las=1,mai=std.inner.mar)

# plot a bunch of polygons
# apply loess smooth to mean surplus production
SP.collapse.mean.smooth<-loess.smooth(-10:10, SP.collapse.mean, span = span, degree = degree,
                                      family = family, evaluation = n.points)
# define bounds of polygone describing +/- SE
bottom.poly<-SP.collapse.mean-SP.collapse.SE
top.poly<-SP.collapse.mean+SP.collapse.SE

bottom.poly.smooth<-loess.smooth(-10:10, bottom.poly, span = span, degree = degree,
                                 family = family, evaluation = n.points*2)
top.poly.smooth<-loess.smooth(-10:10, top.poly, span = span, degree = degree,
                              family = family, evaluation = n.points*2)

x.poly<-c(bottom.poly.smooth$x,rev(top.poly.smooth$x))
y.poly<-c(bottom.poly.smooth$y,rev(top.poly.smooth$y))
par(xpd=F)
polygon(x.poly,y.poly,col=col.2.plot,border=col.2.plot)

lines(SP.collapse.mean.smooth$x,SP.collapse.mean.smooth$y,lwd=2,col=F.colors[1])
label.txt<-expression(paste("or ", tilde(lambda)))
par(las=0)
mtext(text="Fishing Rate",side=2,line=4.5,col=F.colors[length(F.colors)])
mtext(text="Production Rate",side=2,outer=FALSE,line=3,col=F.colors[1])
par(las=1)

##### repeat Analysis for non-collapsed stocks  ################################
std.inner.mar<-c(.7,0.35,0.25,0.3)
par(las=1,mai=std.inner.mar,omi=c(1.5,2.0,1.5,2.0))
nf<-layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=TRUE),widths=c(4,1.5,4,1.5))
txt.mult<-1.25

collapse.threshold<-0.3
collapse.count<-0
# stocks that have at least 25 years of both B and F data

count.years.stocks<-function(x) length(which(!is.na(x)))

n.years.B<-apply(outputB,2,count.years.stocks)
n.years.F<-apply(outputF,2,count.years.stocks)

B.index<-which(n.years.B>=25)
F.index<-which(n.years.F>=25)
stock.index<-intersect(B.index,F.index)


B.index<-which(n.years.B>=25)
F.index<-which(n.years.F>=25)
stock.index<-intersect(B.index,F.index)
years<-rownames(outputB)
stock.nocollapse<-c()
stock.collapse<-c()
year.collapse<-c()
prepost.data<-vector("list",ncol(outputB))
no.collapse.count<-0
for (i in 1:length(stock.index)){
  stock.2.use<-stock.index[i]
  B.tmp<-outputB[,stock.2.use]
  F.tmp<-outputF[,stock.2.use]
  min.B<-min(B.tmp,na.rm=TRUE)
  
  if (min.B>=collapse.threshold){
    no.collapse.count=no.collapse.count+1   
    year.index<-which(B.tmp==min.B)[1]
    first.notNA.index<-which(!is.na(B.tmp))[1]
    year.collapse<-c(year.collapse,years[which(B.tmp==min.B)])
    # only save output if there is at least four years of pre collapse data
    if (year.index>=(first.notNA.index+4)){
      stock.nocollapse<-c(stock.nocollapse,stock.2.use)
      min.year.index<-max(1,(year.index-10)) # flags so that index is never less than 1
      max.year.index<-min(year.index+10,length(years)) # and index does not exceed length of data
      time.list<-seq(min.year.index,max.year.index)-rep(year.index,length(seq(min.year.index,max.year.index)))
      if (max.year.index>nrow(outputSP)){
        tmp.SP<-c(outputSP[min.year.index:(max.year.index-1),stock.2.use],NA)
        
      } else {
        tmp.SP<-c(outputSP[min.year.index:(max.year.index),stock.2.use])
      }
      data.2.use<-cbind(time.list,outputB[min.year.index:max.year.index,stock.2.use], outputF[min.year.index:max.year.index,stock.2.use],tmp.SP)
      prepost.data[[stock.2.use]]<-data.2.use
      
    }
  }
}


# now plot B time series as a function of years since collapse, with color proportional to F
color.palette<-colorRampPalette(rev(brewer.pal(10,"RdYlBu")))
F.index.list<-seq(0,0.91,by=0.01)
F.index.list=round(100*F.index.list)/100
F.colors<-color.palette(length(F.index.list))
# setup first plot

# plot the non-collapses
plot(0,0,xlim=c(-10,10),ylim=c(0,3),type="n",main="",ylab="",xlab="Year Since Minimum Biomass",cex.axis=txt.mult,cex.lab=txt.mult,xaxs="i")
par(las=0)
mtext(text="Standardized",side=2,outer=FALSE,line=4.5)
mtext(text="Biomass",side=2,outer=FALSE,line=3)
par(las=1)
for (i in 1:length(stock.nocollapse)){
  stock.index.plot<-stock.nocollapse[i]
  stock.data<-prepost.data[[stock.index.plot]]
  #Loop through years
  for (t in 2:nrow(stock.data)){
    if (any(is.na(stock.data[(t-1):t,3]))) {
      col.2.use="black"
      lty="dotted"
      segments(stock.data[t-1,1],stock.data[t-1,2],stock.data[t,1],stock.data[t,2],col="black",lty="dotted",lwd=2)
    } else {  
      interp.plot<-interp.F.col(stock.data[t-1,3],stock.data[t,3],100,F.index.list,F.colors,stock.data[t-1,2],stock.data[t,2])
      B.plot<-interp.plot[[1]]
      col.plot<-interp.plot[[2]]
      time.plot<-rep(stock.data[t-1,1],length(B.plot))+interp.plot[[3]]
      for (j in 2:100){
        segments(time.plot[j-1],B.plot[j-1],time.plot[j],B.plot[j],col=col.plot[j],lwd=2,lty="solid")
      }  
    }
    
  }  
}

# Make colormap of F levels
par(mai=c(0.6,0,.5,.3),xpd=TRUE)
plot(0,0,type="n",xlim=c(0,5),ylim=c(0,.9),axes=FALSE,ylab="",xlab="")
# loop through colors, make squares
xmin<-0.5
xmax<-1.5
for (i in 2:length(F.index.list)){
  F.color<-F.colors[i-1]
  x<-c(xmin,xmax,xmax,xmin)
  y<-c(F.index.list[i-1],F.index.list[i-1],F.index.list[i],F.index.list[i])
  polygon(x,y,col=F.color,border=F.color)
}
F.plot.list<-c(0,0.2,0.4,0.6,0.8)
par(xpd=TRUE)
text(x=rep(1.5,length(F.plot.list)),y=F.plot.list,labels=F.plot.list,pos=4,cex=txt.mult)
text(x=0.0,y=1.1,labels="Fishing",pos=4,cex=txt.mult)
text(x=0.0,y=1.0,labels="Rate",pos=4,cex=txt.mult)



# Plot the f distributions as a function of year since minimum biomass
# setup matrix to hold results
year.f.data<-matrix(NA,nrow=length(stock.nocollapse),ncol=21)
for (i in 1:length(stock.nocollapse)){
  stock.f.data<-prepost.data[[stock.nocollapse[[i]]]][,3]
  for (t in 1:length(stock.f.data)){
    year.f.data[i,t]<-stock.f.data[t]
  } 
}

# set up plot 
span=0.2
degree=2
n.points=41
family="gaussian"
ylims<-c(-.15,0.55)
xlims<-c(-10,10)
short.year.list<--10:10
col.2.plot<-add.alpha(F.colors[length(F.colors)],0.5)
# just set up plot location, nothing is in here
par(las=1,mai=std.inner.mar)
f.nocollapse.mean<-apply(year.f.data,2,mean,na.rm=T)
f.smooth=loess.smooth(-10:10, f.nocollapse.mean, span = span, degree = degree,
                      family = family, evaluation = n.points)
x<-f.smooth$x
y<-f.smooth$y
plot(x,y,type="n",lwd=2,main="",xlab="Year Since Minimum Biomass",ylab="Exploitation Rate",xlim=xlims,ylim=ylims,xaxs="i",yaxs="i",axes=FALSE,cex.axis=txt.mult,cex.lab=txt.mult)
box()
axis(1,at=seq(-10,10,by=5),cex.axis=txt.mult)
y.at.lab<-seq(-.0,.5,by=0.25)
y.at<-(y.at.lab)
axis(side=2,at=y.at,labels=y.at.lab,cex.axis=txt.mult)
# plot a 50% range using polygons
f.nocollapse.SE<-apply(year.f.data,2,sd,na.rm=T)/sqrt(apply(year.f.data,2,function(x) length(which(!is.na(x)))))

bottom.poly<-f.nocollapse.mean-f.nocollapse.SE
top.poly<-f.nocollapse.mean+f.nocollapse.SE

bottom.poly.smooth<-loess.smooth(-10:10, bottom.poly, span = span, degree = degree,
                                 family = family, evaluation = n.points*2)
top.poly.smooth<-loess.smooth(-10:10, top.poly, span = span, degree = degree,
                              family = family, evaluation = n.points*2)

apply(year.f.data,2,function(x) length(which(!is.na(x))))
x.poly<-c(bottom.poly.smooth$x,rev(top.poly.smooth$x))
y.poly<-c(bottom.poly.smooth$y,rev(top.poly.smooth$y))
polygon(x.poly,y.poly,col=col.2.plot,border=col.2.plot)
lines(x,y,lwd=2,col=F.colors[length(F.colors)])

### Make plot of Surplus Production
year.SP.data<-matrix(NA,nrow=length(stock.nocollapse),ncol=21)
for (i in 1:length(stock.nocollapse)){
  stock.SP.data<-prepost.data[[stock.nocollapse[[i]]]][,4]
  for (t in 1:length(stock.SP.data)){
    year.SP.data[i,t]<-stock.SP.data[t]
  } 
}

# to improve plotting clarity, remove outliers (highest and lowest observation in each year)
meansd.SP.minmax<-function(x){
  min.x.index<-which(x==min(x,na.rm=T))
  max.x.index<-which(x==max(x,na.rm=T))
  new.x<-x[-c(min.x.index,max.x.index)]
  mean.x<-mean(new.x,na.rm=T)
  mean.se<-sd(new.x,na.rm=T)/sqrt(length(which(!is.na(new.x))))
  return(c(mean.x,mean.se))
  
  return(c(mean.x,mean.se))
  
}

SP.nocollapse.total<-apply(year.SP.data,2,meansd.SP.minmax)
SP.nocollapse.mean<-SP.nocollapse.total[1,]
SP.nocollapse.SE<-SP.nocollapse.total[2,]

col.2.plot<-add.alpha(F.colors[1],0.5)

# set up plot 
par(las=1,mai=std.inner.mar)

# plot a bunch of polygons

SP.nocollapse.mean.smooth<-loess.smooth(-10:10, SP.nocollapse.mean, span = span, degree = degree,
                                     family = family, evaluation = n.points)

bottom.poly<-SP.nocollapse.mean-SP.nocollapse.SE
top.poly<-SP.nocollapse.mean+SP.nocollapse.SE

bottom.poly.smooth<-loess.smooth(-10:10, bottom.poly, span = span, degree = degree,
                                 family = family, evaluation = n.points*2)
top.poly.smooth<-loess.smooth(-10:10, top.poly, span = span, degree = degree,
                              family = family, evaluation = n.points*2)

x.poly<-c(bottom.poly.smooth$x,rev(top.poly.smooth$x))
y.poly<-c(bottom.poly.smooth$y,rev(top.poly.smooth$y))
par(xpd=F)
polygon(x.poly,y.poly,col=col.2.plot,border=col.2.plot)

lines(SP.nocollapse.mean.smooth$x,SP.nocollapse.mean.smooth$y,lwd=2,col=F.colors[1])
par(las=0)
mtext(text="Fishing Rate",side=2,line=4.5,col=F.colors[length(F.colors)])
mtext(text="Production Rate",side=2,outer=FALSE,line=3,col=F.colors[1])
par(las=1)

# End of code