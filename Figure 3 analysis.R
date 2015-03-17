# Randomization / Monte Carlo test to get expected proportions of stocks above and below
# threshold values if stocks are non stationary, constant variance, and have 1/fb variance scaling

# load packages
library(astsa)
library(KernSmooth)
wordir<- # set path to directory containing allforagedata.csv

################################################
# Functions to load

rkernel<-function(p,est.kernel){
  # function to return random draw from a kernel density function
  # generate cumulative probability density
  x<-est.kernel$x
  density<-est.kernel$y
  deltax<-x[2]-x[1]
  cumprob<-rep(NA,length(x))
  cumprob[1]<-0
  for (i in 2:length(x)){
    cumprob[i]<-sum(density[1:i]*deltax)
  }
  # make sure sums to 1
  cumprob<-cumprob/max(cumprob)
  
  output.x<-rep(NA,length(p))
  for (i in 1:length(p)){
    
    over.index<-which(cumprob>=p[i])[1]
    test.x<-x[(over.index-1):over.index]
    test.y<-cumprob[(over.index-1):over.index]
    output.x[i]<-approx(test.y,test.x,xout=p[i])$y
  }
  return(output.x)
}

rkernel2D<-function(N,kernel.est){
  # simulate draws from 2D kernel density estimate, returns N draws of x1 and x2
  x1<-kernel.est$x1
  x2<-kernel.est$x2
  fhat<-sd.beta.kernel$fhat
  delta.x2<-x2[2]-x2[1]
  delta.x1<-x1[2]-x1[1]
  probs<-sd.beta.kernel$fhat*delta.x2*delta.x1
  probs<-probs/sum(probs)
  
  x1.draw<-rep(NA,N)
  x2.draw<-rep(NA,N)
  # get marginal pdf of x1
  x1.marg<-rowSums(probs)
  # get CDF
  x1.cdf<-rep(NA,length(x1.marg))
  x1.cdf<-0
  for (i in 2:length(x1.marg)){
    x1.cdf[i]<-sum(x1.marg[1:i])
  }
  # make sure sums to 1
  x1.cdf<-x1.cdf/max(x1.cdf)
  
  for (n in 1:N){
    # get a draw for x1
    p<-runif(1)
    over.index<-which(x1.cdf>=p)[1]
    test.x<-x1[(over.index-1):over.index]
    test.y<-x1.cdf[(over.index-1):over.index]
    x1.draw[n]<-approx(test.y,test.x,xout=p)$y
    
    ##### Draw x2, conditional on x1
    
    # find x1's on either side of this
    x1.low<-rev(which(x1<=x1.draw[n]))[1]
    x1.high<-which(x1>=x1.draw[n])[1]
    
    
    # cycle thorugh possible x2s, interpolating density at xdraw based on x1.low and x1.high
    x2.dens.interp<-rep(NA,length(x2))
    for (i in 1:length(x2)){
      fhat.low<-fhat[x1.low,i]
      fhat.high<-fhat[x1.high,i]
      test.x<-x1[c(x1.low,x1.high)]
      x2.dens.interp[i]<-approx(test.x,c(fhat.low,fhat.high),xout=x1.draw[n])$y
    }
    # Convert x2.density to probability
    x2.prob<-x2.dens.interp*delta.x2/sum(x2.dens.interp*delta.x2)
    
    # get CDF and draw x2 from this
    
    x2.cdf<-rep(NA,length(x1.marg))
    x2.cdf[1]<-0
    for (i in 2:length(x2)){
      x2.cdf[i]<-sum(x2.prob[1:i])
    }
    x2.cdf<-x2.cdf/max(x2.cdf)
    
    p<-runif(1)
    over.index<-which(x2.cdf>=p)[1]
    test.x<-x2[(over.index-1):over.index]
    test.y<-x2.cdf[(over.index-1):over.index]
    x2.draw[n]<-approx(test.y,test.x,xout=p)$y
  }
  return(cbind(x1.draw,x2.draw))
}
# function to return cdf from an empirical distribution, returned at specific values 

getcdf<-function(x,x.ints){
  CDF<-rep(NA,length(x.ints))
  n.x<-length(x)
  CDF[1]<-0
  for (i in 2:length(x.ints)){
    CDF[i]<-length(which(x<=x.ints[i]))/n.x
  }
  return(CDF)
}
# Function to get estimated proportions above and below threshold values
est.p<-function(ts.data,lower,upper){
  # set up kernel density parameters
  # Trim years that have missing data
  no.na<-which(!apply(ts.data,1,FUN=function(x) any(is.na(x))))
  ts.data<-ts.data[no.na,]
  percentiles<-c(seq(0.0025,0.1,by=0.0025),seq(0.9,0.9975,by=0.0025))
  omega=0.5
  n.stocks=nrow(ts.data)
  n.years<-ncol(ts.data)
  kernel.output<-matrix(NA,nrow=n.years,ncol=length(percentiles))
  kernel.range<-c(0,6)
  p.below.threshold<-rep(NA,n.years)
  p.above.threshold<-rep(NA,n.years)
  for (i in 1:n.years){
    year.data<-ts.data[,i]
    # use bkde kernel density smoother, which automatically fits bandwidth based on variance in x
    #smoothed<-bkde(year.data,kernel="normal",range.x=kernel.range)
    smoothed<-density(year.data,kernel="gaussian",from=0,to=kernel.range[2])
    deltax<-smoothed$x[2]-smoothed$x[1]
    prob<-smoothed$y*deltax # turn density into probabilities
    # adjust so all density is betweeen 0 and 1
    prob<-prob/sum(prob)
    # get cumulative probabilities
    cumprob<-rep(NA,length(prob))
    for (j in 1:length(prob)){
      cumprob[j]<-sum(prob[1:j])
    }
    # get kernel output at specific percentiles
   kernel.output[i,]=approx(cumprob,smoothed$x,percentiles,yleft=0)$y
  }
  # Apply a exponential smoother (omega describes the discount rate of past estimates) for each percentile
  smoothed.kernel<-matrix(NA,nrow=n.years,ncol=length(percentiles))
  smoothed.kernel[1,]<-kernel.output[1,]
  year.rep<-seq(1,n.years,b=1) # just a vector of years
  omega.fract<-(1-omega)/(1+omega)
  for (i in 2:n.years) {
    w=omega.fract*omega^(seq(i-1,0,by=-1))
    w.prime=w/sum(w)
    w.mat<-matrix(w.prime,nrow=i,ncol=length(percentiles),byrow=FALSE)
    smoothed.kernel[i,]<-colSums(w.mat*kernel.output[1:i,])
  }
  
  threshold<-lower
  for (i in 1:n.years){
    year.estimates<-smoothed.kernel[i,]
      threshold.index<-which(year.estimates>=threshold)[1]
      p.below.threshold[i]<-0
      if(threshold.index>1&length(threshold.index)>0){
    y<-percentiles[threshold.index:(threshold.index-1)]
    x<-year.estimates[threshold.index:(threshold.index-1)]
    p.below.threshold[i]<-approx(x=x,y=y,xout=threshold)$y 
    }
  }
  # calculate # of stocks above b>upper threshold
  threshold<-upper
  for (i in 1:n.years){
    year.estimates<-smoothed.kernel[i,]
    threshold.index<-which(year.estimates>=threshold)[1]
    p.above.threshold[i]<-0
    if(!is.na(threshold.index)){
    y<-percentiles[threshold.index:(threshold.index-1)]
    x<-year.estimates[threshold.index:(threshold.index-1)]
    p.above.threshold[i]<-1-approx(x=x,y=y,xout=threshold)$y
    }
  }
  return(cbind(p.below.threshold,p.above.threshold))
}

# function to generate colored time series
color.ts<-function(N,beta){
  # this is based on Kim Cuddington's code.  You generate twice as much time series as you need, and then chop off the bottom half.  This is done to deal with the symmetric nature of the simulated data
  DIM<-N*2
  f<-c(seq(0,floor(DIM/2)),seq((ceiling(DIM/2)-1),1, by=-1)) /DIM
  phases<-runif(N,min=0,max=2*pi)
  af<-(f)^(-beta/2)
  af[af==Inf]<-0 # Replace infinity values with zero
  FFT.inv.coef<-af*(cos(phases)+1i*sin(phases))
  x.sim<-(Re(fft(FFT.inv.coef,inverse=TRUE))/N)[1:N]
  return(x.sim)
}

# Function to apply multiple segment method to calculated beta of 1/f^beta

calc.beta<-function(ts.data){
  test.lengths<-2^seq(3,6)
  n.ts<-length(ts.data)
  
  test.lengths.2.use<-which(test.lengths<n.ts)
  # make life easy and just grow the vectors
  ts.length.list<-c()
  beta.list<-c()
  
  for (i in 1:length(test.lengths.2.use)){
    short.ts.length<-test.lengths[test.lengths.2.use[i]]
    number.segments<-n.ts-short.ts.length-1
    for (j in 1:number.segments){
      tmp.ts<-ts.data[j:(j+short.ts.length-1)]
      tmp.spec<-spectrum(tmp.ts,plot=FALSE)
      freq<-tmp.spec$freq
      spec<-tmp.spec$spec
      tmp.data.frame<-data.frame(freq=log(freq),spec=log(spec))
      lm.out<-lm(spec~freq,data=tmp.data.frame)
      ts.length.list<-c(ts.length.list,short.ts.length)
      beta.list<-c(beta.list,-lm.out$coef[2])
    }
    
  }
  # fit the parameters a and n
  
  inverse.sqrt.n<-1/sqrt(ts.length.list)
  tmp.data.frame<-data.frame(inverse.n=inverse.sqrt.n,g.n=beta.list)
  lm.out<-lm(g.n~inverse.n,data=tmp.data.frame)
  a<-lm.out$coef[1]
  b<-lm.out$coef[2]
  beta<-as.numeric(a)
  CI<-as.numeric(confint(lm.out, level=0.95))
  SE<-as.numeric(coef(summary(lm.out))[1,2])
  return(c(beta,CI[1],CI[3],SE))
}

# fast functon to calculate beta, does not use the multiple segment method
calc.beta.fast<-function(ts.data){
  
  tmp.spec<-spectrum(ts.data,plot=FALSE)
  freq<-tmp.spec$freq
  spec<-tmp.spec$spec
  tmp.data.frame<-data.frame(freq=log(freq),spec=log(spec))
  lm.out<-lm(spec~freq,data=tmp.data.frame)
  beta<--lm.out$coef[2]
  return(beta)
}

getcdf<-function(x,x.ints){
  CDF<-rep(NA,length(x.ints))
  n.x<-length(x)
  CDF[1]<-0
  for (i in 2:length(x.ints)){
    CDF[i]<-length(which(x<=x.ints[i]))/n.x
  }
  return(CDF)
}

min.and.duration<-function(x,y){
  # function to return the minimum of time series x, and the number of years until time series increases by y units
  minB<-min(x)
  year.min.B<-which(x==minB)
  
  recovery.time<-NA
  # get years to recovery if minB is less than y
  if (minB<y){
    # years to recovery
    short.x<-x[-(1:year.min.B)]
    recovery.index<-which(short.x>=y)
    if (length(recovery.index>0)){
      recovery.time<-recovery.index[1]
    }
  }
  return(c(minB,year.min.B,recovery.time))
}

addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}
################ END of FUNCTION ####################


######## Load Data and Get Properties of Real Time series


# initialize data
# Set working directory to the path that contains csv file

setwd(workdir)
datafilename<-"allforagedata.csv"
data<-read.csv(file=datafilename,header=TRUE)

# Initialize Data 
stockid<-unique(data$ASSESSID)
n.stocks<-length(stockid)
min.year<-min(data$Year)
max.year<-max(data$Year)
n.years<-max.year-min.year+1
year.list<-seq(min.year,max.year,by=1)
outputB<-matrix(NA,nrow=n.years,ncol=n.stocks)
outputF<-matrix(NA,nrow=n.years,ncol=n.stocks)
outputr<-matrix(NA,nrow=n.years,ncol=n.stocks)
# go through each stock, get TB/mean(TB) or SSB/mean(SSB) by year and place in matrix
biomass<-tapply(data$TB,data$ASSESSID,mean,na.rm=TRUE)

min.biomass<-rep(NA,n.stocks)
max.biomass<-rep(NA,n.stocks)
first.year.listB<-rep(NA,n.stocks)
last.year.listB<-rep(NA,n.stocks)

# loop through stocks, gather biomass, catch data and exploitation rate
for (i in 1:n.stocks){
  stockindex<-which(data$ASSESSID==stockid[i])
  TB.tmp<-data$TB[stockindex]
  Year.tmp<-data$Year[stockindex]
  n.TB<-length(which(!is.na(TB.tmp)))
  
  #first.year.TB<-Year.tmp[which(!is.na(TB.tmp))[1]]
  last.year<-max(Year.tmp)
  # last.year.listB[i]<-last.year
  first.year<-min(Year.tmp)
  ly.index<-which(year.list==last.year)
  fy.index<-which(year.list==first.year)
  if(n.TB>0){
    B.list<-TB.tmp/mean(TB.tmp,na.rm=TRUE)
    #  first.year.listB[i]<-first.year.TB
    outputB[fy.index:ly.index,i]<-B.list
  }
  ER.tmp<-data$ER[stockindex]
  n.ER<-length(which(!is.na(ER.tmp)))
  #first.year.ER<-Year.tmp[which(!is.na(ER.tmp))[1]]
  last.year<-max(Year.tmp)
  #last.year.listER[i]<-last.year
  first.year<-min(Year.tmp)
  ly.index<-which(year.list==last.year)
  fy.index<-which(year.list==first.year)
  if(n.ER>0){
    #  first.year.listER[i]<-first.year.ER
    outputF[fy.index:ly.index,i]<-ER.tmp
  }
} 
colnames(outputB)<-stockid
rownames(outputB)<-seq(min.year,max.year)



colnames(outputF)<-stockid
rownames(outputF)<-seq(min.year,max.year)

# get timeseries of population growth rate
calc.r<-function(B,F) log(B[-1]/B[-length(B)]+F[-length(F)])
outputr<-matrix(NA,nrow=nrow(outputB)-1,ncol=ncol(outputB))
for (i in 1:ncol(outputB)){
  outputr[,i]<-calc.r(outputB[,i],outputF[,i])
}

yearcount<-function(x) length(which(!is.na(x)))
n.years<-apply(FUN=yearcount,X=outputB,2)
long.index<-which(n.years>=25)
stock.count<-apply(FUN=yearcount,X=outputB[,long.index],1)


# only save output for populations that have 25 years or more of population biomass
outputB<-outputB[,long.index]
outputF<-outputF[,long.index]
outputr<-outputr[,long.index]

B.data.mat<-matrix(NA,nrow=ncol(outputB),ncol=11)
rownames(B.data.mat)<-colnames(outputB)
colnames(B.data.mat)<-c("TS Length","SD.R","Beta.r","lowerCI","upperCI","beta SE","Beta.b","lowerCI","upperCI","beta SE","SD.B")
# outputB saved from B and F summary, includes all data

outputB.all<-outputB

###### Load data, extract total biomass ONLY for calculation of lambdas


# get total number of years of all biomass ts (total biomass and SSB)
nyear.calc<-function(x) length(which(!is.na(x)))
nyears.list<-apply(outputB.all,2,nyear.calc)
# minimum number of years for calculating lambda properties
min.years<-33

# Calculate properties of time series for which we have total biomass and exploitation rate
for ( i in 1:ncol(outputB)){
  Bindex<-intersect(which(!is.na(outputB[,i])),which(!is.na(outputF[,i])))
  ts.data<-outputB[Bindex,i]
  ts.f.data<-outputF[Bindex,i]
  B.data.mat[i,1]<-length(which(!is.na(outputB[,i])))
  if (length(ts.data)>=min.years){
    # Estimate lambda and r
    lambda<- ts.data[-1]/ts.data[-length(ts.data)]+ts.f.data[-length(ts.f.data)]
    r<-log(lambda)
    # get spectral propoerties
    beta.r<-calc.beta(r)
    beta.ts<-calc.beta(ts.data)
    # save output
    B.data.mat[i,2]<-sd((r))
    B.data.mat[i,3:6]<-beta.r
    B.data.mat[i,7:10]<-beta.ts
    B.data.mat[i,11]<-sd(ts.data)
  }
}

#  Get kernel density smoothers of damn near everything
estimated.beta<-B.data.mat[which(!is.na(B.data.mat[,3])),3]
# set negative values to 0
estimated.beta<-replace(estimated.beta,which(estimated.beta<0),0)
estimated.beta.weights<-1/B.data.mat[which(!is.na(B.data.mat[,3])),6]
estimated.beta.weights<-replace(estimated.beta.weights,which(estimated.beta==0),4)
beta.kernel<-density(estimated.beta,kernel="gaussian",weights=estimated.beta.weights/sum(estimated.beta.weights))

estimated.sd.r<-((B.data.mat[which(!is.na(B.data.mat[,2])),2]))
sd.r.kernel<-bkde(estimated.sd.r,kernel="normal",range.x=c(0,1))

estimated.sd.B<-((B.data.mat[which(!is.na(B.data.mat[,11])),11]))
sd.B.kernel<-bkde(estimated.sd.B,kernel="normal",range.x=c(0,2))

estimated.beta.B<-B.data.mat[which(!is.na(B.data.mat[,7])),7]
estimated.beta.B.weights<-1/B.data.mat[which(!is.na(B.data.mat[,3])),10]
beta.B.kernel<-density(estimated.beta.B,kernel="gaussian",weights=estimated.beta.B.weights/sum(estimated.beta.B.weights))

# do a two dimensional KDE for beta and sd of r
d<-2 # the bandwidth dimension
sd.bw<-sqrt(var(estimated.sd.r))*(4/((d+2)*length(estimated.sd.r)))^(1/(d+4))
beta.bw<-sqrt(var(estimated.beta))*(4/((d+2)*length(estimated.beta)))^(1/(d+4))

sd.beta.kernel<-bkde2D(cbind(estimated.beta,estimated.sd.r),bandwidth=c(beta.bw,sd.bw),range.x=list(c(-1.0,3.5),c(0,.8)),truncate=FALSE,gridsize=c(200,200))

######## End of Loading Data  ###########

#########Begin Randomization Test ##############
# setup information
ts.lengths<-as.numeric(nyears.list[which(nyears.list>25)])
n.stocks<-55 # number of stocks with 25 or more years of data
n.sims<-1000

# tolerance on simulated B time series, in terms of beta and sd
beta.B.tol<-c(0.5,4.0)
sd.B.tol<-c(0.3,1.25)
# setup vectors for storing output
minblist<-seq(0,1,by=0.025)
maxblist<-seq(1,5,by=.2)
min.output<-matrix(NA,nrow=n.sims,ncol=length(minblist))
max.output<-matrix(NA,nrow=n.sims,ncol=length(maxblist))
beta.B.bin<-seq(-0.1,4,by=0.05)
sd.B.bin<-seq(0,1.5,by=0.025)
cond.minb.output<-matrix(NA,nrow=n.sims,ncol=1)
beta.output<-matrix(NA,nrow=n.sims,ncol=length(beta.B.bin))
sd.output<-matrix(NA,nrow=n.sims,ncol=length(sd.B.bin))
p.collapse<-matrix(NA,nrow=n.sims,ncol=2)
p.uber.collapse<-matrix(NA,nrow=n.sims,ncol=2)
duration.output<-rep(NA,n.sims)

# Begin simulations
for (sim in 1:n.sims){
  # create a matrix to hold n.stocks time series and their properites
  ts.properties<-matrix(NA,nrow=n.stocks,ncol=3)
  beta.B.sim<-matrix(NA,nrow=n.stocks,ncol=1)
  sd.B.sim<--matrix(NA,nrow=n.stocks,ncol=1)
  # randomly sample ts.lengths from observations, with replacement
  ts.sim.lengths<-sample(ts.lengths,n.stocks,replace=TRUE)
  min.b<-matrix(NA,nrow=n.stocks,ncol=1)
  ts.sims<-matrix(NA,nrow=27,ncol=n.stocks)
  for (i in 1:n.stocks){
    # Get t.series with specified beta within tolerance range (0.5, 3.0)
    flag=0
    # generate a large time series to get scale it down so that it has mean 1
    t.mean.calc<-ts.sim.lengths[i]
    while (flag==0){
      # get random draw from 2D kernel density of beta and sd
      beta.sd<-rkernel2D(1,sd.beta.kernel)
      beta<-max(0,beta.sd[1])
      target.sd<-beta.sd[2]
      # Simulate time series and get it's beta
      r.sim<-(color.ts(N=150,beta=beta))
      test.beta.r<-calc.beta.fast(r.sim)[1]
      # only proceed if test.beta.r is "good", within 0.25 of target
      beta.diff<-abs(test.beta.r-beta)
      if (beta.diff<=0.25) {
        # adjust variance of r.sim to achieve target SD
        r.sim.sd<-sd(r.sim)
        r.sim.adj<-r.sim*(target.sd/r.sim.sd)
        # correct to mean =0
        r.sim.adj<-r.sim.adj-rep(mean(r.sim.adj),150)
        # turn this into a bunch of abundances, with initial population size=1
        ts.sim<-rep(NA,150)
        ts.sim[1]<-1
        for (t in 2:150){
          ts.sim[t]<-exp(r.sim.adj[t-1])*ts.sim[t-1]
    
        }
        # get a sample measured over last t.mean.calc years
        ts.sample<-ts.sim[(150-t.mean.calc+1):150]
        # get a standarized biomass
        adj.sample<-ts.sample/mean(ts.sample)
        test.beta.B<-calc.beta.fast(ts.sim)
        test.sd.B<-sd(ts.sim)/mean(ts.sim)
        # Determine whether generated time series matches biomass time series properites
        beta.flag=(test.beta.B<beta.B.tol[2]&test.beta.B>beta.B.tol[1])
        sd.flag<-(test.sd.B<sd.B.tol[2]&test.sd.B>sd.B.tol[1])
        # if conditions are met, save output in ts.output and move onto next stock
        if (beta.flag&sd.flag){
          flag=1
          ts.properties[i,1]<-min(adj.sample)
          ts.properties[i,2]<-max(adj.sample)
          if (min(adj.sample)<0.25) ts.properties[i,3]<-min.and.duration(adj.sample,1)[3]
          beta.B.sim[i]<-calc.beta.fast(adj.sample)
          sd.B.sim[i]<-sd(adj.sample)
          # save the output into matrix called ts.sims
          # different rules depending on the length of time series
          if (t.mean.calc>=27) {
          ts.sims[,i]<-adj.sample[(t.mean.calc-27+1):t.mean.calc] # take the last 32 years
          } else {
            ts.sims[(27-t.mean.calc+1):27,i]<-adj.sample
          }
        }
      }
    }
  }
  # get CDF of minB
  cdf<-getcdf(ts.properties[,1],minblist)
  min.output[sim,]<-cdf
  cdf<-getcdf(ts.properties[,2],maxblist)
  max.output[sim,]<-cdf
  beta.cdf<-getcdf(beta.B.sim,beta.B.bin)
  sd.cdf<-getcdf(sd.B.sim,sd.B.bin)
  beta.output[sim,]<-beta.cdf
  sd.output[sim,]<-sd.cdf
  duration.output[sim]<-mean(ts.properties[,3],na.rm=T)
  # mean minimum biomass for those stocks that were below 0.25
  cond.minb.output[sim,]<-mean(ts.properties[which(ts.properties<=0.2)])
  # figure out proportion of stocks collapsed by year.  Use only 52 stocks because there are 52 stocks in 2006, the year with the most observed tocks
  n.stocks.year<-apply(ts.sims[,1:52],1,FUN=function(x) length(which(!is.na(x))))
  n.stocks.collapse<-apply(ts.sims[,1:52],1,FUN=function(x) length(which(x<0.25)))
  collapse.bonanza<-cbind(n.stocks.collapse/n.stocks.year,n.stocks.collapse/n.stocks.year)
  p.collapse[sim,1:2]<-c(min(collapse.bonanza[,1]),max(collapse.bonanza[,1]))
  uber.collapse<-apply(ts.sims[,1:52],1,FUN=function(x) length(which(x<0.2)))/n.stocks.year
  p.uber.collapse[sim,1:2]<-c(min(uber.collapse[1]),max(uber.collapse[1]))
}

######End Randomization Algorithm


save(file="RandomizationTestlambda.Rdata",min.output,max.output,beta.output,sd.output,cond.minb.output,duration.output,p.collapse,p.uber.collapse)

#######################End of the randomization test

### Clear memory of all variables except working directory
rm(list=ls()[-which(ls()=="workdir")])

########## load new Functions####################
getcdf<-function(x,x.ints){
  CDF<-rep(NA,length(x.ints))
  n.x<-length(x)
  CDF[1]<-0
  for (i in 2:length(x.ints)){
    CDF[i]<-length(which(x<=x.ints[i]))/n.x
  }
  return(CDF)
}

min.and.duration<-function(x,y){
  # function to return the minimum of time series x, and the number of years until time series increases by y units
  minB<-min(x)
  year.min.B<-which(x==minB)
  
  recovery.time<-NA
  # get years to recovery if minB is less than y
  if (minB<y){
    # years to recovery
    short.x<-x[-(1:year.min.B)]
    recovery.index<-which(short.x>=y)
    if (length(recovery.index>0)){
      recovery.time<-recovery.index[1]
    }
  }
  return(c(minB,year.min.B,recovery.time))
}

addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}
###########End of functions #######
###Reload data
load("B and F summary.Rdata") # This data file is created by running the script "summarize B and F data.r"
load("RandomizationTestlambda.Rdata")

##Compare randomization to data

# Extract properties of Actual biomass timeseries
minB.output<-matrix(NA,nrow=ncol(outputB),ncol=4)
yearlist<-as.numeric(rownames(outputB))
colnames(minB.output)<-c("minB","year of minB","year collapse","duration")
for ( i in 1:ncol(outputB)){
  Bindex<-which(!is.na(outputB[,i]))
  Bindex<-intersect(which(!is.na(outputB[,i])),which(!is.na(outputF[,i])))
  ts.data<-outputB[Bindex,i]
  ts.f.data<-outputF[Bindex,i]
  #B.data.mat[i,1]<-length(ts.data)
  if (length(ts.data)>=25){
    minB.output[i,1]<-min(ts.data)
    minB.output[i,2]<-max(ts.data)
    minB.output[i,3]<-yearlist[which(outputB[,i]==min(ts.data))]
    minB.output[i,4]<-min.and.duration(ts.data,1)[3]
  }
}


# setup plot

par(omi=c(1,1,1,1),las=1,cex.axis=1.5,cex.lab=1.5)
layout(matrix(c(1,2,2),nrow=3,ncol=1))

# get histogram of minimum biomass
minblist<-seq(0,1,by=0.1)
minb.hist<-hist(minB.output[,1],breaks=minblist,right=FALSE,plot=FALSE)
top<-15
par(mai=c(0.2,.75,.75,0.30))
par(xpd=NA)
barplot(minb.hist$counts,axes=FALSE, xaxs="i", xlim=c(0,10.5),ylim=c(0, 10), space=0,col="#74ADD1",)
axis(2,at=c(0,5,10),labels=TRUE)
par(las=0)
mtext(side=2,text="Count",line=3)
par(las=1)

# get observed cumulative distribution of minimum biomass
obs.minblist<-seq(0,1,by=0.01)
min.B.index<-which(!is.na(minB.output[,1]))
obs.minb.cdf<-getcdf(minB.output[min.B.index,1],obs.minblist)
par(mai=c(0.75,0.75,.1,0.46))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),lwd=2,xlab="Minimum Stock Biomass",ylab="Cumulative Probability",main="",yaxs="i",xaxs="i")
# use the ablines to line up histogram with bottom graphic
# get upper /lower percentiles of cumulative distrubitions
lowFUN<-function(x) quantile(x,0.025)
upFUN<-function(x) quantile(x,0.975)
medFUN<-function(x) quantile(x,0.5)
minblist<-seq(0,1,by=0.025)
lower.bmin<-apply(min.output,2,lowFUN)
upper.bmin<-apply(min.output,2,upFUN)

# Generate plots

# Make polygon on 95% Ci
poly.col<-addTrans("#FDAE61",125)
x.poly<-c(minblist,rev(minblist))
y.poly<-c(lower.bmin,rev(upper.bmin))
polygon(x.poly,y.poly,col=poly.col,border=poly.col,lwd=.5)
lines(obs.minblist,obs.minb.cdf,type="l",col="#74ADD1",lwd=2)



# what is mean lower biomass of stocks that dropped below 0.25
mean.min.B<-mean(minB.output[which(minB.output<=0.2),1])
percent.lower<-length(which(cond.minb.output<=mean.min.B))/(length(cond.minb.output))

# fraction of stocks collapsed in 2006 (n=52)

upper.pcollapse<-upFUN(p.collapse)
# What is duration of collapse
collapse.index<-which(minB.output[,1]<=0.25)
mean(minB.output[collapse.index,4])
