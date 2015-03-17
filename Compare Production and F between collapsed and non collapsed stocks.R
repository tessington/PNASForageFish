### This is code Compares the surplus production and exploitation rate between collapsed and non-collapsed stocks
# delete all objects from memory


datafilename<-"allforagedata.csv"
workdir<-# set workdir to folder where allforagedata.csv is stored
setwd(workdir) 
data<-read.csv(file=datafilename,header=TRUE)

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


# stocks that have at least 25 years of both B and F data
count.years.stocks<-function(x) length(which(!is.na(x)))

n.years.B<-apply(outputB,2,count.years.stocks)
n.years.F<-apply(outputF,2,count.years.stocks)


B.index<-which(n.years.B>=25)
F.index<-which(n.years.F>=25)
stock.index<-intersect(B.index,F.index)
years<-rownames(outputB)

# round B data down to nearest 0.01:
outputB<-floor(100*outputB)/100

# get minB for each stock
# find stocks that collapsed, and get collapse year
collapse.year.fun<-function(X) which(X<=.25)[1]

collapse.year<-apply(outputB[,stock.index],2,FUN=collapse.year.fun)
minB.stock<-apply(outputB[,stock.index],2,min,na.rm=T)
minB.years<-apply(outputB[,stock.index],MARGIN=2,FUN=which.min)

mean.SP<-rep(NA,length(minB.years))
mean.F<-rep(NA,length(minB.years))

for (i in 1:length(mean.SP)){
  if (!is.na(collapse.year[i])){
    mean.SP[i]<-mean(outputSP[c(collapse.year[i]-2,collapse.year[i]-1),stock.index[i]])
    mean.F[i]<-mean(outputF[c(collapse.year[i]-2,collapse.year[i]-1),stock.index[i]])
  } else {
    
  mean.SP[i]<-mean(outputSP[c(minB.years[i]-2,minB.years[i]-1),stock.index[i]])
  mean.F[i]<-mean(outputF[c(minB.years[i]-2,minB.years[i]-1),stock.index[i]])
}
}

# get collapsed stocks
collapse.threshold<-0.25
collapse.index<-which(minB.stock<=collapse.threshold)
# This analysis  only those stocks that had enough data to allow plotting in Figure 2 (basically had longer pre- and post-collapse data)
collapse.index<-c(1,  2,  3, 10, 11, 12, 14, 15, 17, 18, 22, 24, 25, 39, 41, 42, 44)

# get non collapsed stocks
no.collapse.threshold<-0.3
no.collapse.index<-which(minB.stock>=no.collapse.threshold)
t.test(mean.SP[collapse.index],mean.SP[no.collapse.index],alterantive="two.sided",var.equal=F)
t.test(mean.F[collapse.index],mean.F[no.collapse.index],alterantive="two.sided",var.equal=F)
