# code to simulate consequences of hypothetical harvest control rule, assuming that production is independent of biomass
# This is the version used in the paper

wordir<- # set path to directory containing allforagedata.csv

setwd(workdir)
datafilename<-"allforagedata.csv"

data<-read.csv(file=datafilename,header=TRUE)
###### Load data, extract total biomass 
stockid<-unique(data$ASSESSID)
n.stocks<-length(stockid)
min.year<-min(data$Year)
max.year<-max(data$Year)
n.years<-max.year-min.year+1
year.list<-seq(min.year,max.year,by=1)
outputB<-matrix(NA,nrow=n.years,ncol=n.stocks)
# go through each stock, get TB/mean(TB) or SSB/mean(SSB) by year and place in matrix
biomass<-tapply(data$TB,data$ASSESSID,mean,na.rm=TRUE)
min.biomass<-rep(NA,n.stocks)
max.biomass<-rep(NA,n.stocks)
first.year.listB<-rep(NA,n.stocks)
last.year.listB<-rep(NA,n.stocks)


for (i in 1:n.stocks){
  stockindex<-which(data$ASSESSID==stockid[i])
  TB.tmp<-data$TB[stockindex]
  Year.tmp<-data$Year[stockindex]
  # figure out whether to use TB.tmp or SSB.smp (whichever is longer)
  n.TB<-length(which(!is.na(TB.tmp)))
  
  first.year.TB<-Year.tmp[which(!is.na(TB.tmp))[1]]
  last.year<-max(Year.tmp)
  last.year.listB[i]<-last.year
  first.year<-min(Year.tmp)
  ly.index<-which(year.list==last.year)
  fy.index<-which(year.list==first.year)
  if(n.TB>25){
    B.list<-TB.tmp/mean(TB.tmp,na.rm=TRUE)
    first.year.listB[i]<-first.year.TB
    outputB[fy.index:ly.index,i]<-B.list
    min.biomass[i]<-min(B.list,na.rm=TRUE)
    max.biomass[i]<-max(B.list,na.rm=TRUE)
  }
  
} 
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

# which stocks had a minimum biomass below threshold
bcollapse<-0.25 #Set collapse biomass threshold
blim=0.5 # Set lower biomass limit where fishing equals zero


# Initialize vectors to store model output
collapse.index<-which(min.biomass<=blim) #Which stocks are vulnerable to Collapse by having a lower Biomass less than .5
lost.yield<-rep(NA,length(collapse.index))
new.min.B<-rep(NA,length(collapse.index))
old.min.B<-rep(NA,length(collapse.index))
lost.yield.overall<-rep(NA,length(collapse.index))
for (i in 1:length(collapse.index)){
  stock.index<-collapse.index[i]
  B.list<-outputB[,stock.index]
  F.list<-outputF[,stock.index]
  C.list<-B.list*F.list
  SP<-B.list[-1]-B.list[-length(B.list)]+C.list[-length(B.list)]
  # find first year of Where biomass Is below The biomass limit
  year.index<-which(B.list<=blim)[1]
  total.catch<-sum(C.list[(year.index):(year.index+9)]) #
  # total catch minus collapse Years
  total.catch.other.years<-sum(C.list[-((year.index):(year.index+9))],na.rm=T)
  new.F.list<-F.list[(year.index):(year.index+9)]
  new.B.list<-rep(NA,10)
  new.B.list[1]<-B.list[(year.index)]
  new.F.list[1]<-0
  new.SP.list<-SP[(year.index):(year.index+9)]
  new.C.list<-C.list[(year.index):(year.index+9)]
  # proceed if there is three years post and pre crash
  if(!any(is.na(new.SP.list[1:9]))){
  for (j in 2:10){
    new.B.list[j]<-new.B.list[j-1]+new.SP.list[j-1]-new.B.list[j-1]*new.F.list[j-1]
    if (new.B.list[j]<=blim) new.F.list[j]=0
  }
  
  total.catch.newF<-sum(new.B.list*new.F.list)
  overall.total.catch.newF<-total.catch.newF+total.catch.other.years
  overall.total.catch<-total.catch+total.catch.other.years
  lost.yield[i]<-(total.catch-total.catch.newF)/total.catch
  lost.yield.overall[i]<-(overall.total.catch-overall.total.catch.newF)/overall.total.catch
  new.min.B[i]<-min(new.B.list)
  old.min.B[i]<-min(B.list[(year.index):(year.index+10)])
  }
}



# Overall lost yield
mean(lost.yield.overall,na.rm=T)

new.min.mean<-mean(new.min.B,na.rm=T)
old.min.mean<-mean(old.min.B,na.rm=T)
mean((new.min.B-old.min.B)/old.min.B,na.rm=T)
mean((new.min.B)/old.min.B,na.rm=T)

# How many Stocks Have a minimum biomass below collapsed threshold
sim.collapse<-length(which(round(new.min.B,2)<bcollapse))
actual.collapse<-length(which(round(old.min.B,2)<=bcollapse))
