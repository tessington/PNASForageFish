# Code to summarize and collate biomass and fishing mortality rate data that are contained in csv file "allforagedata.csv"

wordir<- # set path to directory containing allforagedata.csv
  
# Summarize biomass data (as B:mean B for all stocks in the RAM legacy database, save in matrix called outputB
setwd(workdir)

# load data
#workdir<-"/Users/essington/Desktop/Desktop/Rcode/ForageFish/Forage Status Paper"
datafilename<-"allforagedata.csv"

setwd(workdir)
data<-read.csv(file=datafilename,header=TRUE)

stockid<-unique(data$ASSESSID)
n.stocks<-length(stockid)
min.year<-min(data$Year)
max.year<-max(data$Year)
n.years<-max.year-min.year+1
year.list<-seq(min.year,max.year,by=1)
outputB<-matrix(NA,nrow=n.years,ncol=n.stocks)
# go through each stock, get TB/mean(TB) or SSB/mean(SSB) by year and place in matrix
biomass<-tapply(data$TB,data$ASSESSID,mean,na.rm=TRUE)
SSBbiomass<-tapply(data$SSB,data$ASSESSID,mean,na.rm=TRUE)
min.biomass<-rep(NA,n.stocks)
max.biomass<-rep(NA,n.stocks)
first.year.listB<-rep(NA,n.stocks)
last.year.listB<-rep(NA,n.stocks)
outputP<-matrix(NA,nrow=n.years,ncol=n.stocks)

for (i in 1:n.stocks){
  stockindex<-which(data$ASSESSID==stockid[i])
  TB.tmp<-data$TB[stockindex]
  SSB.tmp<-data$SSB[stockindex]
  Year.tmp<-data$Year[stockindex]
  TC.tmp<-data$TC[stockindex]
  TL.tmp<-data$TC[stockindex]
  TC<-TL.tmp
  if(any(!is.na(TC.tmp))) {
    TC<-TC.tmp
  }
  # figure out whether to use TB.tmp or SSB.smp 
  n.TB<-length(which(!is.na(TB.tmp)))
  n.SSB<-length(which(!is.na(SSB.tmp)))
  first.year.TB<-Year.tmp[which(!is.na(TB.tmp))[1]]
  first.year.SSB<-Year.tmp[which(!is.na(SSB.tmp))[1]]
  last.year<-max(Year.tmp)
  last.year.listB[i]<-last.year
  first.year<-min(Year.tmp)
  ly.index<-which(year.list==last.year)
  fy.index<-which(year.list==first.year)
  if (n.TB==0) {
    B.list<-SSB.tmp/mean(SSB.tmp,na.rm=TRUE)
    first.year.listB[i]<-first.year.SSB
  } else {
    B.list<-TB.tmp/mean(TB.tmp,na.rm=TRUE)
    first.year.listB[i]<-first.year.TB
  }
  outputB[fy.index:ly.index,i]<-B.list
  min.biomass[i]<-min(B.list.na.rm=TRUE)
  max.biomass[i]<-max(B.list,na.rm=TRUE)
} 
colnames(outputB)<-stockid
rownames(outputB)<-seq(min.year,max.year)
## Repeat for Exploitation Rate
outputF<-matrix(NA,nrow=n.years,ncol=n.stocks)
first.year.listF<-rep(NA,n.stocks)
last.year.listF<-rep(NA,n.stocks)
# turn F into an approximate Exploitation Rate 1-exp(-F)
data$F<-rep(1,nrow(data))-exp(-data$F)
for (i in 1:n.stocks){
  stockindex<-which(data$ASSESSID==stockid[i])
  F.tmp<-data$F[stockindex]
  ER.tmp<-data$ER[stockindex]
  
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

save(file="B and F summary.Rdata",outputB,outputF,outputr,first.year.listB,last.year.listB,first.year.listF,last.year.listF,min.biomass,max.biomass)


#start.year<-1975
#end.year<-2006



#n.stocks.continuous<-function(x,start.year,end.year) {
##  start.year.index<-which(rownames(x)==start.year)
#  end.year.index<-which(rownames(x)==end.year)
  
#  continuous.list<-c()
#  for (i in 1:ncol(x)){
#    if (!any(is.na(x[start.year.index:end.year.index,i]))) continuous.list<-c(continuous.list,i)
#  }
#  return(list(n.stocks=length(continuous.list),stock.list=continuous.list))
#}


##start.year.list<-seq(1970,1985)
#end.year.list<-seq(2000,2010)
#n.stocks.matrix<-matrix(NA,nrow=length(start.year.list),ncol=length(end.year.list))

#for (i in 1:length(start.year.list)){
#  for (j in 1:length(end.year.list)){
#    n.stocks.output<-n.stocks.continuous(outputB,start.year=start.year.list[i],end.year=end.year.list[j])
#    n.stocks.matrix[i,j]<-n.stocks.output$n.stocks
#  }
#}

#rownames(n.stocks.matrix)<-start.year.list
#colnames(n.stocks.matrix)<-end.year.list


#### generate a table that lists first and last year of biomass / F data and number of years
#return.first.last<-function(x){
#  yearlist<-1907:2013
# c(yearlist[min(which(!is.na(x)))],yearlist[max(which(!is.na(x)))])
#}
#year.index<-apply(FUN=return.first.last,MARGIN=2,X=outputB)
