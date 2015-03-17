# A series of  analyses on stock collapses  to ask whether the rate of collapses are different by decade or by region


workdir<- # Set working directory to the path that contains csv file

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
# go through each stock, get TB/mean(TB) or SSB/mean(SSB) by year and place in matrix
biomass<-tapply(data$TB,data$ASSESSID,mean,na.rm=TRUE)
SSBbiomass<-tapply(data$SSB,data$ASSESSID,mean,na.rm=TRUE)
min.biomass<-rep(NA,n.stocks)
max.biomass<-rep(NA,n.stocks)
first.year.listB<-rep(NA,n.stocks)
last.year.listB<-rep(NA,n.stocks)

for (i in 1:n.stocks){
  stockindex<-which(data$ASSESSID==stockid[i])
  TB.tmp<-data$TB[stockindex]
  SSB.tmp<-data$SSB[stockindex]
  Year.tmp<-data$Year[stockindex]
  # get total catch (TC) or total landings (TL) and save.  prefentially use TC if available
  TC.tmp<-data$TC[stockindex]
  TL.tmp<-data$TC[stockindex]
  TC<-TL.tmp
  if(any(!is.na(TC.tmp))) {
    TC<-TC.tmp
  }
  # figure out whether to use total biomass or spawning stock biomass (preferntially use total biomass
  
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

yearcount<-function(x) length(which(!is.na(x)))
n.years<-apply(FUN=yearcount,X=outputB,2)
long.index<-which(n.years>=25)
stock.count<-apply(FUN=yearcount,X=outputB[,long.index],1)


# only save output for populations that have 25 years or more of population biomass
outputB<-outputB[,long.index]
outputF<-outputF[,long.index]



############## Conduct analysis of stock collapses #############################
n.stockdata<-function(x) length(which(!is.na(x)))
n.stock.collapse<-function(x) length(x<0.25)
first.year<-function(x) as.numeric(names(x))[which(!is.na(x))[1]]
min.year<-function(x) which(x==min(x,na.rm=T))
min.year.collapse<-function(x) which(x<=.25)[1] # change the value here to get different collapse thresholds (e.g. 0.15, 0.25)
min.years<-apply(outputB,2,n.stockdata)
stock.index<-which(min.years>=25)
outputB<-outputB[,stock.index]
minB<-apply(outputB,2,min,na.rm=T)

n.stocks.collapse.25<-length(which(minB<0.25))
n.stocks.collapse.15<-length(which(minB<0.15))

collapse.thresold=0.25
n.stocks.year<-apply(outputB,1,n.stockdata)
n.stocks.col<-apply(outputB,1,n.stock.collapse,collapse.threshold=collapse.thresold)
years<-as.numeric(names(n.stocks.year))


# how many stocks exhibited lowest biomass in 2000's, 1990's, 1980's, etc.


year.collapse<-as.numeric(rownames(outputB))[apply(outputB,2,min.year.collapse)]


all.years<-as.numeric(rownames(outputB))
# create matrix called data.outputB that has 1 if there is biomass data that year, 0 otherwise
data.outputB<-replace(outputB,!is.na(outputB),1)
data.outputB<-replace(data.outputB,is.na(outputB),0)
index.1950s<-intersect(which(all.years>=1950),which(all.years<1960))
index.1960s<-intersect(which(all.years>=1960),which(all.years<1970))
index.1970s<-intersect(which(all.years>=1970),which(all.years<1980))
index.1980s<-intersect(which(all.years>=1980),which(all.years<1990))
index.1990s<-intersect(which(all.years>=1990),which(all.years<2000))
index.2000s<-which(all.years>=2000)

n.stocks.1950s<-length(which(colSums(data.outputB[index.1950s,])>=5))
n.stocks.1960s<-length(which(colSums(data.outputB[index.1960s,])>=5))
n.stocks.1970s<-length(which(colSums(data.outputB[index.1970s,])>=5))
n.stocks.1980s<-length(which(colSums(data.outputB[index.1980s,])>=5))
n.stocks.1990s<-length(which(colSums(data.outputB[index.1990s,])>=5))
n.stocks.2000s<-length(which(colSums(data.outputB[index.2000s,])>=5))

n.collapse.1950s<-length(intersect(which(year.collapse>=1950),which(year.collapse<1960)))
n.collapse.1960s<-length(intersect(which(year.collapse>=1960),which(year.collapse<1970)))
n.collapse.1970s<-length(intersect(which(year.collapse>=1970),which(year.collapse<1980)))
n.collapse.1980s<-length(intersect(which(year.collapse>=1980),which(year.collapse<1990)))
n.collapse.1990s<-length(intersect(which(year.collapse>=1990),which(year.collapse<2000)))
n.collapse.2000s<-length(intersect(which(year.collapse>=2000),which(year.collapse<2013)))


decade.list<-as.factor(c(1950,1960,1970,1980,1990,2000))

n.stocks.list<-c(n.stocks.1950s,n.stocks.1960s,n.stocks.1970s,n.stocks.1980s,n.stocks.1990s,n.stocks.2000s)
n.collapse.list<-c(n.collapse.1950s,n.collapse.1960s,n.collapse.1970s,n.collapse.1980s,n.collapse.1990s,n.collapse.2000s)

# Run GLM on whether # of collapses : total number of stocks varies by decade
collapse.glm.decade<-glm(cbind(n.collapse.list,n.stocks.list)~decade.list,family=binomial)
anova(collapse.glm.decade,test="Chisq")
summary(collapse.glm.decade)

# Finally, get the number of stocks in collapse for each year
n.collapse.by.year<-apply(FUN=n.stock.collapse,X=outputB,MAR=1)
n.stocks.by.year<-apply(FUN=n.stockdata,X=outputB,MAR=1)
cbind(n.collapse.by.year,n.stocks.by.year)

########################################################################
# Test for differences in collapse frequency by region
# for each stock, get the region.
stock.names<-names(stock.index)
extract.region<-function(x,Reg) Reg[which(x==Reg[,1])[1],2]
region.id<-sapply(stock.names,extract.region,data[,c(1,14)])
region.id<-as.character(region.id)

unique.regions<-unique(region.id)
# figure out how many stocks are in each region
stocks.region<-function(X,Y) length(which(Y==X))
n.stocks.region<-sapply(unique.regions,FUN=stocks.region,Y=region.id)

# how many of these collapsed
collapse.threshold<-0.25
collapse.stocks<-region.id[which(minB<=collapse.threshold)]
n.collapse.region<-sapply(unique.regions,FUN=stocks.region,Y=collapse.stocks)
collapse.region.decade<-glm(cbind(n.collapse.region,n.stocks.region)~unique.regions,family=binomial)
anova(collapse.region.decade,test="Chisq")
summary(collapse.region.decade)


# Repeat with HIGH collapse threshold
collapse.threshold<-0.5
collapse.stocks<-region.id[which(minB<=collapse.threshold)]
n.collapse.region<-sapply(unique.regions,FUN=stocks.region,Y=collapse.stocks)
collapse.region.decade<-glm(cbind(n.collapse.region,n.stocks.region)~unique.regions,family=binomial)
anova(collapse.region.decade,test="Chisq")
summary(collapse.region.decade)

# Repead with LOW collapse threshold
collapse.threshold<-0.15
collapse.stocks<-region.id[which(minB<=collapse.threshold)]
n.collapse.region<-sapply(unique.regions,FUN=stocks.region,Y=collapse.stocks)
collapse.region.decade<-glm(cbind(n.collapse.region,n.stocks.region)~unique.regions,family=binomial)
anova(collapse.region.decade,test="Chisq")
summary(collapse.region.decade)

# Do finer regions, standard collapse threshold
collapse.threshold<-0.25
concatenate.fn<-function(X) paste(X[1],X[2],sep=" ")
data$REGIONEW<-apply(X=cbind(as.character(data$REGION),as.character(data$EW)),1,FUN=concatenate.fn)
region.id<-sapply(stock.names,extract.region,data[,c(1,16)])
unique.regions<-unique(region.id)
stocks.region<-function(X,Y) length(which(Y==X))
n.stocks.region<-sapply(unique.regions,FUN=stocks.region,Y=region.id)
collapse.stocks<-region.id[which(minB<=collapse.threshold)]
n.collapse.region<-sapply(unique.regions,FUN=stocks.region,Y=collapse.stocks)
collapse.region.glm<-glm(cbind(n.collapse.region,n.stocks.region)~unique.regions,family=binomial)
anova(collapse.region.glm,test="Chisq")
summary(collapse.region.glm)

