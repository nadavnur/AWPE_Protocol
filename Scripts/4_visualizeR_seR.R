# TODO: Add comment
# 
# Author: lsalas
###############################################################################

library(plyr);library(ggplot2)
pth<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/IandMR8/AWPE/results/"

################
# get the empirical estimates of confidence intervals
getFvisitsData<-function(pattern,pth){
	fls<-list.files(pth,pattern=pattern)
	Fvisit<-data.frame()
	for(ff in fls){
		if(pattern=="1visit"){
			wv<-"W2";fsv<-"F1";nsv<-"N1"
		}else if(pattern=="2visit"){
			wc<-as.numeric(regexpr("week",ff))
			wv<-paste("W",substr(ff,wc-1,wc-1),sep="")
			fsv<-"F2";nsv<-"N2"
		}else if(pattern=="3visit"){
			wv<-"W2";fsv<-"F3";nsv<-"N3"
		}else{}
		if(substr(ff,12,12)=="N"){
			nsv<-substr(ff,12,13)
		}
		if(substr(ff,14,14)=="F"){
			fsv<-substr(ff,14,15)
		}
		yrv<-as.numeric(substr(ff,6,6))
		sc<-as.numeric(regexpr("site",ff))
		svt<-substr(ff,sc-2,sc-1)
		if(svt=="10"){
			sv<-10
		}else if(svt=="18"){
			sv<-18
		}else{
			sv<-as.numeric(substr(svt,2,2))
		}
		load(paste(pth,ff,sep=""))
		if(sv==1){
			simresdf<-subset(simresdf,siteId=="Site1")
		}
		simresdf$sim<-paste("Y",simresdf$Year,"_",simresdf$siteId,sep="")
		tdf<-data.frame()
		for(ss in unique(simresdf$sim)){
			df<-subset(simresdf,sim==ss)
			df$simCount<-1:(nrow(df))
			tdf<-rbind(tdf,df)
		}
		
		meanEstR<-aggregate(estR~simCount,tdf,mean)
		qvals<-as.numeric(quantile(meanEstR$estR,probs=c(0.05,0.95)))
		
		resdf<-data.frame(years=paste("Y",yrv,sep=""),sites=paste("S",sv,sep=""),Nsample=nsv,prob05=qvals[1],prob95=qvals[2],weeks=wv,Fsample=fsv,Nrep=paste("R",yrv*sv,sep=""))
		Fvisit<-rbind(Fvisit,resdf)
	}
	return(Fvisit)
	
}

empiricaldf<-ldply(.data=c("1visit","2visit","3visit"),.fun=getFvisitsData,pth=pth)
empiricaldf<-subset(empiricaldf,years!="YNA")
empiricaldf$scenario<-paste0(empiricaldf$years,empiricaldf$sites,empiricaldf$Nrep,empiricaldf$Nsample,empiricaldf$Fsample,empiricaldf$weeks)
empiricaldf<-empiricaldf[,c("years","sites","Nrep","weeks","Nsample","Fsample","scenario","prob05","prob95")]
names(empiricaldf)<-gsub("prob","emp_prob",names(empiricaldf))
empiricaldf$coverage<-empiricaldf$emp_prob95-empiricaldf$emp_prob05


#################
#calculate the theoretical estimates of confidence intervals
respth<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/IandMR8/AWPE/results/"
fls<-list.files(respth,pattern="AWPE")
theoreticaldf<-data.frame()
for(ff in fls){
	load(file=paste(respth,ff,sep=""))
	bsdf<-results[[1]];ovdf<-results[[2]]
	
	theoreticaldf<-rbind(theoreticaldf,ovdf)
}

#add numYears, numSites, numVisitsAdult, numVisitsJuveniles, numWeeks
theoreticaldf$numYears<-ifelse(grepl("2years",theoreticaldf$Scenario),2,ifelse(grepl("3years",theoreticaldf$Scenario),3,ifelse(grepl("5years",theoreticaldf$Scenario),5,6)))
theoreticaldf$numSites<-ifelse(grepl("1site",theoreticaldf$Scenario),1,ifelse(grepl("4site",theoreticaldf$Scenario),4,ifelse(grepl("10site",theoreticaldf$Scenario),10,18)))
theoreticaldf$numWeeks<-ifelse(grepl("2weeks",theoreticaldf$Scenario),2,4)
theoreticaldf$numVisitsAd<-ifelse(grepl("N1F",theoreticaldf$Scenario),1,ifelse(grepl("N3F1",theoreticaldf$Scenario),3,ifelse(grepl("N2F1",theoreticaldf$Scenario),2,
				ifelse(grepl("3visit",theoreticaldf$Scenario),3,ifelse(grepl("2visit",theoreticaldf$Scenario),2,1)))))
theoreticaldf$numVisitsJv<-ifelse(grepl("N1F3",theoreticaldf$Scenario),3,ifelse(grepl("N3F1",theoreticaldf$Scenario),1,ifelse(grepl("N1F2",theoreticaldf$Scenario),2,
				ifelse(grepl("3visit",theoreticaldf$Scenario),3,ifelse(grepl("2visit",theoreticaldf$Scenario),2,1)))))
theoreticaldf$Nmodel<-"Additive"


scenvals<-lapply(1:nrow(theoreticaldf),function(x,odf){
			yrv<-odf[x,"numYears"];stv<-theoreticaldf[x,"numSites"];repv<-yrv*stv;wkv<-odf[x,"numWeeks"]
			scen<-odf[x,"Scenario"]
			if(grepl("1visit",scen)){
				fsv<-"F1";nsv<-"N1"
			}else if(grepl("2visit",scen)){
				fsv<-"F2";nsv<-"N2"
			}else if(grepl("3visit",scen)){
				fsv<-"F3";nsv<-"N3"
			}else{}
			if(substr(scen,12,12)=="N"){
				nsv<-substr(scen,12,13)
			}
			if(substr(scen,14,14)=="F"){
				fsv<-substr(scen,14,15)
			}
			sceval<-paste0("Y",yrv,"S",stv,"R",repv,nsv,fsv,"W",wkv)
			return(sceval)
		},odf=theoreticaldf)

theoreticaldf$scenario<-unlist(scenvals)
theoreticaldf<-subset(theoreticaldf,Parameter=="EstR")
theoreticaldf<-subset(theoreticaldf,Scenario!="AWPE_T6T2years1visits2weeks18sites")
theoreticaldf<-theoreticaldf[,c("scenario","Lower05","Upper95")]
names(theoreticaldf)<-c("scenario","the_prob05","the_prob95")
theoreticaldf$thcoverage<-theoreticaldf$the_prob95-theoreticaldf$the_prob05
#######################

## Both tables should be of the same size...
nrow(empiricaldf)==nrow(theoreticaldf)

## Merging and calculating disparities
sedf<-merge(empiricaldf,theoreticaldf,by="scenario")
sedf$mismatch<-sedf$the_prob05-sedf$emp_prob05

#Prepare to plot these
plotdf<-sedf[,c("scenario","emp_prob05","emp_prob95","Nsample","Fsample","Nrep","coverage","mismatch")]
plotdf$IntervalType<-"Empirical"
tdf<-sedf[,c("scenario","the_prob05","the_prob95","Nsample","Fsample","Nrep","thcoverage","mismatch")]
tdf$IntervalType<-"Theoretical"
names(tdf)<-names(plotdf)
plotdf<-rbind(plotdf,tdf)
names(plotdf)<-gsub("emp_prob","prob",names(plotdf))
plotdf$Nrep<-as.integer(substr(plotdf$Nrep,2,4))

####################################
## 
## Visual 1
# When num replicates is low, confidence interval is large, and empirical interval is larger than theoretical
pt1<-ggplot(plotdf,aes(x=scenario,ymin=prob05,ymax=prob95)) + geom_errorbar(aes(color=IntervalType),position="dodge") + coord_flip() + facet_wrap(~Fsample,ncol=3,scales="free")

## Visual 2
# Compare coverage (span of inner 95th percentile of values) vs number of samples (i.e., numYears * numSites) - note lower values (i.e., bias) in theoretical estimates at lower sample sizes
pt2a<-ggplot(plotdf,aes(x=Nrep,y=coverage)) + geom_point(aes(color=IntervalType)) + facet_grid(Fsample~Nsample) + labs(x="Number of samples (years x sites)", y="Coverage (span of inner 95th percentile of values)",color="Interval estimation")

# Bias is evident if plotted as ratio of coverage: empirical/theoretical
pdf<-sedf[,c("scenario","Nrep","thcoverage","coverage","Nsample","Fsample")]
pdf$coverage_ratio<-pdf$coverage/pdf$thcoverage
pdf$Nrep<-as.integer(substr(pdf$Nrep,2,4))

pt2b<-ggplot(pdf,aes(x=Nrep,y=coverage_ratio)) + geom_point() + facet_grid(Fsample~Nsample) + geom_hline(yintercept=1)

## Visual 3
## We model it and prove it's significant
pdf$lgNrep<-log(pdf$Nrep)
g<-lm(coverage_ratio~lgNrep*Nsample+lgNrep*Fsample,pdf)
summary(g)
pdf$predicted<-predict(g)

pt3<-ggplot(pdf,aes(x=Nrep,y=coverage_ratio)) + geom_point() + geom_line(aes(y=predicted),color="red") + facet_grid(Fsample~Nsample) + geom_hline(yintercept=1)

