# TODO: Add estF and SD of estF.
# 
# Author: lsalas
###############################################################################


## This code takes as input a simulations file and specifications from a yaml on how to process it, and produces: 
#
#	Mean, SE, 5th and 95th for each site, averaged over all years of diffR and estR 
#	Mean, SE, 5th and 95th averaged over all sites and years of diffR and estR
#
# Here is what it does:
# For each year of the simulation it first determines which sites were surveyed and how many times for N or F. 
# For each survey simulation (i.e., each simulation, each year, set of sites, and survey type [N/F]), 
# the value of N or F is calculated at day = 0 (i.e., the top of the curvature) if there was more than 
# one survey to the site each year, or if surveyed only once, taken as the value at the day surveyed for 
# each site. Then we fitted a regression for the set of sites to estimate the count for each site in each 
# simulation. In the case of count of adults (N), it was count Adults = siteId + Year + Year^2, if the 
# sites were surveyed more than once, or without the quadrature if sampled only once. In the case of 
# counts of juveniles (F), we included the interaction of siteId and year to permit each site to have its 
# own productivity: count Juvs = siteId x Year + Year^2. Again, the quadrature was not used if the sites 
# were surveyed only once. After all the values were estimated for all simulations, we calculated R for 
# each simulation, and then summarized across simulations by obtaining means, standard errors, and 
# confidence intervals.

## Dependencies
library(yaml);library(plyr)


## Functions we will need

# getEstimateR: estimates R from a simulation
# dfN is the set of values of N from a simulation set
# dfF is the set of values of F from a simulation set
# trJuvVisits is the number of surveys conducted where juveniles/chicks were counted
# trAdVisits is the number of surveys conducted where adults were counted
getEstimateR<-function(dfN,dfF,trJuvVisits,trAdVisits){
	dfN$Year<-as.character(dfN$YearSampled)
	if(trAdVisits>1){		#trAdVisits==3
		fmlad<-as.formula("sample~siteId+Year+I(day^2)")
	}else{
		#dfN<-aggregate(sample~Year+siteId+areaMean+siteMean,dfN,max)
		fmlad<-as.formula("sample~siteId+Year")
	}
	nest<-lm(fmlad,dfN)
	predata<-unique(dfN[,c("Year","siteId","areaMean","siteMean")])
	predata$day<-0
	est<-predict(nest,predata,se=T)
	predata$estimateN<-est$fit
	predata$SEestimateN<-est$se.fit	#get rid of this
	
	dfF$Year<-as.character(dfF$YearSampled)
	if(trJuvVisits==3){
		fmljv<-as.formula("lgFsampled~siteId*Year+I(Fday^2)")
	}else{
		dfF<-aggregate(lgFsampled~Year+siteId+sampR+logF,dfF,max)
		fmljv<-as.formula("lgFsampled~siteId*Year")
	}
	fest<-lm(fmljv,dfF)
	prefdata<-unique(dfF[,c("Year","siteId","sampR","logF")])
	prefdata$Fday<-0
	estf<-predict(fest,prefdata,se=T)
	prefdata$estimateF<-estf$fit
	prefdata$SEestimateF<-estf$se.fit	#get rid of this too
	
	estdata<-merge(predata,prefdata,by=c("Year","siteId"))
	estdata$estR<-exp(estdata$estimateF)/exp(estdata$estimateN)
	return(estdata)
}

# validateInputs: validates the inputs of an analysis definition .yaml file against the data
# trSites is the number of sites surveyed in a simulation
# trYears is the number of years surveyed in a simulation
# trVisitsJv is the number of surveys conducted where juveniles/chicks were counted
# trVisitsAd is the number of surveys conducted where adults were counted
# data is the simulated data for adult surveys
# Fsimdata is the simulated data for juvenile/chick surveys
validateInputs<-function(trSites,trYears,trVisitsJv,trVisitsAd,data,Fsimdata){	#
	report<-""
	if(!trYears %in% c(2,3,5,6)){
		report<-paste(report,"The number of years specified in the analysis(",trYears,") is invalid - must be 2, 3 or 5; ",sep="")
	}
	if(NROW(unique(Fsimdata$siteId))<trSites){
		report<-paste(report,"The simulation file contains fewer sites (",NROW(unique(Fsimdata$siteId)),") than specified in this analysis; ",sep="")
	}
	if(trVisitsJv!= max(unique(Fsimdata$FsampId))){
		report<-paste(report,"The number of visits for juveniles in the simulation file (",max(unique(Fsimdata$FsampId)),") does not match what is specified in this analysis - must be the same number of visits as in the simulations.",sep="")
	}
	if(trVisitsAd!= max(unique(data$replicateId))){
		report<-paste(report,"The number of visits for juveniles in the simulation file (",max(unique(data$replicateId)),") does not match what is specified in this analysis - must be the same number of visits as in the simulations.",sep="")
	}
	
	return(report)
}

# validateTotalSites validates the total number of sites surveyd for juveniles in the simulation data
# tsites is the expected total number of sites
# Fsimdata is the simulated data for juvenile/chick surveys 
validateTotalSites<-function(tsites,Fsimdata){
	report<-""
	if(NROW(unique(Fsimdata$siteId))<tsites){
		report<-paste(report,"The simulation file contains fewer sites (",NROW(unique(Fsimdata$siteId)),") than specified in total, across all transects, in this analysis; ",sep="")
	}
	return(report)
}

# getResdf is a wrapper function that calls getEstimateR
# Ndf is the set of values of N from a simulation set
# Fdf is the set of values of F from a simulation set
# trJuvVisits is the number of surveys conducted where juveniles/chicks were counted
# trAdVisits is the number of surveys conducted where adults were counted
getResdf<-function(Ndf,Fdf,trJuvVisits,trAdVisits){
	resdf<-data.frame()
	for(ii in unique(Ndf$SetId)){
		dfN<-subset(Ndf,SetId==ii)
		dfF<-subset(Fdf,SetId==ii)
		estdf<-getEstimateR(dfN,dfF,trJuvVisits,trAdVisits)
		resdf<-rbind(resdf,estdf)
	}
	resdf$diffR<-resdf$estR-resdf$sampR
	resdf$absDiffR<-abs(resdf$diffR)
	return(resdf)
}

# getMeansQuants calculates the mean and standard deviation of results for the "var" parameter by site
# srdf is the simulated data to be summarized into mean and confidence interval values for each site
# var is the name of the variable to summarize
# snyrs is an auxiliary table with metadata on years surveyed
getMeansQuants<-function(srdf,var,snyrs){
	fml<-paste(var,"~siteId+TransectId",sep="");mvn<-paste("mean_",var,sep="");svn<-paste("SD_",var,sep="")
	evn<-paste("SE_",var,sep="");lvn<-paste("Lower05_",var,sep="");uvn<-paste("Upper05_",var,sep="")
	mdf<-aggregate(as.formula(fml),srdf,mean);names(mdf)<-c("siteId","TransectId",mvn)
	sdf<-aggregate(as.formula(fml),srdf,sd);names(sdf)<-c("siteId","TransectId",svn)
	tdf<-merge(mdf,sdf,by=c("siteId","TransectId"));tdf<-merge(tdf,snyrs,by=c("siteId","TransectId"))
	tdf<-tdf[,c("siteId","TransectId","NumYears",mvn,svn)]
	tdf[,evn]<-tdf[,svn]/sqrt(tdf$NumYears)
	tdf[,lvn]<-apply(tdf[,c(mvn,evn)],1,function(x){qnorm(p=0.05,mean=x[1],sd=x[2])})
	tdf[,uvn]<-apply(tdf[,c(mvn,evn)],1,function(x){qnorm(p=0.95,mean=x[1],sd=x[2])})
	tdf<-tdf[order(tdf$TransectId,tdf$siteId),]
	return(tdf)
}

# calculateEstsBySite calculates N and R for each site in the data
# srdf is the data to use for the calculations
# scenario is a variable passing the name of the analysis scenario being processed
calculateEstsBySite<-function(srdf,scenario){
	snyrs<-aggregate(Year~siteId+TransectId,srdf,FUN=function(x){y<-NROW(unique(x));return(y)});names(snyrs)<-c("siteId","TransectId","NumYears")
	smnN<-aggregate(estimateN~siteId+TransectId,srdf,mean);names(smnN)<-c("siteId","TransectId","lgmean_estN")
	smnNsd<-aggregate(estimateN~siteId+TransectId,srdf,sd);names(smnNsd)<-c("siteId","TransectId","SD_lgestN")
	smnN<-merge(smnN,smnNsd,by=c("siteId","TransectId"))
	smnN$mean_estN<-exp(smnN$lgmean_estN)
	
	sestR<-getMeansQuants(srdf=srdf,var="estR",snyrs=snyrs)
	sdifR<-getMeansQuants(srdf=srdf,var="diffR",snyrs=snyrs)
	sabdR<-getMeansQuants(srdf=srdf,var="absDiffR",snyrs=snyrs)
	sestF<-getMeansQuants(srdf=srdf,var="estimateF",snyrs=snyrs)
	
	bysitedf<-merge(sestR,sdifR,by=c("siteId","TransectId","NumYears"))
	bysitedf<-merge(bysitedf,sabdR,by=c("siteId","TransectId","NumYears"))
	bysitedf<-merge(bysitedf,sestF,by=c("siteId","TransectId","NumYears"))
	bysitedf<-merge(bysitedf,smnN[,c("siteId","TransectId","mean_estN","SD_lgestN")],by=c("siteId","TransectId"))
	bysitedf<-bysitedf[order(bysitedf$TransectId,bysitedf$siteId),]
	bysitedf$Scenario<-scenario
	return(bysitedf)
}

# calculateEstsOverall summarizes results for the entire set of sites
# bysitedf is the file with estimates by site
# scenario is a variable passing the name of the analysis scenario being processed
calculateEstsOverall<-function(bysitedf,scenario){
	odf<-bysitedf[,c("mean_estR","mean_diffR","mean_absDiffR","mean_estN","SE_estR","SE_diffR","SE_absDiffR")]
	odf$wmeR<-odf$mean_estR*sqrt(odf$mean_estN);odf$wdiR<-odf$mean_diffR*sqrt(odf$mean_estN);odf$wabR<-odf$mean_absDiffR*sqrt(odf$mean_estN)
	odf$wSE_estR<-odf$SE_estR*sqrt(odf$mean_estN);odf$wSE_diffR<-odf$SE_diffR*sqrt(odf$mean_estN);odf$wSE_absDiffR<-odf$SE_absDiffR*sqrt(odf$mean_estN)
	sumW<-sum(sqrt(odf$mean_estN))
	wom<-adply(.data=odf[,c("wmeR","wdiR","wabR")],.margins=2,.fun=function(x,sumW){
				y<-sum(x)/sumW;
				nv<-names(x);
				tdf<-data.frame(weightedMean=y);
				return(tdf)},sumW=sumW)
	wom$Parameter<-c("EstR","DiffR","AbsDiffR")
	wom$Mean<-c(mean(odf$mean_estR),mean(odf$mean_diffR),mean(odf$mean_absDiffR))
	wom$SD<-c(mean(bysitedf$SD_estR),mean(bysitedf$SD_diffR),mean(bysitedf$SD_absDiffR))
	wom$SE<-c(mean(bysitedf$SE_estR)/sqrt(nrow(bysitedf)),mean(bysitedf$SE_diffR)/sqrt(nrow(bysitedf)),mean(bysitedf$SE_absDiffR)/sqrt(nrow(bysitedf)))
	wom$wghtSE<-c(sum(odf$wSE_estR)/(sumW*sqrt(nrow(odf))),sum(odf$wSE_diffR)/(sumW*sqrt(nrow(odf))),sum(odf$wSE_absDiffR)/(sumW*sqrt(nrow(odf))))
	wom<-wom[,c("Parameter","Mean","weightedMean","SD","SE","wghtSE")]
	wom$Lower05<-unlist(lapply(X=1:3,FUN=function(x,wom){y<-qnorm(p=0.05,mean=wom[x,"Mean"],sd=wom[x,"SE"]);return(y)},wom=wom))
	wom$Upper95<-unlist(lapply(X=1:3,FUN=function(x,wom){y<-qnorm(p=0.95,mean=wom[x,"Mean"],sd=wom[x,"SE"]);return(y)},wom=wom))
	wom$wghtLower05<-unlist(lapply(X=1:3,FUN=function(x,wom){y<-qnorm(p=0.05,mean=wom[x,"weightedMean"],sd=wom[x,"wghtSE"]);return(y)},wom=wom))
	wom$wghtUpper95<-unlist(lapply(X=1:3,FUN=function(x,wom){y<-qnorm(p=0.95,mean=wom[x,"weightedMean"],sd=wom[x,"wghtSE"]);return(y)},wom=wom))
	wom$Scenario<-scenario
	return(wom)
}

# calculateStats is a wrapper function for calculateEstsBySite
# srdf is the data to use for the calculations
# scenario is a variable passing the name of the analysis scenario being processed
# single site is a flag indicating if statistics should be calculated for a single site or overall
calculateStats<-function(srdf,scenario,singlesite=FALSE){
	reslst<-list()
	
	bysitedf<-calculateEstsBySite(srdf=srdf,scenario=scenario)
	if(singlesite==TRUE){
		bysitedf<-subset(bysitedf,siteId=="Site1")	
	}
	
	reslst[["bysite"]]<-bysitedf
	
	##################################
	
	wom<-calculateEstsOverall(bysitedf=bysitedf,scenario=scenario)
	reslst[["overall"]]<-wom
	
	return(reslst)
}


## Inputs:
defyamlpth<-"C:/Users/lsalas/git/sparklemotion/AWPE/AnalysisDefinitions/"
anaDefs<-list.files(defyamlpth,pattern="6Y")


## Run the analyses by analysis definition...
for(adf in anaDefs){
	#read the analysis definition
	andefs<-try(yaml.load_file(paste(defyamlpth,adf,sep="")),silent=FALSE)
	sfpth<-andefs$SimulationsFile
	print(sfpth)
	load(sfpth)	#loads data, Fsimdata
	
	nsites<-c(4,10,18)
	nyears<-c(2,3,5,6)
	ntrans<-NROW(andefs$Transects)
	
	if(ntrans==1){
		for(yy in nyears){
			
			andefs$Transects[[1]]$NumYears<-yy
			for(ss in nsites){
				andefs$Transects[[1]]$NumSites<-as.integer(ss)
				vimp<-ldply(.data=c(1:NROW(andefs$Transects)),.fun=function(x,andefs,data,Fsimdata){
							trId<-andefs$Transects[[x]]$TransectId;
							trSites<-as.integer(andefs$Transects[[x]]$NumSites);
							trYears<-as.integer(andefs$Transects[[x]]$NumYears);
							trVisitsJv<-as.integer(andefs$Transects[[x]]$NumVisitsJuvenileSurveys);
							trVisitsAd<-as.integer(andefs$Transects[[x]]$NumVisitsAdultSurveys);
							val<-validateInputs(trSites=trSites,trYears=trYears,trVisitsJv=trVisitsJv,trVisitsAd=trVisitsAd,data=data,Fsimdata=Fsimdata);
							tdf<-data.frame(transectId=x,report=val,numSites=trSites);
							return(tdf)
						},andefs=andefs,data=data,Fsimdata=Fsimdata)
				if(paste(vimp$report,collapse="")!="")stop(paste(vimp$report,collapse=""))
				val<-validateTotalSites(tsites=sum(vimp$numSites),Fsimdata=Fsimdata)
				eds<-cumsum(vimp$numSites)
				sts<-1
				if(NROW(eds)>1){
					for(e in 1:(NROW(eds)-1)){sts<-c(sts,eds[e]+1)}
				}
				if(ss==4){ #Randomize the sites selected
					randsts<-sample(0:14,1);sts<-sts+randsts;eds<-eds+randsts
				}else if(ss==10){
					randsts<-sample(0:8,1);sts<-sts+randsts;eds<-eds+randsts
				}
				vimp$startSite<-sts;vimp$endSite<-eds			
				simresdf<-ldply(.data=c(1:NROW(vimp)), .fun=function(x,vimp,andefs,data,Fsimdata){
							tt<-vimp$transectId[x]
							trId<-andefs$Transects[[tt]]$TransectId
							trSites<-as.integer(andefs$Transects[[tt]]$NumSites)
							trYears<-as.integer(andefs$Transects[[tt]]$NumYears)
							trJuvVisits<-as.integer(andefs$Transects[[tt]]$NumVisitsJuvenileSurveys)
							trAdVisits<-as.integer(andefs$Transects[[tt]]$NumVisitsAdultSurveys)
							
							stsite<-vimp$startSite[x];endsite<-vimp$endSite[x]
							seqsites<-c(stsite:endsite);selsites<-paste("Site",seqsites,sep="")
							if(trYears==2){
								selyears<-c(1,4)
							}else if(trYears==3){
								selyears<-c(1,3,5)
							}else if(trYears==5){
								selyears<-c(1:5)
							}else{
								selyears<-c(0:5)
							}
							
							Ndf<-subset(data,(siteId %in% selsites) & (YearSampled %in% selyears))
							Fdf<-subset(Fsimdata,(siteId %in% selsites) & (YearSampled %in% selyears))
							
							resdf<-getResdf(Ndf,Fdf,trJuvVisits,trAdVisits)
							resdf$TransectId<-trId
							return(resdf)
						},vimp=vimp,andefs=andefs,data=data,Fsimdata=Fsimdata)	
				
				nvis<-ifelse(grepl("3visits",andefs$SimulationName),"3visits",ifelse(grepl("2visits",andefs$SimulationName),"2visits",
								ifelse(grepl("1visit",andefs$SimulationName),"1visit",ifelse(grepl("N3F1",andefs$SimulationName),"N3F1visits",
										ifelse(grepl("N1F3",andefs$SimulationName),"N1F3visits",ifelse(grepl("N2F1",andefs$SimulationName),"N2F1visits","N1F2visits"))))))
				nwks<-ifelse(grepl("2weeks",andefs$SimulationName),"2weeks","4weeks")
				scenario<-paste("AWPE_",yy,"years",nvis,nwks,ss,"sites",sep="")
				
				results<-calculateStats(srdf=simresdf,scenario=scenario)
				results$Scenario<-scenario
				save(simresdf,results,file=paste("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/IandMR8/AWPE/results/",scenario,".RData",sep=""))
				
				if(ss==18){
					scenario<-paste("AWPE_",yy,"years",nvis,nwks,"1site",sep="")
					results<-calculateStats(srdf=simresdf,scenario=scenario,singlesite=TRUE)
					results$Scenario<-scenario
					save(simresdf,results,file=paste("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/IandMR8/AWPE/results/",scenario,".RData",sep=""))
				}
				
			}
		}
			
	}else{
		
		#Each transect definition has its own number of years to process, per numYears
		
		vimp<-ldply(.data=c(1:NROW(andefs$Transects)),.fun=function(x,andefs,data,Fsimdata){
					trId<-andefs$Transects[[x]]$TransectId;
					trSites<-as.integer(andefs$Transects[[x]]$NumSites);
					trYears<-as.integer(andefs$Transects[[x]]$NumYears);
					trVisitsJv<-as.integer(andefs$Transects[[x]]$NumVisitsJuvenileSurveys);
					trVisitsAd<-as.integer(andefs$Transects[[x]]$NumVisitsAdultSurveys);
					val<-validateInputs(trSites=trSites,trYears=trYears,trVisitsJv=trVisitsJv,trVisitsAd=trVisitsAd,data=data,Fsimdata=Fsimdata);
					tdf<-data.frame(transectId=x,report=val,numSites=trSites);
					return(tdf)
				},andefs=andefs,data=data,Fsimdata=Fsimdata)
		if(paste(vimp$report,collapse="")!="")stop(paste(vimp$report,collapse=""))
		val<-validateTotalSites(tsites=sum(vimp$numSites),Fsimdata=Fsimdata)
		eds<-cumsum(vimp$numSites)
		sts<-1
		if(NROW(eds)>1){
			for(e in 1:(NROW(eds)-1)){sts<-c(sts,eds[e]+1)}
		}
		vimp$startSite<-sts;vimp$endSite<-eds			
		simresdf<-ldply(.data=c(1:NROW(vimp)), .fun=function(x,vimp,andefs,data,Fsimdata){
					tt<-vimp$transectId[x]
					trId<-andefs$Transects[[tt]]$TransectId
					trSites<-as.integer(andefs$Transects[[tt]]$NumSites)
					trYears<-as.integer(andefs$Transects[[tt]]$NumYears)
					trJuvVisits<-as.integer(andefs$Transects[[tt]]$NumVisitsJuvenileSurveys)
					trAdVisits<-as.integer(andefs$Transects[[tt]]$NumVisitsAdultSurveys)
					
					stsite<-vimp$startSite[x];endsite<-vimp$endSite[x]
					seqsites<-c(stsite:endsite);selsites<-paste("Site",seqsites,sep="")
					if(trYears==2){
						selyears<-c(1,4)
					}else if(trYears==3){
						selyears<-c(1,3,5)
					}else if(trYears==5){
						selyears<-c(1:5)
					}else{
						selyears<-c(0:5)
					}
					
					Ndf<-subset(data,(siteId %in% selsites) & (YearSampled %in% selyears))
					Fdf<-subset(Fsimdata,(siteId %in% selsites) & (YearSampled %in% selyears))
					
					resdf<-getResdf(Ndf,Fdf,trJuvVisits,trAdVisits)
					resdf$TransectId<-trId
					return(resdf)
				},vimp=vimp,andefs=andefs,data=data,Fsimdata=Fsimdata)	
		
		scenario<-andefs$SimulationName
		results<-calculateStats(srdf=simresdf,scenario=scenario)
		results$Scenario<-scenario
		save(simresdf,results,file=paste("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/IandMR8/AWPE/results/",scenario,".RData",sep=""))
		
	}
	
}

q(save="no")




