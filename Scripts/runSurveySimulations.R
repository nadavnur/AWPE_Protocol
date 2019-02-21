# TODO: Add comment
# 
# Author: lsalas
###############################################################################


#source all needed files
#create a desired configuration of the yaml control file
#read the yaml
#execute

#dependencies
library(yaml); library(plyr)

# Run the function to get the F samples
samplesFromCurve<-function(F,curv,seasonStart=-30,seasonEnd=30,numsamples=3,sampInt=15,sdFerr=0.2){
	if(numsamples==1){
		sampdays<-round(runif(1,min=seasonStart,max=seasonEnd))
	}else{
		firstsampend<-seasonEnd-((numsamples-1)*sampInt) - 1
		firstsamp<-round(runif(1,min=seasonStart,max=firstsampend))
		sampdays<-firstsamp
		for(ss in 1:(numsamples-1)){
			sampdays<-c(sampdays,firstsamp+(ss*sampInt))
		}
	}
	df<-data.frame(day=sampdays)
	df$lgFday<-(curv*df$day^2)+F
	df$lgFday<-ifelse(df$lgFday<0,0,df$lgFday)	#this should never happen
	df$lgFsampled<-apply(df,1,FUN=function(x,sdv){y<-ifelse(x[2]==0,0,
						rnorm(n=1,mean=x[2],sd=sdv));return(y)},sdv=sdFerr)
	df$lgFsampled<-ifelse(df$lgFsampled<0,0,df$lgFsampled)
	return(df)
}

classpth<-"C:/Users/lsalas/git/sparklemotion/PowerAnalysis/General/Classes/"

source("C:/Users/lsalas/git/sparklemotion/PowerAnalysis/General/getYamlContents.R")
source(paste(classpth,"SurveySiteClassDefinition.R",sep=""))
source(paste(classpth,"TransectClassDefinition.R",sep=""))
source(paste(classpth,"StudyAreaClassDefinition.R",sep=""))
source(paste(classpth,"PopulationClassDefinition.R",sep=""))
source(paste(classpth,"PopulationSetClassDefinition.R",sep=""))
source(paste(classpth,"SimulationClassDefinition.R",sep=""))
source(paste(classpth,"SimulationDispatchClassDefinition.R",sep=""))

nboot<-1000

scenarios<-c("1visit","2visits2weeks","2visits4weeks","3visits2weeks")
#for(ss in scenarios){
	ss="N1F2"
	ssa="1visit"
	print(ss)
	yamlpth<-paste("C:/Users/lsalas/git/sparklemotion/AWPE/SimDefinitions/N",ssa,".yaml",sep="")
	simdefs<-yaml.load_file(yamlpth)
	simdefnam<-simdefs$SimulationName
		
	test<-new("SimulationDispatch")
	SimulationId(test)<-simdefnam
	SimulationName(test)<-"AWPE simulation of N 1 sample"
	YamlFilePath(test)<-yamlpth
	NumSimulations(test)<-as.integer(nboot)
	VaryingParameter(test)<-"MuTrend"
	VaryingValues(test)<-c(0)		#   c(0,0)
	test<-getSynthData(test)
	data<-SimulatedData(test)
	
	##################################
	## Now the Juveniles
	ssj="2visits2weeks"
	fsimdefs<-yaml.load_file(paste("C:/Users/lsalas/git/sparklemotion/AWPE/SimDefinitions/F",ssj,".yaml",sep=""))
	ftvars<-c("SetId","SimulationId","PopulationId","YearSampled","areaMean","TransectId","siteId","siteMean")
	fdata<-unique(data[,ftvars])
	fdata$trueN<-exp(fdata$siteMean)
	fdn<-names(fdata)
	
	muR<-fsimdefs$Rmu;sdR<-fsimdefs$Rsd
	
	if(fsimdefs$FixSiteEffects==FALSE){
		linkflds<-c("siteId","YearSampled")
	}else{
		linkflds<-c("SetId","siteId","YearSampled")
	}
	siteyears<-unique(fdata[,linkflds])
	rsv<-rnorm(nrow(siteyears),muR,sdR) 
	siteyears$sampR<-rsv
	siteyears$sampR<-ifelse(siteyears$sampR<0.05,0.05,siteyears$sampR)
	fdata<-merge(fdata,siteyears,by=c(linkflds))
	fdata<-fdata[,c(fdn,"sampR")]
	
	# Calculate the true F from the sampled R, then logF
	fdata$trueF<-fdata$trueN*fdata$sampR
	fdata$logF<-log(fdata$trueF)
	linkvars<-c("SetId","SimulationId","PopulationId","YearSampled","TransectId","siteId","logF")
	
	#Now loop through the populations
	Fsimdata<-data.frame()
	numpops<-NROW(fsimdefs$Populations)
	for(pp in 1:numpops){
		popId<-fsimdefs$Populations[[pp]]$PopulationId
		curv<-fsimdefs$Populations[[1]]$Cseason
		seastart<-fsimdefs$Populations[[1]]$SeasonStart
		seaend<-fsimdefs$Populations[[1]]$SeasonEnd
		numsamples<-fsimdefs$Populations[[1]]$MaxReplicates
		repspan<-fsimdefs$Populations[[1]]$ReplicateSpan
		pdet<-fsimdefs$Populations[[1]]$ProbDetection	#not using this for now...
		samperr<-fsimdefs$Populations[[1]]$SigSample
		sampFdf<-adply(.data=fdata[,linkvars],.margins=1,.fun=function(x,curv,numsamples,seasonStart,seasonEnd,sampInt,sdFerr){
					F<-x$logF;
					fsv<-samplesFromCurve(F=F,curv=curv,numsamples=numsamples,seasonStart=seasonStart,seasonEnd=seasonEnd,sampInt=sampInt,sdFerr=sdFerr);
					rdf<-data.frame();
					for(nn in 1:numsamples){
						rdf<-rbind(rdf,x);
					}
					rdf$FpopulationId<-popId;rdf$FsampId<-1:numsamples;rdf$Fday<-fsv$day;rdf$lgFday<-fsv$lgFday;rdf$lgFsampled<-fsv$lgFsampled;
					return(rdf);
				},curv=curv,numsamples=numsamples,seasonStart=seastart,seasonEnd=seaend,sampInt=repspan,sdFerr=samperr)
		
		fsd<-merge(fdata,sampFdf,by=linkvars,all.x=T)
		fsd<-fsd[,c("SetId","SimulationId","PopulationId","YearSampled","TransectId","siteId","areaMean","siteMean","trueN","sampR","trueF","logF","FpopulationId","FsampId","Fday","lgFday","lgFsampled")]
		Fsimdata<-rbind(Fsimdata,fsd)
	}
	
	###########################################
	#save
	filepth<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/IandMR8/AWPE/simulations/"
	#save(data,fdata,Fsimdata,file=paste(filepth,nboot,"sim_",substr(simdefnam,3,nchar(simdefnam)),".RData",sep=""))
	save(data,fdata,Fsimdata,file=paste(filepth,nboot,"sim_",ss,".RData",sep=""))
#}

q(save="no")



#test...
library(ggplot2)
ds1d<-subset(data,SetId=="AWPEsurveyDesign:AWPEtest:1" & SimulationId=="MuTrend:0")
#it's the same mean for all areas...
p<-ggplot(ds1d,aes(x=YearSampled,y=areaMean)) + geom_line(aes(color=AreaId))
print(p)
#each area gets its top mean (day 0 mean)
p<-ggplot(ds1d,aes(x=YearSampled,y=siteMean)) + geom_line(aes(color=AreaId))
#each area gets its sampled values - from which we'd construct curvatures
p<-ggplot(ds1d,aes(x=YearSampled,y=siteMean)) + geom_line(aes(color=AreaId)) + geom_point(aes(y=sample,color=AreaId))


####################



library(ggplot2)
#plot in one graph the mean of means across simulations and sites of true N, F, and R
trueNsimsite<-aggregate(siteMean~siteId+Year,resdf,mean)
trueFsimsite<-aggregate(logF~siteId+Year,resdf,mean)
truesimsite<-merge(trueNsimsite,trueFsimsite,by=c("siteId","Year"))
truesimsite$realR<-exp(truesimsite$logF)/exp(truesimsite$siteMean)
truesiteN<-aggregate(siteMean~Year,truesimsite,mean)
truesiteF<-aggregate(logF~Year,truesimsite,mean)
truesiteR<-aggregate(realR~Year,truesimsite,mean)
truesite<-merge(truesiteN,truesiteF,by="Year")
truesite<-merge(truesite,truesiteR,by="Year")



Nsimsite<-aggregate(estimateN~siteId+Year,resdf,mean)
Fsimsite<-aggregate(estimateF~siteId+Year,resdf,mean)
Rsimsite<-aggregate(estR~siteId+Year,resdf,mean)
simsite<-merge(Nsimsite,Fsimsite,by=c("siteId","Year"))
simsite<-merge(simsite,Rsimsite,by=c("siteId","Year"))
simsite$estRb<-exp(simsite$estimateF)/exp(simsite$estimateN)
siteN<-aggregate(estimateN~Year,simsite,mean)
siteF<-aggregate(estimateF~Year,simsite,mean)
siteR<-aggregate(estR~Year,simsite,mean)
siteRb<-aggregate(estRb~Year,simsite,mean)
aggsite<-merge(siteN,siteF,by="Year")
aggsite<-merge(aggsite,siteR,by="Year")
aggsite<-merge(aggsite,siteRb,by="Year")
#plot in 18 facets the mean across simulations of true N, F, and R

p<-ggplot(data=simsite,aes(x=as.integer(Year),y=exp(estimateF))) + geom_line(aes(color=siteId))

