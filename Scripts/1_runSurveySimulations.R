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

## samplesFromCurve: This function takes F samples from a peak value of F (provided) and a curvature
# F is the maximum value of fecundity (at the peak of the curvature)
# curv is the quadrature of the curve (always negative)
# seasonStart and seasonEnd are the number of days before and after the peak defining the breeding season. 
#   So, seasonStart is always negative, and seasonEnd always positive. (Default: -30, 30)
# numsamples is the number of samples to take in the season: a value 1-3 (Default: 3)
# sampInt is the minimum interval between samples. (Default: 15)
# sdFerr is the error in estimating F in a sample, a standard deviation around the F value obtained from the curve. (Default: 0.2)
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

scenarios<-c("1visit","2visits2weeks","2visits4weeks","3visits2weeks","N1F2visits2weeks","N1F3visits2weeks","N2F1visits2weeks","N3F1visits2weeks")
for(ss in scenarios){
	print(ss)
	ssa<-ifelse(grepl("N1",ss),"1visit",ifelse(grepl("N2",ss),"2visits2weeks",ifelse(grepl("N3",ss),"3visits2weeks",ss)))
	yamlpth<-paste("C:/Users/lsalas/git/sparklemotion/AWPE/SimulationDefinitions/N",ssa,".yaml",sep="")
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
	#ssj="2visits2weeks"
	ssj<-ifelse(grepl("F1",ss),"1visit",ifelse(grepl("F2",ss),"2visits2weeks",ifelse(grepl("F3",ss),"3visits2weeks",ss)))
	fsimdefs<-yaml.load_file(paste("C:/Users/lsalas/git/sparklemotion/AWPE/SimulationDefinitions/F",ssj,".yaml",sep=""))
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
	
	# Calculate the true F from the sampled R and trueN, then logF
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
	save(data,fdata,Fsimdata,file=paste(filepth,nboot,"sim_",ss,".RData",sep=""))
}

q(save="no")



