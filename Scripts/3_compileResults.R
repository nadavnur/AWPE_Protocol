# TODO: Add comment
# 
# Author: lsalas
###############################################################################


###############
# The code below compiles all the data from the processing of simulations into csv files
###############

respth<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/IandMR8/AWPE/results/"
fls<-list.files(respth,pattern="AWPE")
bysitedf<-data.frame()
overalldf<-data.frame()
for(ff in fls){
	load(file=paste(respth,ff,sep=""))
	bsdf<-results[[1]];ovdf<-results[[2]]
	bysitedf<-rbind(bysitedf,bsdf)
	overalldf<-rbind(overalldf,ovdf)
}

#add numYears, numSites, numVisitsAdult, numVisitsJuveniles, numWeeks
overalldf$numYears<-ifelse(grepl("2years",overalldf$Scenario),2,ifelse(grepl("3years",overalldf$Scenario),3,ifelse(grepl("5years",overalldf$Scenario),5,6)))
overalldf$numSites<-ifelse(grepl("1site",overalldf$Scenario),1,ifelse(grepl("4site",overalldf$Scenario),4,ifelse(grepl("10site",overalldf$Scenario),10,18)))
overalldf$numWeeks<-ifelse(grepl("2weeks",overalldf$Scenario),2,4)
overalldf$numVisitsAd<-ifelse(grepl("N1F",overalldf$Scenario),1,ifelse(grepl("N3F1",overalldf$Scenario),3,ifelse(grepl("N2F1",overalldf$Scenario),2,
						ifelse(grepl("3visit",overalldf$Scenario),3,ifelse(grepl("2visit",overalldf$Scenario),2,1)))))
overalldf$numVisitsJv<-ifelse(grepl("N1F3",overalldf$Scenario),3,ifelse(grepl("N3F1",overalldf$Scenario),1,ifelse(grepl("N1F2",overalldf$Scenario),2,
						ifelse(grepl("3visit",overalldf$Scenario),3,ifelse(grepl("2visit",overalldf$Scenario),2,1)))))

bysitedf$numSites<-ifelse(grepl("1site",bysitedf$Scenario),1,ifelse(grepl("4site",bysitedf$Scenario),4,ifelse(grepl("10site",bysitedf$Scenario),10,18)))
bysitedf$numWeeks<-ifelse(grepl("2weeks",bysitedf$Scenario),2,4)
bysitedf$numVisitsAd<-ifelse(grepl("N1F3",bysitedf$Scenario),1,ifelse(grepl("N3F1",bysitedf$Scenario),3,
				ifelse(grepl("3visit",bysitedf$Scenario),3,ifelse(grepl("2visit",bysitedf$Scenario),2,1))))
bysitedf$numVisitsJv<-ifelse(grepl("N1F3",bysitedf$Scenario),3,ifelse(grepl("N3F1",bysitedf$Scenario),1,
				ifelse(grepl("3visit",bysitedf$Scenario),3,ifelse(grepl("2visit",bysitedf$Scenario),2,1))))

overalldf$Nmodel<-"Additive"
bysitedf$Nmodel<-"Additive"

#######################
scenvals<-lapply(1:nrow(overalldf),function(x,odf){
			yrv<-odf[x,"numYears"];stv<-overalldf[x,"numSites"];repv<-yrv*stv;wkv<-odf[x,"numWeeks"]
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
		},odf=overalldf)

overalldf$scenario<-unlist(scenvals)

#######################
scenvalsbs<-lapply(1:nrow(bysitedf),function(x,odf){
			yrv<-odf[x,"numYears"];stv<-overalldf[x,"numSites"];repv<-yrv*stv;wkv<-odf[x,"numWeeks"]
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
		},odf=bysitedf)

bysitedf$scenario<-unlist(scenvalsbs)
#######################
write.csv(bysitedf,file="c:/users/lsalas/desktop/bysite_additive.csv")
write.csv(overalldf,file="c:/users/lsalas/desktop/overall_additive.csv")


