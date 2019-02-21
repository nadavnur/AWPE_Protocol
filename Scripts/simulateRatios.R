# TODO: Add comment
# 
# Author: lsalas
###############################################################################

## !!!!!!!!!!! THIS CODE IS INTEGRATED INTO runSurveySimulations.R  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

# This simulation starts with the simulation of the Adult abundances
# We take the true max adult abundance for an area and year and, with the ratio and the curvature for F, we have the function defining F

#This shall generate the simulated value(s) of surveyed F for each population and year

# We then set 3 dates for surveys for each year, and then SHIFT these by some random quantity (sampled from a distribution) 
# to simulate the error due to missing the peak fledgling production
# With these dates shifted, we obtain the true value for F on these shifted dates, and sample with some sampling error (from a distribution)

# Then...
# If there was only 1 sample for F, then that becomes the estimated F
# If there were 2 samples, then take the mean or max (report both)
# If there were 3 or more, then fit a parabola and get the max of the parabola

# With the estimated F value for that year and the estimated N value for that year (from the simulation of adults), calculate the simulated ratio value for that year
# F is the true number of fledglings
# season is the span in days of the breeding season
# Returns: the curvature of the parabola
getCurvature<-function(F,season){
	curv<- -1*(F/((season^2)/4))
	return(curv)
}


# F is the true number of fledglings
# curv is the curvature of the parabola, estimated through the getCurvature function, always negative
# seasonStart is the start of the survey season in days
# seasonEnd is the end of the survey season in days
# numsamples is the number of survey samples to take
# sampInt is the number of days between samples
# sdFerr is the standard deviation of the survey error
# Returns: a vector of the true value of F at the sampling dates, and the selected survey dates
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

