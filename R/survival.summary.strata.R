
survival.summary.strata<-function(study, time, n.risk, surv.rate, confidence, strata)
{

#.data <- data.frame(study, time, n.risk, surv.rate)
.data <- data.frame(study, time, n.risk, surv.rate, strata)
.data <- .data[order(study, time),]
study <- .data$study
time <- .data$time
n.risk <- .data$n.risk
surv.rate <- .data$surv.rate
strata <- factor(.data$strata)

#### Determination of the time horizon for which each stratum contains at least two studies

# Number of studies
IndiceS<-sort(unique(study))
NbS<-length(IndiceS)
IndiceT<-sort(unique(time))
NbT<-length(IndiceT)
TimeM<-max(IndiceT)


verif <- function(x)
{
ref.study <- sort(time[study==x])
ref.data  <- IndiceT[IndiceT<=max(ref.study)]
return(1*(sum(ref.study!=ref.data)==0))
}

verif.data <- data.frame(study=IndiceS, check=sapply(IndiceS, FUN="verif"))

if (sum(verif.data$check==0)>0) { 
return(list(
verif.data=verif.data,
summary.fixed=NA,
summary.random=NA,
median.fixed=NA,
mean.fixed=NA,
heterogeneity=NA,
p.value=NA))
}

else {

	#strata<-factor(strata)

	Temp<-table(time, strata)
	TimeBound<-max(sort(unique(time))[apply(Temp,1,min)>=2])
	Indic<-((time <= TimeBound) & !is.na(strata))

	study<-study[Indic]
	n.risk<-n.risk[Indic]
	surv.rate<-surv.rate[Indic]
	strata<-strata[Indic]
	time<-time[Indic]

#### Calculus of the conditional survivals

	CondSurv<-rep(-99,length(surv.rate))
	for (i in 1:length(surv.rate))
	{
		if (time[i]==min(time)) {CondSurv[i]<-surv.rate[i]}
		if (time[i]!=min(time)) {CondSurv[i]<-surv.rate[i]/surv.rate[i-1]}
	}

#### Number of strata

NbStrata<-length(unique(strata))


SurvivalFE<-list(NULL)
SurvivalRE<-list(NULL)
MedianTime<-data.frame(matrix(,ncol=7,nrow=(NbStrata+1)))
names(MedianTime)<-c("Strata","MedianFE","Lower95%FE","Upper95%FE","MedianRE","Lower95%RE","Upper95%RE")
MeanTime<-data.frame(matrix(,ncol=7,nrow=(NbStrata+1)))
names(MeanTime)<-c("Strata","MeanFE","Lower95%FE","Upper95%FE","MeanRE","Lower95%RE","Upper95%RE")
Heterogeneity<-data.frame(matrix(,ncol=4,nrow=(NbStrata+1)))
names(Heterogeneity)<-c("Strata","Qstatistic","H index","I-squared%")


####  Meta-analysis in each stratum and in all studies

for (g in 1:(NbStrata+1))
{

	if (g<=NbStrata)  # only studies in the stratum g
	{
		IndicStrata<-strata==levels(strata)[g]
		StudyTemp<-study[IndicStrata]
	}
	if (g==(NbStrata+1))  ## all studies
	{
		IndicStrata<-rep(TRUE,length(strata))
		StudyTemp<-study[IndicStrata]
	}

#### Mise en forme des data

	# Arc sine tranformation of the conditional survival (Equations page 7)
	ArcSineCondSurv<-asin(sqrt((n.risk[IndicStrata]*CondSurv[IndicStrata]+0.25)/(n.risk[IndicStrata]+0.5)))
	VarArcSineCondSurv<-1/(4*(n.risk[IndicStrata]+0.5))

	# Number of studies
	IndiceStudies<-sort(unique(study[IndicStrata]))
	NbStudies<-length(IndiceStudies)

	# Number of points in time (number of pooled conditional survival to be estimated)
	IndiceTimes<-sort(unique(time[IndicStrata]))
	NbTimes<-length(IndiceTimes)
	TimeMax<-max(IndiceTimes)

	
	ys<-matrix(nrow=NbStudies,ncol=NbTimes)
	for (i in 1:NbStudies) ys[i,1:sum(StudyTemp==IndiceStudies[i])]<-ArcSineCondSurv[StudyTemp==IndiceStudies[i]]
	for (i in 1:NbTimes)  ys[is.na(ys[,i]),i]<-mean(ys[!is.na(ys[,i]),i])

	vars<-matrix(nrow=NbStudies,ncol=NbTimes)
	for (i in 1:NbStudies) vars[i,1:sum(StudyTemp==IndiceStudies[i])]<-VarArcSineCondSurv[StudyTemp==IndiceStudies[i]]
	for (i in 1:NbTimes)  vars[is.na(vars[,i]),i]<-10^12



#### Application of the multivariate DerSimonian and Laird methodology (R code from Dan Jackson)
#### Equations page 7
#### Estimation of the pooled arc sine transformed conditional survival in each stratum


	mat<-matrix(0, nrow=NbTimes, ncol=NbTimes)
	maty<-matrix(0, nrow=NbTimes, ncol=1)
	matfixedeffect<-matrix(0, nrow=NbTimes, ncol=NbTimes)
	matyfixedeffect<-matrix(0, nrow=NbTimes, ncol=1)
	esttrun<-matrix(0, nrow=NbTimes, ncol=NbTimes)
	estmat<-matrix(0, nrow=NbTimes, ncol=NbTimes)
	count<-1
	
	for(i in 1:(NbTimes-1))
	{
		for(j in (i+1):NbTimes) 
		{
			crossses<-(vars[,i]*vars[,j])^0.5
			mean1a<-sum(ys[,i]/crossses)/sum(1/crossses)
			mean2a<-sum(ys[,j]/crossses)/sum(1/crossses)
			Q<-sum((ys[,i]-mean1a)*(ys[,j]-mean2a)/crossses)
			bb<-sum(1/crossses)-sum(1/(crossses^2))/sum(1/crossses)
			est<-Q/bb
			estmat[i,j]<-est
		}
	}
	estmat<-estmat+t(estmat)
	for(i in 1:NbTimes)
	{
		cross_ses<-vars[,i]
		mean1a<-sum(ys[,i]/vars[,i])/sum(1/vars[,i])
		Q<-sum((ys[,i]-mean1a)*(ys[,i]-mean1a)/vars[,i])
		rhos<-rep(1,NbStudies)
		check<-(1-1*(vars[,i]==10^12))
		rhos<-rhos*check
		aa<-sum(rhos)-sum(rhos/vars[,i])/sum(1/vars[,i])
		bb<-sum(1/vars[,i])-sum(1/(vars[,i]^2))/sum(1/vars[,i])
		est<-(Q-aa)/bb
		estmat[i,i]<-est
	}
	eig<-eigen(estmat)
	for(i in 1:NbTimes) esttrun<-esttrun+max(0, eig$values[i])*eig$vectors[,i]%*%t(eig$vectors[,i])
	for(k in 1:NbStudies)
	{
		within<-diag(vars[k,])
		invwithin<-diag(1/vars[k,])
		matfixedeffect<-matfixedeffect+invwithin 
		mat1<-within+esttrun
		mat1<-qr.solve(mat1,diag(rep(1,length(mat1[,1]))))
		mat<-mat+mat1
		y<-matrix(ys[k,],nrow=NbTimes,ncol=1)
		mat2<-mat1%*%y
		maty<-maty+mat2
		matyfixedeffect<-matyfixedeffect+invwithin %*%y
	}
	covrandomeffect<-solve(mat)					## Variance matrix of the pooled arcsine transformed cond. surv. with random effects
	betahatrandomeffect<-covrandomeffect%*%maty     	## Pooled estimates of arcsine transformed cond. surv. with random effects
	covfixedeffect<-solve(matfixedeffect)			## Variance matrix of the pooled arcsine transformed cond. surv. with fixed effects
	betahatfixedeffect<-covfixedeffect%*%matyfixedeffect  ## Pooled estimates of arcsine transformed cond. surv. with fixed effects


#### Heterogeneity (section 2.3. of the paper)

	HeterogeneityQ<-0
	for (i in 1:NbStudies) 
	{
		IndicHeterogeneity<-vars[i,]!=10^12
		HeterogeneityQ<-HeterogeneityQ+t(ys[i,IndicHeterogeneity]-betahatfixedeffect[IndicHeterogeneity])%*%(diag(1/vars[i,]))[IndicHeterogeneity,IndicHeterogeneity]%*%(ys[i,IndicHeterogeneity]-betahatfixedeffect[IndicHeterogeneity])
	}
	HSquared<-HeterogeneityQ/(sum(as.vector(vars)!=10^12)-NbStudies)
	ISquared<-max(0,100*(HSquared-1)/HSquared)  ## in percentage
	HeterogeneityResults<-c(HeterogeneityQ,HSquared,ISquared)



#### Summarized survival probabilities

	PooledSurvivalFE<-cumprod((sin(betahatfixedeffect))^2)   	## Fixed effects
	PooledSurvivalRE<-cumprod((sin(betahatrandomeffect))^2)   	## Random effects


#### 95% confidence intervals around the summarized survival probabilities


	if (confidence=="Greenwood")
	{
		# Fixed effects
		VarLogPooledSurvivalFE<-cumsum(4*(cos(betahatfixedeffect)/sin(betahatfixedeffect))^2*diag(covfixedeffect))
		PooledSurvivalICinfFE<-exp(log(PooledSurvivalFE)-1.96*sqrt(VarLogPooledSurvivalFE))
		PooledSurvivalICsupFE<-exp(log(PooledSurvivalFE)+1.96*sqrt(VarLogPooledSurvivalFE))
		CIPooledSurvivalFE<-cbind(PooledSurvivalICinfFE,PooledSurvivalICsupFE)

		# Random effects
		DerivativeOfTransformation<-cos(betahatrandomeffect)/sin(betahatrandomeffect)
		VarLogPooledConditionalSurvivalProbaRE<-(DerivativeOfTransformation%*%t(DerivativeOfTransformation))*covrandomeffect
		MatriceTriangular<-diag(rep(1,length(betahatrandomeffect)))
		MatriceTriangular[upper.tri(MatriceTriangular)]<-1
		VarLogPooledSurvivalRE<-4*(t(MatriceTriangular)%*%VarLogPooledConditionalSurvivalProbaRE%*%MatriceTriangular)
		PooledSurvivalICinfRE<-exp(log(PooledSurvivalRE)-1.96*sqrt(diag(VarLogPooledSurvivalRE)))
		PooledSurvivalICsupRE<-exp(log(PooledSurvivalRE)+1.96*sqrt(diag(VarLogPooledSurvivalRE)))
		CIPooledSurvivalRE<-cbind(PooledSurvivalICinfRE,PooledSurvivalICsupRE)
}


	SimulatedArcSinCondSurvFE<-rmvnorm(n=10000,mean=betahatfixedeffect,sigma = covfixedeffect)
	SimulatedArcSinCondSurvRE<-rmvnorm(n=10000,mean=betahatrandomeffect,sigma = covrandomeffect)
	SimulatedPooledSurvFE<-apply((sin(SimulatedArcSinCondSurvFE))^2,1,cumprod)
	SimulatedPooledSurvRE<-apply((sin(SimulatedArcSinCondSurvRE))^2,1,cumprod)

	if (confidence=="MonteCarlo")
	{
		# Fixed effects

		CIPooledSurvivalFE<-t(apply(SimulatedPooledSurvFE,1,quantile,probs=c(0.025,0.975)))

		# Random effects

		CIPooledSurvivalRE<-t(apply(SimulatedPooledSurvRE,1,quantile,probs=c(0.025,0.975)))
	}




#### Median survival time from the summarized survival curve


	# Fixed effects
	if (PooledSurvivalFE[length(PooledSurvivalFE)]>0.50) MedianTimeFE<-NA
	if (PooledSurvivalFE[length(PooledSurvivalFE)]==0.50) MedianTimeFE<-IndiceTimes[length(IndiceTimes)]
	if (PooledSurvivalFE[length(PooledSurvivalFE)]<0.50)
	{
		IndexMin<-max((1:length(PooledSurvivalFE))[PooledSurvivalFE>0.5])
		IndexMax<-min((1:length(PooledSurvivalFE))[PooledSurvivalFE<0.5])
		MedianTimeFE<-IndiceTimes[IndexMax]-(IndiceTimes[IndexMax]-IndiceTimes[IndexMin])*(PooledSurvivalFE[IndexMax]-0.5)/(PooledSurvivalFE[IndexMax]-PooledSurvivalFE[IndexMin])
	}
	MedianTimeFEbtp<-rep(NA,10000)
	for (i in 1:10000)
	{
		if (SimulatedPooledSurvFE[length(PooledSurvivalFE),i]<0.50)
		{
			IndexMin<-max((1:length(SimulatedPooledSurvFE[,i]))[SimulatedPooledSurvFE[,i]>0.5])
			IndexMax<-min((1:length(SimulatedPooledSurvFE[,i]))[SimulatedPooledSurvFE[,i]<0.5])
			MedianTimeFEbtp[i]<-IndiceTimes[IndexMax]-(IndiceTimes[IndexMax]-IndiceTimes[IndexMin])*(SimulatedPooledSurvFE[IndexMax,i]-0.5)/(SimulatedPooledSurvFE[IndexMax,i]-SimulatedPooledSurvFE[IndexMin,i])
		}
	}
	if (sum(is.na(MedianTimeFEbtp))==0) CIMedianTimeFE<-quantile(MedianTimeFEbtp,probs=c(0.025,0.975))
	if (sum(is.na(MedianTimeFEbtp))>0) CIMedianTimeFE<-c(NA,NA)



	# Random effects
	if (PooledSurvivalRE[length(PooledSurvivalRE)]>0.50) MedianTimeRE<-NA
	if (PooledSurvivalRE[length(PooledSurvivalRE)]==0.50) MedianTimeRE<-IndiceTimes[length(IndiceTimes)]
	if (PooledSurvivalRE[length(PooledSurvivalRE)]<0.50)
	{
		IndexMin<-max((1:length(PooledSurvivalRE))[PooledSurvivalRE>0.5])
		IndexMax<-min((1:length(PooledSurvivalRE))[PooledSurvivalRE<0.5])
		MedianTimeRE<-IndiceTimes[IndexMax]-(IndiceTimes[IndexMax]-IndiceTimes[IndexMin])*(PooledSurvivalRE[IndexMax]-0.5)/(PooledSurvivalRE[IndexMax]-PooledSurvivalRE[IndexMin])
	}
	MedianTimeREbtp<-rep(NA,10000)
	for (i in 1:10000)
	{
		if (SimulatedPooledSurvRE[length(PooledSurvivalRE),i]<0.50)
		{
			IndexMin<-max((1:length(SimulatedPooledSurvRE[,i]))[SimulatedPooledSurvRE[,i]>0.5])
			IndexMax<-min((1:length(SimulatedPooledSurvRE[,i]))[SimulatedPooledSurvRE[,i]<0.5])
			MedianTimeREbtp[i]<-IndiceTimes[IndexMax]-(IndiceTimes[IndexMax]-IndiceTimes[IndexMin])*(SimulatedPooledSurvRE[IndexMax,i]-0.5)/(SimulatedPooledSurvRE[IndexMax,i]-SimulatedPooledSurvRE[IndexMin,i])
		}
	}
	if (sum(is.na(MedianTimeREbtp))==0) CIMedianTimeRE<-quantile(MedianTimeREbtp,probs=c(0.025,0.975))
	if (sum(is.na(MedianTimeREbtp))>0) CIMedianTimeRE<-c(NA,NA)




#### Mean survival time from the summarized survival curve

	DureeIntervalle<-IndiceTimes[-1]-IndiceTimes[-length(IndiceTimes)]

	#MeanTimeFE<-sum(DureeIntervalle*(PooledSurvivalFE[-length(PooledSurvivalFE)]+PooledSurvivalFE[-1])/2)
	
	MeanTimeFE <- IndiceTimes[1]*(1+PooledSurvivalFE[1])/2+sum(DureeIntervalle * (PooledSurvivalFE[-length(PooledSurvivalFE)] +PooledSurvivalFE[-1])/2)
	
	Temp<-(SimulatedPooledSurvFE[-length(PooledSurvivalFE),]+SimulatedPooledSurvFE[-1,])/2
	for (i in 1:(length(PooledSurvivalFE)-1)) Temp[i,]<-Temp[i,]*DureeIntervalle[i]
	#CIMeanTimeFE<-quantile(apply(Temp,2,sum),probs=c(0.025,0.975))
    CIMeanTimeFE <- quantile(apply(Temp, 2, sum)+IndiceTimes[1]*(SimulatedPooledSurvFE[1,]+rep(1,10000))/2, probs = c(0.025,0.975))
	
	#MeanTimeRE<-sum(DureeIntervalle*(PooledSurvivalRE[-length(PooledSurvivalRE)]+PooledSurvivalRE[-1])/2)
	
	MeanTimeRE <- IndiceTimes[1]*(1+PooledSurvivalRE[1])/2+sum(DureeIntervalle * (PooledSurvivalRE[-length(PooledSurvivalRE)] +PooledSurvivalRE[-1])/2)
	
	MeanTimeRE <- IndiceTimes[1]*(1+PooledSurvivalRE[1])/2+sum(DureeIntervalle * (PooledSurvivalRE[-length(PooledSurvivalRE)] +PooledSurvivalRE[-1])/2)
	
	Temp<-(SimulatedPooledSurvRE[-length(PooledSurvivalRE),]+SimulatedPooledSurvRE[-1,])/2
	for (i in 1:(length(PooledSurvivalRE)-1)) Temp[i,]<-Temp[i,]*DureeIntervalle[i]
	# CIMeanTimeRE<-quantile(apply(Temp,2,sum),probs=c(0.025,0.975))
    CIMeanTimeRE <- quantile(apply(Temp, 2, sum)+IndiceTimes[1]*(SimulatedPooledSurvRE[1,]+rep(1,10000))/2, probs = c(0.025,0.975))

#### Results per stratum

	SummarySurvivalFE<-data.frame(cbind(IndiceTimes,round(PooledSurvivalFE,3),round(CIPooledSurvivalFE,3)))
	names(SummarySurvivalFE)<-c("time","SurvivalFE","Lower CI95%","Upper CI95%")

	MedianSurvivalTimeFE<-round(c(MedianTimeFE,CIMedianTimeFE),2)
	MeanSurvivalTimeFE<-round(c(MeanTimeFE,CIMeanTimeFE),2)

	SummarySurvivalRE<-data.frame(cbind(IndiceTimes,round(PooledSurvivalRE,3),round(CIPooledSurvivalRE,3)))
	names(SummarySurvivalRE)<-c("time","SurvivalRE","Lower CI95%","Upper CI95%")

	MedianSurvivalTimeRE<-round(c(MedianTimeRE,CIMedianTimeRE),2)
	MeanSurvivalTimeRE<-round(c(MeanTimeRE,CIMeanTimeRE),2)


#### Put results of each stratum in the global results

	SurvivalFE[[g]]<-SummarySurvivalFE
	SurvivalRE[[g]]<-SummarySurvivalRE
	if (g<=NbStrata)
	{
		MedianTime[g,]<-c(levels(strata)[g],MedianSurvivalTimeFE,MedianSurvivalTimeRE)
		MeanTime[g,]<-c(levels(strata)[g],MeanSurvivalTimeFE,MeanSurvivalTimeRE)
		Heterogeneity[g,]<-c(levels(strata)[g],round(HeterogeneityResults,2))
	}
	if (g==(NbStrata+1))
	{
		MedianTime[g,]<-c("All strata",MedianSurvivalTimeFE,MedianSurvivalTimeRE)
		MeanTime[g,]<-c("All strata",MeanSurvivalTimeFE,MeanSurvivalTimeRE)
		Heterogeneity[g,]<-c("All strata",round(HeterogeneityResults,2))
	}
}

names(SurvivalFE)<-c(levels(strata),"All strata")
names(SurvivalRE)<-c(levels(strata),"All strata")


BetweenStudyQStatistic<-as.numeric(Heterogeneity$Qstatistic[NbStrata+1])-sum(as.numeric(Heterogeneity$Qstatistic[1:NbStrata]))
DegreeOfFreedom<-NbTimes*(NbStrata-1)
Pvalue<-1-pchisq(BetweenStudyQStatistic,df=DegreeOfFreedom)

return(list(
verif.data=verif.data,
summary.fixed=SurvivalFE,
summary.random=SurvivalRE,
median.fixed=MedianTime,
mean.fixed=MeanTime,
heterogeneity=Heterogeneity,
p.value=Pvalue))

}
}
