
markov.4states.rsadd <- function(times1, times2, sequences, weights=NULL, dist, cuts.12=NULL,cuts.13=NULL,cuts.14=NULL,cuts.23=NULL,cuts.24=NULL,ini.dist.12=NULL, ini.dist.13=NULL, ini.dist.14=NULL, ini.dist.23=NULL, ini.dist.24=NULL, cov.12=NULL, init.cov.12=NULL, names.12=NULL, cov.13=NULL, init.cov.13=NULL, names.13=NULL, cov.14=NULL, init.cov.14=NULL, names.14=NULL, cov.23=NULL, init.cov.23=NULL, names.23=NULL, cov.24=NULL, init.cov.24=NULL, names.24=NULL, p.age, p.sex, p.year, p.rate.table, conf.int=TRUE, silent=TRUE, precision=10^(-6))
{

#check conditions
if (missing(times1)) 
        stop("Argument 'times1' is missing with no default")
if (missing(times2)) 
        stop("Argument 'times2' is missing with no default")
if (missing(sequences)) 
        stop("Argument 'sequences' is missing with no default")
if (missing(dist)) 
        stop("Argument 'dist' is missing with no default")
if (missing(p.age)) 
        stop("Argument 'p.age' is missing with no default")
if (missing(p.sex)) 
        stop("Argument 'p.sex' is missing with no default")
if (missing(p.year)) 
        stop("Argument 'p.year' is missing with no default")        
if (missing(p.rate.table)) 
        stop("Argument 'p.rate.table' is missing with no default") 
               
if (!is.vector(times1) | !is.numeric(times1))
        stop("Argument 'times1' must be a numeric vector")
if (min(times1,na.rm=T)<0)
		stop("Negative values for 'times1' are not allowed")
if (is.na(min(times1)))
		warning("individuals with missing values for 'times1' are deleted \n")
if (!is.vector(times2) | !is.numeric(times2))
        stop("Argument 'times2' must be a numeric vector")
if (min(times2,na.rm=T)<0)
		stop("Negative values for 'times2' are not allowed")
		
if (!is.vector(sequences) | !is.numeric(sequences) | (min(names(table(sequences)) %in% c(1,12,13,14,123,124))==0) )
        stop("Argument 'sequences' must be a numeric vector with values 1, 12, 13, 14, 123, or 124")
if (min( c(1,12,13,14,123,124) %in% names(table(sequences))) ==0)
        warning("all sequencess (1, 12, 13, 14, 123, 124) are not present \n ")
        
if (min(length(times1),length(times2),length(sequences)) != max(length(times1),length(times2),length(sequences)))
        stop("Arguments 'times1', 'times2', and 'sequences' need to have the same number of rows")
if (!all(is.na(times2[which(sequences==1 | sequences==13 | sequences==14)])))		
		stop("Arguments 'times2' should be NA for right-censored individuals in X=1 or individuals who directly transited from X=1 to X=3 or from X=1 to X=4")
if (min(times2-times1,na.rm=T)<=0)
		stop("Arugment 'times2' have to higher than 'times1'")
                
if(!is.null(weights))
{		
	if (!is.vector(weights) | !is.numeric(weights))
        stop("Argument 'weights' must be a numeric vector")
	if (min(weights,na.rm=T)<0)
		stop("Negative values for 'weights' are not allowed")
	if (is.na(min(weights)))
		warning("individuals with missing values for 'weights' will be removed from the analysis \n")
}

if(!(dist[1] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 12")} 
if(!(dist[2] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 13")} 
if(!(dist[3] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 14")} 
if(!(dist[4] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 23")} 
if(!(dist[5] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 24")} 
 
if(dist[1]!="PE" & (!is.null(cuts.12)))
 {stop("Arguments 'cuts.12' is only allowed for piecewise exponential distribution (PE for the first argument in 'dist')")} 
if(dist[2]!="PE" & (!is.null(cuts.13)))
 {stop("Arguments 'cuts.13' is only allowed for piecewise exponential distribution (PE for the second argument in 'dist')")}
if(dist[3]!="PE" & (!is.null(cuts.14)))
 {stop("Arguments 'cuts.14' is only allowed for piecewise exponential distribution (PE for the third argument in 'dist')")}  
if(dist[4]!="PE" & (!is.null(cuts.23)))
 {stop("Arguments 'cuts.23' is only allowed for piecewise exponential distribution (PE for the fourth argument in 'dist')")}
if(dist[5]!="PE" & (!is.null(cuts.24)))
 {stop("Arguments 'cuts.24' is only allowed for piecewise exponential distribution (PE for the fifth argument in 'dist')")}  
 
if(dist[1]=="PE" & !is.null(cuts.12))
 {
 if (!all(is.numeric(cuts.12)) | !all(!is.na(cuts.12)) | !all(cuts.12>0) | !all(is.finite(cuts.12)) | is.unsorted(cuts.12)) 
 {stop("Arguments 'cuts.12' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
if(dist[1]=="PE" & !is.null(cuts.12))
{
 if (max(cuts.12)>=max(times1,na.rm=T)) 
 {stop("Arguments 'cuts.12': check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for times1)")}
}
if(dist[2]=="PE" & !is.null(cuts.13))
 {
 if (!all(is.numeric(cuts.13)) | !all(!is.na(cuts.13)) | !all(cuts.13>0) | !all(is.finite(cuts.13)) | is.unsorted(cuts.13)) 
 {stop("Arguments 'cuts.13' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
if(dist[2]=="PE" & !is.null(cuts.13))
{
 if (max(cuts.13)>=max(times1,na.rm=T)) 
 {stop("Arguments 'cuts.13': check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for times1)")}
}
if(dist[3]=="PE" & !is.null(cuts.14))
 {
 if (!all(is.numeric(cuts.14)) | !all(!is.na(cuts.14)) | !all(cuts.14>0) | !all(is.finite(cuts.14)) | is.unsorted(cuts.14)) 
 {stop("Arguments 'cuts.14' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
if(dist[3]=="PE" & !is.null(cuts.14))
{
 if (max(cuts.14)>=max(times1,na.rm=T)) 
 {stop("Arguments 'cuts.14': check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for times1)")}
}
if(dist[4]=="PE" & !is.null(cuts.23))
 {
 if (!all(is.numeric(cuts.23)) | !all(!is.na(cuts.23)) | !all(cuts.23>0) | !all(is.finite(cuts.23)) | is.unsorted(cuts.23)) 
 {stop("Arguments 'cuts.23' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
if(dist[4]=="PE" & !is.null(cuts.23))
{
 if (max(cuts.23)>=max(times1,na.rm=T)) 
 {stop("Arguments 'cuts.23': check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for times1)")}
}
if(dist[5]=="PE" & !is.null(cuts.24))
 {
 if (!all(is.numeric(cuts.24)) | !all(!is.na(cuts.24)) | !all(cuts.24>0) | !all(is.finite(cuts.24)) | is.unsorted(cuts.24)) 
 {stop("Arguments 'cuts.24' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
if(dist[5]=="PE" & !is.null(cuts.24))
{
 if (max(cuts.24)>=max(times1,na.rm=T)) 
 {stop("Arguments 'cuts.24': check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for times1)")}
}
 
if(!is.null(ini.dist.12) & !is.numeric(ini.dist.12))
 {stop("Argument 'ini.dist.12' must be a numeric vector (default is NULL)")} 
if(!is.null(ini.dist.13) & !is.numeric(ini.dist.13))
 {stop("Argument 'ini.dist.13' must be a numeric vector (default is NULL)")} 
if(!is.null(ini.dist.14) & !is.numeric(ini.dist.13))
 {stop("Argument 'ini.dist.14' must be a numeric vector (default is NULL)")} 
if(!is.null(ini.dist.23) & !is.numeric(ini.dist.23))
 {stop("Argument 'ini.dist.23' must be a numeric vector (default is NULL)")}  
if(!is.null(ini.dist.24) & !is.numeric(ini.dist.23))
 {stop("Argument 'ini.dist.24' must be a numeric vector (default is NULL)")}  

if(dist[1]=="PE" & !is.null(ini.dist.12) & length(ini.dist.12)!=(length(cuts.12)+1))
 {stop("Incorrect number of parameters initialized for transition 12 (piecewise model)")} 
if(dist[2]=="PE" & !is.null(ini.dist.13) & length(ini.dist.13)!=(length(cuts.13)+1))
 {stop("Incorrect number of parameters initialized for transition 13 (piecewise model)")}
if(dist[3]=="PE" & !is.null(ini.dist.14) & length(ini.dist.14)!=(length(cuts.14)+1))
 {stop("Incorrect number of parameters initialized for transition 14 (piecewise model)")}
if(dist[4]=="PE" & !is.null(ini.dist.23) & length(ini.dist.23)!=(length(cuts.23)+1))
 {stop("Incorrect number of parameters initialized for transition 23 (piecewise model)")}
if(dist[5]=="PE" & !is.null(ini.dist.24) & length(ini.dist.24)!=(length(cuts.24)+1))
 {stop("Incorrect number of parameters initialized for transition 24 (piecewise model)")} 

if( (dist[1]=="E" & is.null(cuts.12) & !is.null(ini.dist.12) & length(ini.dist.12)!=1) )
 {stop("Exponential distribution (transition 12) needs initialization of one parameter")} 
if( (dist[1]=="W" & is.null(cuts.12) & !is.null(ini.dist.12) & length(ini.dist.12)!=2) )
 {stop("Weibull distribution (transition 12) needs initialization of two parameters")} 
if( (dist[1]=="WG" & is.null(cuts.12) & !is.null(ini.dist.12) & length(ini.dist.12)!=3) )
 {stop("Generalized Weibull distribution (transition 12) needs initialization of three parameters")}  

if( (dist[2]=="E" & is.null(cuts.13) & !is.null(ini.dist.13) & length(ini.dist.13)!=1) )
 {stop("Exponential distribution (transition 13) needs initialization of one parameter")} 
if( (dist[2]=="W" & is.null(cuts.13) & !is.null(ini.dist.13) & length(ini.dist.13)!=2) )
 {stop("Weibull distribution (transition 13) needs initialization of two parameters")} 
if( (dist[2]=="WG" & is.null(cuts.13) & !is.null(ini.dist.13) & length(ini.dist.13)!=3) )
 {stop("Generalized Weibull distribution (transition 13) needs initialization of three parameters")}
 
 if( (dist[3]=="E" & is.null(cuts.14) & !is.null(ini.dist.14) & length(ini.dist.14)!=1) )
 {stop("Exponential distribution (transition 14) needs initialization of one parameter")} 
if( (dist[3]=="W" & is.null(cuts.14) & !is.null(ini.dist.14) & length(ini.dist.14)!=2) )
 {stop("Weibull distribution (transition 14) needs initialization of two parameters")} 
if( (dist[3]=="WG" & is.null(cuts.14) & !is.null(ini.dist.14) & length(ini.dist.14)!=3) )
 {stop("Generalized Weibull distribution (transition 14) needs initialization of three parameters")}
 
 if( (dist[4]=="E" & is.null(cuts.23) & !is.null(ini.dist.23) & length(ini.dist.23)!=1) )
 {stop("Exponential distribution (transition 23) needs initialization of one parameter")} 
if( (dist[4]=="W" & is.null(cuts.23) & !is.null(ini.dist.23) & length(ini.dist.23)!=2) )
 {stop("Weibull distribution (transition 23) needs initialization of two parameters")} 
if( (dist[4]=="WG" & is.null(cuts.23) & !is.null(ini.dist.23) & length(ini.dist.23)!=3) )
 {stop("Generalized Weibull distribution (transition 23) needs initialization of three parameters")}

 if( (dist[5]=="E" & is.null(cuts.24) & !is.null(ini.dist.24) & length(ini.dist.24)!=1) )
 {stop("Exponential distribution (transition 24) needs initialization of one parameter")} 
if( (dist[5]=="W" & is.null(cuts.24) & !is.null(ini.dist.24) & length(ini.dist.24)!=2) )
 {stop("Weibull distribution (transition 24) needs initialization of two parameters")} 
if( (dist[5]=="WG" & is.null(cuts.24) & !is.null(ini.dist.24) & length(ini.dist.24)!=3) )
 {stop("Generalized Weibull distribution (transition 24) needs initialization of three parameters")}
 

if(!is.null(cov.12))
{
if ((!is.vector(cov.12) & !is.data.frame(cov.12) & !is.matrix(cov.12)) | !all(sapply(cov.12,is.numeric)))
 {stop("Argument 'cov.12' must be a numeric matrix or data.frame (default is NULL)")} 
if (nrow(data.frame(cov.12))!=length(times1))
 {stop("Argument 'cov.12' needs to have the same number of rows than 'times1'")}
if (sum(apply(sapply(data.frame(cov.12),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.12),is.na),1,sum))," individuals with missing values on 'cov.12' will be removed from the analysis \n")
if(!is.null(init.cov.12))
	{
	if (!is.numeric(init.cov.12))
	{stop("Argument 'init.cov.12' must be a numeric vector (default is NULL)")} 
	if (ncol(data.frame(cov.12))!=length(init.cov.12))
	{stop("Argument 'init.cov.12' needs to have the same length than number of columns of 'cov.12'")}
	}
if (!is.null(names.12))
	{
	if (!is.character(names.12))
	{stop("Argument 'names.12' must be a character vector (default is NULL)")} 
	if (ncol(data.frame(cov.12))!=length(names.12))
	{stop("Argument 'names.12' needs to have the same length than number of columns of 'cov.12'")}
	}
}

if(!is.null(cov.13))
{
if ((!is.vector(cov.13) & !is.data.frame(cov.13) & !is.matrix(cov.13)) | !all(sapply(cov.13,is.numeric)))
 {stop("Argument 'cov.13' must be a numeric matrix or data.frame (default is NULL)")} 
if (nrow(data.frame(cov.13))!=length(times1))
 {stop("Argument 'cov.13' needs to have the same number of rows than 'times1'")}
if (sum(apply(sapply(data.frame(cov.13),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.13),is.na),1,sum))," individuals with missing values on 'cov.13' will be removed from the analysis \n")
if(!is.null(init.cov.13))
	{
	if (!is.numeric(init.cov.13))
	{stop("Argument 'init.cov.13' must be a numeric vector (default is NULL)")} 
	if (ncol(data.frame(cov.13))!=length(init.cov.13))
	{stop("Argument 'init.cov.13' needs to have the same length than number of columns of 'cov.13'")}
	}
if (!is.null(names.13))
	{
	if (!is.character(names.13))
	{stop("Argument 'names.13' must be a character vector (default is NULL)")} 
	if (ncol(data.frame(cov.13))!=length(names.13))
	{stop("Argument 'names.13' needs to have the same length than number of columns of 'cov.13'")}
	}
}

if(!is.null(cov.14))
{
if ((!is.vector(cov.14) & !is.data.frame(cov.14) & !is.matrix(cov.14)) | !all(sapply(cov.14,is.numeric)))
 {stop("Argument 'cov.14' must be a numeric matrix or data.frame (default is NULL)")} 
if (nrow(data.frame(cov.14))!=length(times1))
 {stop("Argument 'cov.14' needs to have the same number of rows than 'times1'")}
if (sum(apply(sapply(data.frame(cov.14),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.14),is.na),1,sum))," individuals with missing values on 'cov.14' will be removed from the analysis \n")
if(!is.null(init.cov.14))
	{
	if (!is.numeric(init.cov.14))
	{stop("Argument 'init.cov.14' must be a numeric vector (default is NULL)")} 
	if (ncol(data.frame(cov.14))!=length(init.cov.14))
	{stop("Argument 'init.cov.14' needs to have the same length than number of columns of 'cov.14'")}
	}
if (!is.null(names.14))
	{
	if (!is.character(names.14))
	{stop("Argument 'names.14' must be a character vector (default is NULL)")} 
	if (ncol(data.frame(cov.14))!=length(names.14))
	{stop("Argument 'names.14' needs to have the same length than number of columns of 'cov.14'")}
	}
}

if(!is.null(cov.23))
{
if ((!is.vector(cov.23) & !is.data.frame(cov.23) & !is.matrix(cov.23)) | !all(sapply(cov.23,is.numeric)))
 {stop("Argument 'cov.23' must be a numeric matrix or data.frame (default is NULL)")} 
if (nrow(data.frame(cov.23))!=length(times2))
 {stop("Argument 'cov.23' needs to have the same number of rows than 'times2'")}
if (sum(apply(sapply(data.frame(cov.23),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.23),is.na),1,sum))," individuals with missing values on 'cov.23' will be removed from the analysis \n")
if(!is.null(init.cov.23))
	{
	if (!is.numeric(init.cov.23))
	{stop("Argument 'init.cov.23' must be a numeric vector (default is NULL)")} 
	if (ncol(data.frame(cov.23))!=length(init.cov.23))
	{stop("Argument 'init.cov.23' needs to have the same length than number of columns of 'cov.23'")}
	}
if (!is.null(names.23))
	{
	if (!is.character(names.23))
	{stop("Argument 'names.23' must be a character vector (default is NULL)")} 
	if (ncol(data.frame(cov.23))!=length(names.23))
	{stop("Argument 'names.23' needs to have the same length than number of columns of 'cov.23'")}
	}
} 
 
if(!is.null(cov.24))
{
if ((!is.vector(cov.24) & !is.data.frame(cov.24) & !is.matrix(cov.24)) | !all(sapply(cov.24,is.numeric)))
 {stop("Argument 'cov.24' must be a numeric matrix or data.frame (default is NULL)")} 
if (nrow(data.frame(cov.24))!=length(times2))
 {stop("Argument 'cov.24' needs to have the same number of rows than 'times1'")}
if (sum(apply(sapply(data.frame(cov.24),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.24),is.na),1,sum))," individuals with missing values on 'cov.24' will be removed from the analysis \n")
if(!is.null(init.cov.24))
	{
	if (!is.numeric(init.cov.24))
	{stop("Argument 'init.cov.24' must be a numeric vector (default is NULL)")} 
	if (ncol(data.frame(cov.24))!=length(init.cov.24))
	{stop("Argument 'init.cov.24' needs to have the same length than number of columns of 'cov.24'")}
	}
if (!is.null(names.24))
	{
	if (!is.character(names.24))
	{stop("Argument 'names.24' must be a character vector (default is NULL)")} 
	if (ncol(data.frame(cov.24))!=length(names.24))
	{stop("Argument 'names.24' needs to have the same length than number of columns of 'cov.24'")}
	}
}   
  
if(!is.numeric(p.age))
 {stop("Argument 'p.age' must be a numeric vector")} 
if(max(p.age)<365.24)
 {warning("some ages at the baseline are very low (<1 year), check that 'p.age' is given in days")} 
if (is.na(min(p.age)))
		warning("individuls with missing values for 'p.age' will be removed from the analysis \n")
		
if(!is.character(p.sex))
 {stop("Argument 'p.sex' must be a character vector")} 
if(min(names(table(p.sex)) %in% c("female","male"))==0)
 {stop("Argument 'p.sex' must be 'male' or 'female'")}
if (is.na(min(p.sex)))
		warning("individuls with missing values for 'p.sex' will be removed from the analysis \n")
		
if(!is.numeric(p.year))
 {stop("Argument 'p.year' must be a numeric vector")} 
if(as.numeric(Sys.Date()- as.Date(max(p.year), origin = "1960-01-01"))<0)
 {stop("Some dates are greater than current date: check that 'p.year' is a date format (number of days since 01.01.1960)")}
if (is.na(min(p.year)))
		warning("individuls with missing values for 'p.year' will be removed from the analysis \n")
 
if(!is.ratetable(p.rate.table))
 {stop("p.rate.table must be a ratetable object with the expected mortality rates by age, sex, and cohort year")}

if(attributes(p.rate.table)$dimid[1]!=c("age") | attributes(p.rate.table)$dimid[2]!=c("sex") | attributes(p.rate.table)$dimid[3]!=c("year"))
 {print("check that expected mortality rates in p.rate.table are given by age, sex, and cohort year")}

if(round(abs(min(attributes(p.rate.table)$cutpoints[[1]][-1]-attributes(p.rate.table)$cutpoints[[1]][-length(attributes(p.rate.table)$cutpoints[[1]])])-365.241),3)>0.001 | round(abs(max(attributes(p.rate.table)$cutpoints[[1]][-1]-attributes(p.rate.table)$cutpoints[[1]][-length(attributes(p.rate.table)$cutpoints[[1]])])-365.241),3)>0.001)
{stop("p.rate.table must have one-year age groups (365.24 or 365.241 days)")}

if(min(attributes(p.rate.table)$cutpoints[[3]][-1]-attributes(p.rate.table)$cutpoints[[3]][-length(attributes(p.rate.table)$cutpoints[[3]])])!=365 | max(attributes(p.rate.table)$cutpoints[[3]][-1]-attributes(p.rate.table)$cutpoints[[3]][-length(attributes(p.rate.table)$cutpoints[[3]])])!=366)
 {stop("p.rate.table must have one-year time intervals (365 or 366 days)")}

if(min(p.year,na.rm=T)<min(attributes(p.rate.table)$cutpoints[[3]]))
 {stop("Some dates are lower than starting year in the table of event rates: check that 'p.year' is a date format (number of days since 01.01.1960)")}

if(!(conf.int %in% c("TRUE","FALSE")))
 {stop("Argument 'conf.int' must be TRUE or FALSE (default is TRUE)")} 

if(!is.null(precision))
 {
 if(!is.numeric(precision))
 {stop("Argument 'precision' must be numeric (default is 0)")} 
 if(precision<0)
 {stop("Argument 'precision' must be greater or equal to 0 (default is 0)")}
 }
 
if(!(silent %in% c("TRUE","FALSE")))
 {stop("Argument 'silent' must be TRUE or FALSE (default is TRUE)")} 
 
warning("The table of event rates stops at:",as.character(max(attributes(p.rate.table)$cutpoints[[3]])),"\n") 
 
coef12<-NULL
sigma12<-NULL
nu12<-NULL
theta12<-NULL

coef13<-NULL
sigma13<-NULL
nu13<-NULL
theta13<-NULL

coef23<-NULL
nu23<-NULL
theta23<-NULL
sigma23<-NULL

coef14<-NULL
sigma14<-NULL
nu14<-NULL
theta14<-NULL

coef24<-NULL
nu24<-NULL
theta24<-NULL
sigma24<-NULL

#sojourn time distributions
if(dist[1]=="WG" | dist[1]=="W" | (dist[1]=="E" & is.null(cuts.12)))
 {
 H12<-function(t,z,cuts) { exp(as.matrix(z) %*% coef12) * ((((1+(t/sigma12)^nu12))^(1/theta12))-1) }
 log.h12<-function(t,z,cuts) { (as.matrix(z) %*% coef12) - log(theta12) + ((1/theta12)-1) * log1p((t/sigma12)^nu12) + log(nu12) + (nu12-1)*log(t) - nu12*log(sigma12) }
 }

if(dist[1]=="PE" & !is.null(cuts.12))
 {
cuts.12 <- sort(cuts.12)/365.24
if ((cuts.12[1] <= 0) || (cuts.12[length(cuts.12)] == Inf)) 
   stop("'cuts.12' must be positive and finite.")
cuts.12 <- c(0, cuts.12, Inf)

H12<-function(t,z,cuts) {
 H<-rep(0,length(t))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(t>=cuts[i]))*exp(as.matrix(z) %*% coef12)*((pmin(cuts[i+1],t)-cuts[i])/sigma12[i])
  }
return(H)
rm(H)
 }
log.h12<-function(t,z,cuts) {
 log.h<-rep(0,length(t))
 for (i in (1:(length(cuts)-1)))
  {
  log.h<-log.h+(1*(t>=cuts[i])*(t<cuts[i+1]))*(as.matrix(z) %*% coef12-log(sigma12[i]))
  }
 return(log.h)
 rm(log.h)
}
}


if(dist[2]=="WG" | dist[2]=="W" | (dist[2]=="E" & is.null(cuts.13)))
 {
 H13<-function(t,z,cuts) { exp(as.matrix(z) %*% coef13) * ((((1+(t/sigma13)^nu13))^(1/theta13))-1) }
 log.h13<-function(t,z,cuts) { (as.matrix(z) %*% coef13) - log(theta13) + ((1/theta13)-1) * log1p((t/sigma13)^nu13) + log(nu13) + (nu13-1)*log(t) - nu13*log(sigma13) }
 }

if(dist[2]=="PE" & !is.null(cuts.13))
 {
cuts.13 <- sort(cuts.13)/365.24
if ((cuts.13[1] <= 0) || (cuts.13[length(cuts.13)] == Inf)) 
   stop("'cuts.13' must be positive and finite.")
cuts.13 <- c(0, cuts.13, Inf)

H13<-function(t,z,cuts) {
 H<-rep(0,length(t))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(t>=cuts[i]))*exp(as.matrix(z) %*% coef13)*((pmin(cuts[i+1],t)-cuts[i])/sigma13[i])
  }
return(H)
rm(H)
 }
log.h13<-function(t,z,cuts) {
 log.h<-rep(0,length(t))
 for (i in (1:(length(cuts)-1)))
  {
  log.h<-log.h+(1*(t>=cuts[i])*(t<cuts[i+1]))*(as.matrix(z) %*% coef13-log(sigma13[i]))
  }
 return(log.h)
 rm(log.h)
}
}


if(dist[3]=="WG" | dist[3]=="W" | (dist[3]=="E" & is.null(cuts.14)))
 {
 H14<-function(t,z,cuts) { exp(as.matrix(z) %*% coef14) * ((((1+(t/sigma14)^nu14))^(1/theta14))-1) }
 log.h14<-function(t,z,cuts) { (as.matrix(z) %*% coef14) - log(theta14) + ((1/theta14)-1) * log1p((t/sigma14)^nu14) + log(nu14) + (nu14-1)*log(t) - nu14*log(sigma14) }
 }

if(dist[3]=="PE" & !is.null(cuts.14))
 {
cuts.14 <- sort(cuts.14)/365.24
if ((cuts.14[1] <= 0) || (cuts.14[length(cuts.14)] == Inf)) 
   stop("'cuts.14' must be positive and finite.")
cuts.14 <- c(0, cuts.14, Inf)

H14<-function(t,z,cuts) {
 H<-rep(0,length(t))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(t>=cuts[i]))*exp(as.matrix(z) %*% coef14)*((pmin(cuts[i+1],t)-cuts[i])/sigma14[i])
  }
return(H)
rm(H)
 }
log.h14<-function(t,z,cuts) {
 log.h<-rep(0,length(t))
 for (i in (1:(length(cuts)-1)))
  {
  log.h<-log.h+(1*(t>=cuts[i])*(t<cuts[i+1]))*(as.matrix(z) %*% coef14-log(sigma14[i]))
  }
 return(log.h)
 rm(log.h)
}
}

if(dist[4]=="WG" | dist[4]=="W" | (dist[4]=="E" & is.null(cuts.23)))
 {
 H23<-function(t,u,z,cuts) { exp(as.matrix(z) %*% coef23) * ( ((((1+(u/sigma23)^nu23))^(1/theta23))-1) - ((((1+(t/sigma23)^nu23))^(1/theta23))-1)) }
 log.h23<-function(t,z,cuts) { (as.matrix(z) %*% coef23) - log(theta23) + ((1/theta23)-1) * log1p((t/sigma23)^nu23) + log(nu23) + (nu23-1)*log(t) - nu23*log(sigma23) }
 }
 
if(dist[4]=="PE" & !is.null(cuts.23))
 {
cuts.23 <- sort(cuts.23)/365.24
if ((cuts.23[1] <= 0) || (cuts.23[length(cuts.23)] == Inf)) 
   stop("'cuts.23' must be positive and finite.")
cuts.23 <- c(0, cuts.23, Inf)

H23<-function(t,u,z,cuts) {
 H<-rep(0,length(t))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(u>=cuts[i])*(t<=cuts[i+1]))*exp(as.matrix(z) %*% coef23)*((pmin(cuts[i+1],u)-pmax(cuts[i],t))/sigma23[i])
  }
return(H)
rm(H)
 }
log.h23<-function(t,z,cuts) {
 log.h<-rep(0,length(t))
 for (i in (1:(length(cuts)-1)))
  {
  log.h<-log.h+(1*(t>=cuts[i])*(t<cuts[i+1]))*(as.matrix(z) %*% coef23-log(sigma23[i]))
  }
 return(log.h)
 rm(log.h)
}
}

if(dist[4]=="WG" | dist[4]=="W" | (dist[4]=="E" & is.null(cuts.24)))
 {
 H24<-function(t,u,z,cuts) { exp(as.matrix(z) %*% coef24) * ( ((((1+(u/sigma24)^nu24))^(1/theta24))-1) - ((((1+(t/sigma24)^nu24))^(1/theta24))-1)) }
 log.h24<-function(t,z,cuts) { (as.matrix(z) %*% coef24) - log(theta24) + ((1/theta24)-1) * log1p((t/sigma24)^nu24) + log(nu24) + (nu24-1)*log(t) - nu24*log(sigma24) }
 }
 
if(dist[4]=="PE" & !is.null(cuts.24))
 {
cuts.24 <- sort(cuts.24)/365.24
if ((cuts.24[1] <= 0) || (cuts.24[length(cuts.24)] == Inf)) 
   stop("'cuts.24' must be positive and finite.")
cuts.24 <- c(0, cuts.24, Inf)

H24<-function(t,u,z,cuts) {
 H<-rep(0,length(t))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(u>=cuts[i])*(t<=cuts[i+1]))*exp(as.matrix(z) %*% coef24)*((pmin(cuts[i+1],u)-pmax(cuts[i],t))/sigma24[i])
  }
return(H)
rm(H)
 }
log.h24<-function(t,z,cuts) {
 log.h<-rep(0,length(t))
 for (i in (1:(length(cuts)-1)))
  {
  log.h<-log.h+(1*(t>=cuts[i])*(t<cuts[i+1]))*(as.matrix(z) %*% coef24-log(sigma24[i]))
  }
 return(log.h)
 rm(log.h)
}
}

#contributions to the log-likelihood
c1<-function(times1, z12, z13, z14, cutimes12, cutimes13, cutimes14, H4P.times1)
 {return( -H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13)- H14(times1, z14, cutimes14)- H4P.times1)}

c12<-function(times1, times2, z12, z13, z14, z23, z24, cutimes12, cutimes13, cutimes14, cutimes23, cutimes24, H4P.times2)
 {return(log.h12(times1, z12, cutimes12) - H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13) - H14(times1, z14, cutimes14) - H23(times1, times2, z23, cutimes23) - H24(times1, times2, z24, cutimes24)- H4P.times2)}

c13<-function(times1, z12, z13, z14, cutimes12, cutimes13, cutimes14, H4P.times1)
 {return(log.h13(times1, z13, cutimes13) - H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13) - H14(times1, z14, cutimes14)- H4P.times1)}

c14<-function(times1, z12, z13, z14, cutimes12, cutimes13, cutimes14, H4P.times1, h4P.times1)
 {return(log(exp(log.h14(times1, z14, cutimes14))+ h4P.times1) - H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13) - H14(times1, z14, cutimes14)- H4P.times1)}

c123<-function(times1, times2, z12, z13, z14, z23, z24, cutimes12, cutimes13, cutimes14, cutimes23, cutimes24, H4P.times2)
 {return(log.h12(times1, z12, cutimes12) - H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13) - H14(times1, z14, cutimes14) + log.h23(times2, z23, cutimes23) - H23(times1, times2, z23, cutimes23) - H24(times1, times2, z24, cutimes24) - H4P.times2)}

c124<-function(times1, times2, z12, z13, z14, z23, z24, cutimes12, cutimes13, cutimes14, cutimes23, cutimes24, H4P.times2, h4P.times2)
 {return(log.h12(times1, z12, cutimes12) - H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13) - H14(times1, z14, cutimes14) + log(exp(log.h24(times2, z24, cutimes24))+ h4P.times2) - H23(times1, times2, z23, cutimes23) - H24(times1, times2, z24, cutimes24) - H4P.times2)}


d1<-times1
d2<-times2
p.age<-p.age

#non missing data
.D <- cbind(d1, cov.12, cov.13, cov.14, cov.23, cov.24, p.age, p.sex, p.year)
.na <- (pmin(apply(!is.na(.D),FUN="min",MARGIN=1))==1)

#initialization of the parameters
if (is.null(cov.12)) {cov.12.mat <- cbind(rep(0, length(d1))); n.12 <- NULL} else { cov.12.mat <- cbind(cov.12); n.12 <- paste("covariate(s) on trans. 12:  num", 1:ncol(data.frame(cov.12))); if(!is.null(names.12)) {n.12 <- names.12} }

if (is.null(cov.13)) {cov.13.mat <- cbind(rep(0, length(d1))); n.13 <- NULL} else { cov.13.mat <- cbind(cov.13); n.13 <- paste("covariate(s) on trans. 13:  num", 1:ncol(data.frame(cov.13))); if(!is.null(names.13)) {n.13 <- names.13} }

if (is.null(cov.14)) {cov.14.mat <- cbind(rep(0, length(d1))); n.14 <- NULL} else { cov.14.mat <- cbind(cov.14); n.14 <- paste("covariate(s) on trans. 14:  num", 1:ncol(data.frame(cov.14))); if(!is.null(names.14)) {n.14 <- names.14} }

if (is.null(cov.23)) {cov.23.mat <- cbind(rep(0, length(d1))); n.23 <- NULL} else { cov.23.mat <- cbind(cov.23); n.23 <- paste("covariate(s) on trans. 23:  num", 1:ncol(data.frame(cov.23))); if(!is.null(names.23)) {n.23 <- names.23} }

if (is.null(cov.24)) {cov.24.mat <- cbind(rep(0, length(d1))); n.24 <- NULL} else { cov.24.mat <- cbind(cov.24); n.24 <- paste("covariate(s) on trans. 24:  num", 1:ncol(data.frame(cov.24))); if(!is.null(names.24)) {n.24 <- names.24} }


if (is.null(ini.dist.12)) {i.12.dist<-rep(0, 1*(dist[1]=="E" & is.null(cuts.12)) + 2*(dist[1]=="W") + 3*(dist[1]=="WG") + 1*(dist[1]=="PE" & !is.null(cuts.12))*(length(cuts.12)-1))}
 else {i.12.dist<-ini.dist.12}

if (is.null(ini.dist.13)) {i.13.dist<-rep(0, 1*(dist[2]=="E" & is.null(cuts.13)) + 2*(dist[2]=="W") + 3*(dist[2]=="WG") + 1*(dist[2]=="PE" & !is.null(cuts.13))*(length(cuts.13)-1))}
 else {i.13.dist<-ini.dist.13}

 if (is.null(ini.dist.14)) {i.14.dist<-rep(0, 1*(dist[3]=="E" & is.null(cuts.14)) + 2*(dist[3]=="W") + 3*(dist[3]=="WG") + 1*(dist[3]=="PE" & !is.null(cuts.14))*(length(cuts.14)-1))}
 else {i.14.dist<-ini.dist.14}

if (is.null(ini.dist.23)) {i.23.dist<-rep(0, 1*(dist[4]=="E" & is.null(cuts.23)) + 2*(dist[4]=="W") + 3*(dist[4]=="WG") + 1*(dist[4]=="PE" & !is.null(cuts.23))*(length(cuts.23)-1))}
 else {i.23.dist<-ini.dist.23}

 if (is.null(ini.dist.24)) {i.24.dist<-rep(0, 1*(dist[5]=="E" & is.null(cuts.24)) + 2*(dist[5]=="W") + 3*(dist[5]=="WG") + 1*(dist[5]=="PE" & !is.null(cuts.24))*(length(cuts.24)-1))}
 else {i.24.dist<-ini.dist.24}


if (!is.null(init.cov.12)) {i.12<-init.cov.12}
if (is.null(init.cov.12) & is.null(cov.12)) {i.12<-NULL}
if (is.null(init.cov.12) & !is.null(cov.12)) {i.12<-rep(0, ncol(data.frame(cov.12)))}

if (!is.null(init.cov.13)) {i.13<-init.cov.13}
if (is.null(init.cov.13) & is.null(cov.13)) {i.13<-NULL}
if (is.null(init.cov.13) & !is.null(cov.13)) {i.13<-rep(0, ncol(data.frame(cov.13)))}

if (!is.null(init.cov.14)) {i.14<-init.cov.14}
if (is.null(init.cov.14) & is.null(cov.14)) {i.14<-NULL}
if (is.null(init.cov.14) & !is.null(cov.14)) {i.14<-rep(0, ncol(data.frame(cov.14)))}

if (!is.null(init.cov.23)) {i.23<-init.cov.23}
if (is.null(init.cov.23) & is.null(cov.23)) {i.23<-NULL}
if (is.null(init.cov.23) & !is.null(cov.23)) {i.23<-rep(0, ncol(data.frame(cov.23)))}

if (!is.null(init.cov.24)) {i.24<-init.cov.24}
if (is.null(init.cov.24) & is.null(cov.24)) {i.24<-NULL}
if (is.null(init.cov.24) & !is.null(cov.24)) {i.24<-rep(0, ncol(data.frame(cov.24)))}

ini <- c(i.12.dist, i.13.dist, i.14.dist, i.23.dist, i.24.dist, i.12, i.13, i.14, i.23, i.24)

{w <- rep(1, length(d1))}

#instantaneous risk of death in general population
h4P<-function(rate.table,age,sex,year,t){
age.year<-age/365.24
t.year<-t/365.24
year.year<- as.numeric( format(as.Date(year, origin = "1960-01-01"), "%Y"))
maxyear.ratetable<-max(as.numeric(attributes(rate.table)$dimnames[[3]]))
minage.year.ratetable<-round(min(as.numeric(attributes(rate.table)$dimnames[[1]])/365.24))
return(mapply(FUN=function(age,sex,year) {rate.table[age,sex,year]},trunc(trunc(age.year)+t.year)-minage.year.ratetable+1,sex,as.character(pmin(maxyear.ratetable,trunc(year.year+t.year)))))
}

log.h4P<-function(rate.table,age,sex,year,t){log(h4P(rate.table,age,sex,year,t))}

#cumulative risk between age at year 'year' and age at year ('year' + 't')
H4P<-function(rate.table,age,sex,year,t){
age.year<-age/365.24
t.year<-t/365.24
  sumcum<-rep(0,length(age))
  j<-0   
  i<-1   
  while (i<=length(age)){
  while (j<trunc(t.year[i])) {
    sumcum[i]<-sumcum[i]+365.24*h4P(rate.table,age.year[i]*365.24,sex[i],year[i],j*365.24)
    j<-j+1
    } 
   sumcum[i]<-sumcum[i]+(t.year[i]-trunc(t.year[i]))*h4P(rate.table,age[i],sex[i],year[i],j)
   i<-i+1
   j<-0
  } #end of subject loop
  return(sumcum)
}


#parameters for contributions associated to each transition
.w1 <- w[(sequences==1 & .na)]
.d1.1 <- d1[(sequences==1 & .na)]
.c1.12 <- cov.12.mat[(sequences==1 & .na),]
.c1.13 <- cov.13.mat[(sequences==1 & .na),]
.c1.14 <- cov.14.mat[(sequences==1 & .na),]
.c1.p.age <- p.age[(sequences==1 & .na)]
.c1.p.sex <- p.sex[(sequences==1 & .na)]
.c1.p.year <- p.year[(sequences==1 & .na)]
.c1.H4P.times1<-H4P(rate.table=p.rate.table,age=.c1.p.age,sex=.c1.p.sex,year=.c1.p.year,.d1.1)

.w12 <- w[(sequences==12 & .na)]
.d12.1 <- d1[(sequences==12 & .na)]
.d12.2 <- d2[(sequences==12 & .na)]
.c12.12 <- cov.12.mat[(sequences==12 & .na),]
.c12.13 <- cov.13.mat[(sequences==12 & .na),]
.c12.14 <- cov.14.mat[(sequences==12 & .na),]
.c12.23 <- cov.23.mat[(sequences==12 & .na),]
.c12.24 <- cov.24.mat[(sequences==12 & .na),]
.c12.p.age <- p.age[(sequences==12 & .na)]
.c12.p.sex <- p.sex[(sequences==12 & .na)]
.c12.p.year <- p.year[(sequences==12 & .na)]
.c12.H4P.times2<-H4P(rate.table=p.rate.table,age=.c12.p.age,sex=.c12.p.sex,year=.c12.p.year,.d12.2)

.w13 <- w[(sequences==13 & .na)]
.d13.1 <- d1[(sequences==13 & .na)]
.c13.12 <- as.matrix(cov.12.mat[(sequences==13 & .na),])
.c13.13 <- as.matrix(cov.13.mat[(sequences==13 & .na),])
.c13.14 <- as.matrix(cov.14.mat[(sequences==13 & .na),])
.c13.p.age <- p.age[(sequences==13 & .na)]
.c13.p.sex <- p.sex[(sequences==13 & .na)]
.c13.p.year <- p.year[(sequences==13 & .na)]
.c13.H4P.times1<-H4P(rate.table=p.rate.table,age=.c13.p.age,sex=.c13.p.sex,year=.c13.p.year,.d13.1)

.w14 <- w[(sequences==14 & .na)]
.d14.1 <- d1[(sequences==14 & .na)]
.c14.12 <- cov.12.mat[(sequences==14 & .na),]
.c14.13 <- cov.13.mat[(sequences==14 & .na),]
.c14.14 <- cov.14.mat[(sequences==14 & .na),]
.c14.p.age <- p.age[(sequences==14 & .na)]
.c14.p.sex <- p.sex[(sequences==14 & .na)]
.c14.p.year <- p.year[(sequences==14 & .na)]
.c14.H4P.times1<-H4P(rate.table=p.rate.table,age=.c14.p.age,sex=.c14.p.sex,year=.c14.p.year,.d14.1)
.c14.h4P.times1<-h4P(rate.table=p.rate.table,age=.c14.p.age,sex=.c14.p.sex,year=.c14.p.year,.d14.1)

.w123 <- w[(sequences==123 & .na)]
.d123.1 <- d1[(sequences==123 & .na)]
.d123.2 <- d2[(sequences==123 & .na)]
.c123.12 <- cov.12.mat[(sequences==123 & .na),]
.c123.13 <- cov.13.mat[(sequences==123 & .na),]
.c123.14 <- cov.14.mat[(sequences==123 & .na),]
.c123.23 <- cov.23.mat[(sequences==123 & .na),]
.c123.24 <- cov.24.mat[(sequences==123 & .na),]
.c123.p.age <- p.age[(sequences==123 & .na)]
.c123.p.sex <- p.sex[(sequences==123 & .na)]
.c123.p.year <- p.year[(sequences==123 & .na)]
.c123.H4P.times2<-H4P(rate.table=p.rate.table,age=.c123.p.age,sex=.c123.p.sex,year=.c123.p.year,.d123.2)

.w124 <- w[(sequences==124 & .na)]
.d124.1 <- d1[(sequences==124 & .na)]
.d124.2 <- d2[(sequences==124 & .na)]
.c124.12 <- cov.12.mat[(sequences==124 & .na),]
.c124.13 <- cov.13.mat[(sequences==124 & .na),]
.c124.14 <- cov.14.mat[(sequences==124 & .na),]
.c124.23 <- cov.23.mat[(sequences==124 & .na),]
.c124.24 <- cov.24.mat[(sequences==124 & .na),]
.c124.p.age <- p.age[(sequences==124 & .na)]
.c124.p.sex <- p.sex[(sequences==124 & .na)]
.c124.p.year <- p.year[(sequences==124 & .na)]
.c124.H4P.times2<-H4P(rate.table=p.rate.table,age=.c124.p.age,sex=.c124.p.sex,year=.c124.p.year,.d124.2)                                            
.c124.h4P.times2<-h4P(rate.table=p.rate.table,age=.c124.p.age,sex=.c124.p.sex,year=.c124.p.year,.d124.2)

#log-likelihood
logV<-function(x)
{
if (dist[1]=="E" & is.null(cuts.12)) {assign("sigma12",exp(x[1]), inherits = TRUE); assign("nu12",1, inherits = TRUE); assign("theta12",1, inherits = TRUE); i<-1}
if (dist[1]=="W") { assign("sigma12",exp(x[1]), inherits = TRUE); assign("nu12",exp(x[2]), inherits = TRUE); assign("theta12",1, inherits = TRUE); i<-2}
if (dist[1]=="WG") {assign("sigma12",exp(x[1]), inherits = TRUE); assign("nu12",exp(x[2]), inherits = TRUE); assign("theta12",exp(x[3]), inherits = TRUE); i<-3}
if (dist[1]=="PE" & !is.null(cuts.12)) {assign("sigma12",exp(x[1:(length(cuts.12)-1)]), inherits = TRUE); i<-(length(cuts.12)-1)}

if (dist[2]=="E" & is.null(cuts.13)) {assign("sigma13",exp(x[i+1]), inherits = TRUE); assign("nu13",1, inherits = TRUE); assign("theta13",1, inherits = TRUE); i<-i+1}
if (dist[2]=="W") { assign("sigma13",exp(x[i+1]), inherits = TRUE); assign("nu13",exp(x[i+2]), inherits = TRUE); assign("theta13",1, inherits = TRUE); i<-i+2}
if (dist[2]=="WG") {assign("sigma13",exp(x[i+1]), inherits = TRUE); assign("nu13",exp(x[i+2]), inherits = TRUE); assign("theta13",exp(x[i+3]), inherits = TRUE); i<-i+3}
if (dist[2]=="PE" & !is.null(cuts.13)) {assign("sigma13",exp(x[(i+1):(i+length(cuts.13)-1)]), inherits = TRUE); i<-(i+length(cuts.13)-1)}

if (dist[3]=="E" & is.null(cuts.14)) {assign("sigma14",exp(x[i+1]), inherits = TRUE); assign("nu14",1, inherits = TRUE); assign("theta14",1, inherits = TRUE); i<-i+1}
if (dist[3]=="W") { assign("sigma14",exp(x[i+1]), inherits = TRUE); assign("nu14",exp(x[i+2]), inherits = TRUE); assign("theta14",1, inherits = TRUE); i<-i+2}
if (dist[3]=="WG") {assign("sigma14",exp(x[i+1]), inherits = TRUE); assign("nu14",exp(x[i+2]), inherits = TRUE); assign("theta14",exp(x[i+3]), inherits = TRUE); i<-i+3}
if (dist[3]=="PE" & !is.null(cuts.14)) {assign("sigma14",exp(x[(i+1):(i+length(cuts.14)-1)]), inherits = TRUE); i<-(i+length(cuts.14)-1)}

if (dist[4]=="E" & is.null(cuts.23)) {assign("sigma23",exp(x[i+1]), inherits = TRUE); assign("nu23",1, inherits = TRUE); assign("theta23",1, inherits = TRUE); i<-i+1}
if (dist[4]=="W") { assign("sigma23",exp(x[i+1]), inherits = TRUE); assign("nu23",exp(x[i+2]), inherits = TRUE); assign("theta23",1, inherits = TRUE); i<-i+2}
if (dist[4]=="WG") {assign("sigma23",exp(x[i+1]), inherits = TRUE); assign("nu23",exp(x[i+2]), inherits = TRUE); assign("theta23",exp(x[i+3]), inherits = TRUE); i<-i+3}
if (dist[4]=="PE" & !is.null(cuts.23)) {assign("sigma23",exp(x[(i+1):(i+length(cuts.23)-1)]), inherits = TRUE); i<-(i+length(cuts.23)-1)}

if (dist[5]=="E" & is.null(cuts.24)) {assign("sigma24",exp(x[i+1]), inherits = TRUE); assign("nu24",1, inherits = TRUE); assign("theta24",1, inherits = TRUE); i<-i+1}
if (dist[5]=="W") { assign("sigma24",exp(x[i+1]), inherits = TRUE); assign("nu24",exp(x[i+2]), inherits = TRUE); assign("theta24",1, inherits = TRUE); i<-i+2}
if (dist[5]=="WG") {assign("sigma24",exp(x[i+1]), inherits = TRUE); assign("nu24",exp(x[i+2]), inherits = TRUE); assign("theta24",exp(x[i+3]), inherits = TRUE); i<-i+3}
if (dist[5]=="PE" & !is.null(cuts.24)) {assign("sigma24",exp(x[(i+1):(i+length(cuts.24)-1)]), inherits = TRUE); i<-(i+length(cuts.24)-1)}

if (is.null(cov.12)) {assign("coef12",0, inherits = TRUE)}
 else {assign("coef12",x[(i+1):(i+ncol(data.frame(cov.12)))], inherits = TRUE); i <-i+ncol(data.frame(cov.12))}
 if (is.null(cov.13)) {assign("coef13",0, inherits = TRUE)}
 else {assign("coef13",x[(i+1):(i+ncol(data.frame(cov.13)))], inherits = TRUE); i <-i+ncol(data.frame(cov.13))}
if (is.null(cov.14)) {assign("coef14",0, inherits = TRUE)}
 else {assign("coef14",x[(i+1):(i+ncol(data.frame(cov.14)))], inherits = TRUE); i <-i+ncol(data.frame(cov.14))}
if (is.null(cov.23)) {assign("coef23",0, inherits = TRUE)}
 else {assign("coef23",x[(i+1):(i+ncol(data.frame(cov.23)))], inherits = TRUE); i <-i+ncol(data.frame(cov.23))}
if (is.null(cov.24)) {assign("coef24",0, inherits = TRUE)}
 else {assign("coef24",x[(i+1):(i+ncol(data.frame(cov.24)))], inherits = TRUE); i <-i+ncol(data.frame(cov.24))}

return( -1*(
 sum( .w1 * c1(.d1.1, .c1.12, .c1.13, .c1.14, cuts.12, cuts.13, cuts.14, .c1.H4P.times1) ) +
 sum( .w12 * c12(.d12.1, .d12.2, .c12.12, .c12.13, .c12.14, .c12.23, .c12.24, cuts.12, cuts.13, cuts.14, cuts.23, cuts.24, .c12.H4P.times2) ) +
 sum( .w13 * c13(.d13.1, .c13.12, .c13.13, .c13.14, cuts.12, cuts.13, cuts.14, .c13.H4P.times1) ) +
 sum( .w14 * c14(.d14.1, .c14.12, .c14.13, .c14.14, cuts.12, cuts.13, cuts.14, .c14.H4P.times1, .c14.h4P.times1) ) +
 sum( .w123 * c123(.d123.1, .d123.2, .c123.12, .c123.13, .c123.14, .c123.23, .c123.24, cuts.12, cuts.13, cuts.14, cuts.23, cuts.24, .c123.H4P.times2) ) +
 sum( .w124 * c124(.d124.1, .d124.2, .c124.12, .c124.13, .c124.14, .c124.23, .c124.24, cuts.12, cuts.13, cuts.14, cuts.23, cuts.24, .c124.H4P.times2, .c124.h4P.times2) ) ) )

}

#warning("logV(ini)=",logV(ini), "\n")


#first maximum likelihood optimization
n<-1
res<-tryCatch(optim(ini, logV, hessian=conf.int, control=list(maxit=100000)))

if(inherits(res, "error"))  {
warning("Maximum likelihood optimization fails to converge", "\n")
  } else  {
if(silent==FALSE) {warning(-1*res$value, "\n")}  

#further maximum likelihood optimizations
if(is.null(precision)) {delta <- 10^(-6)} else {delta <-precision}

while (n<=2 & !(inherits(res, "error"))) {
temp.value<-res$value
res<-tryCatch(optim(res$par, logV, hessian=conf.int, control=list(maxit=100000)))

if (!(inherits(res, "error"))) {
   n<-1*((temp.value-res$value)>delta) + (n+1)*((temp.value-res$value)<=delta)
   if(silent==FALSE) {warning(-1*res$value, "\n")} }
   }
if(inherits(res, "error")) {
warning("Maximum likelihood optimization fails to converge", "\n")
  } else {

#output
if (conf.int==TRUE) {
  if (max(!is.na(tryCatch(solve(res$hessian), error=function(e) NA)),na.rm=F)==1){
  table.res <- data.frame(Estimate = round(res$par, 4),
  SE = round(sqrt(diag(solve(res$hessian))), 4),
  Wald = round(res$par/sqrt(diag(solve(res$hessian))), 4),
  Pvalue = round(2*(1-pnorm(abs(res$par/sqrt(diag(solve(res$hessian)))), 0, 1)) , 4) )
  names(table.res)<-c("Estimate","Std.Error","t.value","Pr(>|t|)")
  table.covariance<-res$hessian
  }
  else {
  table.res <- data.frame(Estimate = round(res$par, 4) )
  table.covariance<-NULL
  warning("Hessian matrix not defined", "\n")
  } #end else for hessian matrix condition
}

if (conf.int==FALSE) {
table.res <- data.frame(Estimate = round(res$par, 4) )
table.covariance<-NULL
}

if (dist[1]=="E" & is.null(cuts.12))  { lab12<-c("log(sigma) on trans. 12")}
if (dist[1]=="W" & is.null(cuts.12))  { lab12<-c("log(sigma) on trans. 12", "log(nu) on trans. 12")}
if (dist[1]=="WG" & is.null(cuts.12)) { lab12<-c("log(sigma) on trans. 12", "log(nu) on trans. 12", "log(theta) on trans. 12")}
if (dist[1]=="PE" & !is.null(cuts.12)) {
 lab12<-rep("",length(cuts.12)-1)
 for (i in (1:(length(cuts.12)-1)))
  {
  lab12[i]<-paste("log(sigma) on trans. 12, interval [",round(cuts.12[i],3),";",round(cuts.12[i+1],3),"[",sep="")
  }
 }

if (dist[2]=="E" & is.null(cuts.13))  { lab13<-c("log(sigma) on trans. 13")}
if (dist[2]=="W" & is.null(cuts.13))  { lab13<-c("log(sigma) on trans. 13", "log(nu) on trans. 13")}
if (dist[2]=="WG" & is.null(cuts.13)) { lab13<-c("log(sigma) on trans. 13", "log(nu) on trans. 13", "log(theta) on trans. 13")}
if (dist[2]=="PE" & !is.null(cuts.13)) {
 lab13<-rep("",length(cuts.13)-1)
 for (i in (1:(length(cuts.13)-1)))
  {
  lab13[i]<-paste("log(sigma) on trans. 13, interval [",round(cuts.13[i],3),";",round(cuts.13[i+1],3),"[",sep="")
  }
 }

if (dist[3]=="E" & is.null(cuts.14))  { lab14<-c("log(sigma) on trans. 1E")}
if (dist[3]=="W" & is.null(cuts.14))  { lab14<-c("log(sigma) on trans. 1E", "log(nu) on trans. 1E")}
if (dist[3]=="WG" & is.null(cuts.14)) { lab14<-c("log(sigma) on trans. 1E", "log(nu) on trans. 1E", "log(theta) on trans. 1E")}
if (dist[3]=="PE" & !is.null(cuts.14)) {
 lab14<-rep("",length(cuts.14)-1)
 for (i in (1:(length(cuts.14)-1)))
  {
  lab14[i]<-paste("log(sigma) on trans. 1E, interval [",round(cuts.14[i],3),";",round(cuts.14[i+1],3),"[",sep="")
  }
 }

if (dist[4]=="E" & is.null(cuts.23))  { lab23<-c("log(sigma) on trans. 23")}
if (dist[4]=="W" & is.null(cuts.23))  { lab23<-c("log(sigma) on trans. 23", "log(nu) on trans. 23")}
if (dist[4]=="WG" & is.null(cuts.23)) { lab23<-c("log(sigma) on trans. 23", "log(nu) on trans. 23", "log(theta) on trans. 23")}
if (dist[4]=="PE" & !is.null(cuts.23)) {
 lab23<-rep("",length(cuts.23)-1)
 for (i in (1:(length(cuts.23)-1)))
  {
  lab23[i]<-paste("log(sigma) on trans. 23, interval [",round(cuts.23[i],3),";",round(cuts.23[i+1],3),"[",sep="")
  }
 }

if (dist[5]=="E" & is.null(cuts.24))  { lab24<-c("log(sigma) on trans. 2E")}
if (dist[5]=="W" & is.null(cuts.24))  { lab24<-c("log(sigma) on trans. 2E", "log(nu) on trans. 2E")}
if (dist[5]=="WG" & is.null(cuts.24)) { lab24<-c("log(sigma) on trans. 2E", "log(nu) on trans. 2E", "log(theta) on trans. 2E")}
if (dist[5]=="PE" & !is.null(cuts.24)) {
 lab24<-rep("",length(cuts.24)-1)
 for (i in (1:(length(cuts.24)-1)))
  {
  lab24[i]<-paste("log(sigma) on trans. 2E, interval [",round(cuts.24[i],3),";",round(cuts.24[i+1],3),"[",sep="")
  }
 }

lab<-c(lab12, lab13, lab14, lab23,lab24, n.12, n.13, n.14, n.23, n.24)

rownames(table.res) <- paste(1:length(lab), lab)

warning("\n Number of data rows:",nrow(.D))
warning("Number of data rows with missing values (deleted):",nrow(.D)-sum(.na),"\n")

return(list(
object="m4 (4-state relative survival markov model with additive risks)",
dist=dist,
cuts.12=cuts.12,
cuts.13=cuts.13,
cuts.14=cuts.14,
cuts.23=cuts.23,
cuts.24=cuts.23,
covariates=c( max(0, length(n.12)), max(0, length(n.13)), max(0, length(n.14)), max(0, length(n.23)), max(0, length(n.24)) ),
table=table.res,
cov.matrix=table.covariance,
LogLik=(-1*res$value),
AIC=2*length(res$par)-2*(-1*res$value)))
  }
  }
}


