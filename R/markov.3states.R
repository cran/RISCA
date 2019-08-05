
markov.3states <- function(times1, times2, sequences, weights=NULL, dist, cuts.12=NULL,cuts.13=NULL,cuts.23=NULL,ini.dist.12=NULL, ini.dist.13=NULL, ini.dist.23=NULL, cov.12=NULL, init.cov.12=NULL, names.12=NULL, cov.13=NULL, init.cov.13=NULL, names.13=NULL, cov.23=NULL, init.cov.23=NULL, names.23=NULL, conf.int=TRUE, silent=TRUE, precision=10^(-6))

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
        
if (!is.vector(times1) | !is.numeric(times1))
        stop("Argument 'times1' must be a numeric vector")
if (min(times1,na.rm=T)<0)
		stop("Negative values for 'times1' are not allowed")
if (is.na(min(times1)))
		warning("individuals with missing values for 'times1' will be removed from the analysis \n")
if (!is.vector(times2) | !is.numeric(times2))
        stop("Argument 'times2' must be a numeric vector")
if (min(times2,na.rm=T)<0)
		stop("Negative values for 'times2' are not allowed")
		
if (!is.vector(sequences) | !is.numeric(sequences) | (min(names(table(sequences)) %in% c(1,12,13,123))==0) )
        stop("Argument 'sequences' must be a numeric vector with values 1, 12, 13, or 123")
if (min( c(1,12,13,123) %in% names(table(sequences)))==0)
        warning("all possible sequencess (1, 12, 13, 123) are not present \n")
        
if (min(length(times1),length(times2),length(sequences)) != max(length(times1),length(times2),length(sequences)))
        stop("Arguments 'times1', 'times2', and 'sequences' need to have the same number of rows")
if (!all(is.na(times2[which(sequences==1 | sequences==13)])))		
	stop("Argument 'times2' have to be NA for right-censoring in X=1 or individuals who directly transited from X=1 to X=3")	
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

if(length(dist)!=3)
 {stop("Argument 'dist' have to contain 3 values")} 
 
if(!(dist[1] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 12")} 
if(!(dist[2] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 13")} 
if(!(dist[3] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 23")} 
  
if(dist[1]!="PE" & (!is.null(cuts.12)))
 {stop("Arguments 'cuts.12' is only allowed for piecewise exponential distribution (PE for the first argument in 'dist')")} 
if(dist[2]!="PE" & (!is.null(cuts.13)))
 {stop("Arguments 'cuts.13' is only allowed for piecewise exponential distribution (PE for the second argument in 'dist')")} 
if(dist[3]!="PE" & (!is.null(cuts.23)))
 {stop("Arguments 'cuts.23' is only allowed for piecewise exponential distribution (PE for the third argument in 'dist')")} 
 
if(dist[1]=="PE" & !is.null(cuts.12))
 {
 if (!all(is.numeric(cuts.12)) | !all(!is.na(cuts.12)) | !all(cuts.12>0) | !all(is.finite(cuts.12)) | is.unsorted(cuts.12)) 
 {stop("Arguments 'cuts.12' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
if(dist[1]=="PE" & !is.null(cuts.12))
{
 if (max(cuts.12)>=max(times1,na.rm=T)) 
 {stop("Arguments 'cuts.12' check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for times1)")}
}
if(dist[2]=="PE" & !is.null(cuts.13))
 {
 if (!all(is.numeric(cuts.13)) | !all(!is.na(cuts.13)) | !all(cuts.13>0) | !all(is.finite(cuts.13)) | is.unsorted(cuts.13)) 
 {stop("Arguments 'cuts.13' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
if(dist[2]=="PE" & !is.null(cuts.13))
{
 if (max(cuts.13)>=max(times1,na.rm=T)) 
 {stop("Arguments 'cuts.13' check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for times1)")}
}
if(dist[3]=="PE" & !is.null(cuts.23))
 {
 if (!all(is.numeric(cuts.23)) | !all(!is.na(cuts.23)) | !all(cuts.23>0) | !all(is.finite(cuts.23)) | is.unsorted(cuts.23)) 
 {stop("Arguments 'cuts.23' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
if(dist[3]=="PE" & !is.null(cuts.23))
{
 if (max(cuts.23)>=max(times1,na.rm=T)) 
 {stop("Arguments 'cuts.23' check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for times1)")}
}
 
if(!is.null(ini.dist.12) & !is.numeric(ini.dist.12))
 {stop("Argument 'ini.dist.12' must be a numeric vector (default is NULL)")} 
if(!is.null(ini.dist.13) & !is.numeric(ini.dist.13))
 {stop("Argument 'ini.dist.13' must be a numeric vector (default is NULL)")} 
if(!is.null(ini.dist.23) & !is.numeric(ini.dist.23))
 {stop("Argument 'ini.dist.23' must be a numeric vector (default is NULL)")}  
 
if(dist[1]=="PE" & !is.null(ini.dist.12) & length(ini.dist.12)!=(length(cuts.12)+1))
 {stop("Incorrect number of parameters initialized for transition 12 (piecewise model)")} 
if(dist[2]=="PE" & !is.null(ini.dist.13) & length(ini.dist.13)!=(length(cuts.13)+1))
 {stop("Incorrect number of parameters initialized for transition 13 (piecewise model)")}
if(dist[3]=="PE" & !is.null(ini.dist.23) & length(ini.dist.23)!=(length(cuts.23)+1))
 {stop("Incorrect number of parameters initialized for transition 23 (piecewise model)")}
 
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
 
if( (dist[3]=="E" & is.null(cuts.23) & !is.null(ini.dist.23) & length(ini.dist.23)!=1) )
 {stop("Exponential distribution (transition 23) needs initialization of one parameter")} 
if( (dist[3]=="W" & is.null(cuts.23) & !is.null(ini.dist.23) & length(ini.dist.23)!=2) )
 {stop("Weibull distribution (transition 23) needs initialization of two parameters")}
if( (dist[3]=="WG" & is.null(cuts.23) & !is.null(ini.dist.23) & length(ini.dist.23)!=3) )
 {stop("Generalized Weibull distribution (transition 23) needs initialization of three parameters")}
 
 

if(!is.null(cov.12))
{
if ((!is.vector(cov.12) & !is.data.frame(cov.12) & !is.matrix(cov.12)) | !all(sapply(cov.12, is.numeric)))
 {stop("Argument 'cov.12' must be a numeric matrix or data.frame (default is NULL)")} 
if (nrow(data.frame(cov.12))!=length(times1))
 {stop("Argument 'cov.12' needs to have the same number of rows than 'times1'")}
if (sum(apply(sapply(data.frame(cov.12),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.12),is.na),1,sum))," individuals with missing values on 'cov.12' will be removed  from the analysis \n")
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

if(!is.null(cov.23))
{
if ((!is.vector(cov.23) & !is.data.frame(cov.23) & !is.matrix(cov.23)) | !all(sapply(cov.23,is.numeric)))
 {stop("Argument 'cov.23' must be a numeric matrix or data.frame (default is NULL)")} 
if (nrow(data.frame(cov.23))!=length(times2))
 {stop("Argument 'cov.23' needs to have the same number of rows than 'times2'")}
if (sum(apply(sapply(data.frame(cov.23),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.23),is.na),1,sum))," individuals with missing values on 'cov.23' will be removed from tha analysis \n")
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

 
 if(!(conf.int %in% c(TRUE,FALSE)))
 {stop("Argument 'conf.int' must be TRUE or FALSE (default is TRUE)")} 

 if(!is.null(precision))
 {
 if(!is.numeric(precision))
 {stop("Argument 'precision' must be numeric (default is 0)")} 
 if(precision<0)
 {stop("Argument 'precision' must be greater or equal to 0 (default is 0)")}
 }
 
  if(!(silent %in% c(TRUE,FALSE)))
 {stop("Argument 'silent' must be TRUE or FALSE (default is TRUE)")} 
 
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


if(dist[3]=="WG" | dist[3]=="W" | (dist[3]=="E" & is.null(cuts.23)))
 {
H23<-function(t,u,z,cuts) { exp(as.matrix(z) %*% coef23) * ( ((((1+(u/sigma23)^nu23))^(1/theta23))-1) - ((((1+(t/sigma23)^nu23))^(1/theta23))-1)) }
 log.h23<-function(t,z,cuts) { (as.matrix(z) %*% coef23) - log(theta23) + ((1/theta23)-1) * log1p((t/sigma23)^nu23) + log(nu23) + (nu23-1)*log(t) - nu23*log(sigma23) }
 }

if(dist[3]=="PE" & !is.null(cuts.23))
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
h23<-function(t,z,cuts) {
h<-rep(0,length(t))
for (i in (1:(length(cuts)-1)))
  {
  h<-h+(1*(t>=cuts[i])*(t<cuts[i+1]))* exp(as.matrix(z) %*% coef23)*(1/sigma23[i])
}
return(h)
rm(h)
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


#contributions to the log-likelihood
c1<-function(times1, z12, z13, cutimes12, cutimes13)
 {return( -H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13) ) }

c13<-function(times1, z12, z13, cutimes12, cutimes13)
 {return(log.h13(times1, z13, cutimes13) - H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13))}

c123<-function(times1, times2, z12, z13, z23, cutimes12, cutimes13, cutimes23)
 {return(log.h12(times1, z12, cutimes12) - H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13) + log.h23(times2, z23, cutimes23) - H23(times1, times2, z23, cutimes23))}

c12<-function(times1, times2, z12, z13, z23, cutimes12, cutimes13, cutimes23)
 {return(log.h12(times1, z12, cutimes12) - H12(times1, z12, cutimes12) - H13(times1, z13, cutimes13) - H23(times1, times2, z23, cutimes23))}

#chronological time
d1<-times1
d2<-times2

#missing data
.D <- cbind(d1, cov.12, cov.13, cov.23)
.na <- (!is.na(apply(.D, FUN=sum, MARGIN=1)))

#initialization of the parameters
if (is.null(cov.12)) {cov.12.mat <- cbind(rep(0, length(d1))); n.12 <- NULL} else { cov.12.mat <- cbind(cov.12); n.12 <- paste("covariate(s) on trans. 12:  num", 1:ncol(data.frame(cov.12))); if(!is.null(names.12)) {n.12 <- names.12} }

if (is.null(cov.13)) {cov.13.mat <- cbind(rep(0, length(d1))); n.13 <- NULL} else { cov.13.mat <- cbind(cov.13); n.13 <- paste("covariate(s) on trans. 13:  num", 1:ncol(data.frame(cov.13))); if(!is.null(names.13)) {n.13 <- names.13} }

if (is.null(cov.23)) {cov.23.mat <- cbind(rep(0, length(d1))); n.23 <- NULL} else { cov.23.mat <- cbind(cov.23); n.23 <- paste("covariate(s) on trans. 23:  num", 1:ncol(data.frame(cov.23))); if(!is.null(names.23)) {n.23 <- names.23} }

if (is.null(ini.dist.12)) {i.12.dist<-rep(0, 1*(dist[1]=="E" & is.null(cuts.12)) + 2*(dist[1]=="W") + 3*(dist[1]=="WG") + 1*(dist[1]=="PE" & !is.null(cuts.12))*(length(cuts.12)-1))}
 else {i.12.dist<-ini.dist.12}

if (is.null(ini.dist.13)) {i.13.dist<-rep(0, 1*(dist[2]=="E" & is.null(cuts.13)) + 2*(dist[2]=="W") + 3*(dist[2]=="WG") + 1*(dist[2]=="PE" & !is.null(cuts.13))*(length(cuts.13)-1))}
 else {i.13.dist<-ini.dist.13}

if (is.null(ini.dist.23)) {i.23.dist<-rep(0, 1*(dist[3]=="E" & is.null(cuts.23)) + 2*(dist[3]=="W") + 3*(dist[3]=="WG") + 1*(dist[3]=="PE" & !is.null(cuts.23))*(length(cuts.23)-1))}
 else {i.23.dist<-ini.dist.23}

if (!is.null(init.cov.12)) {i.12<-init.cov.12}
if (is.null(init.cov.12) & is.null(cov.12)) {i.12<-NULL}
if (is.null(init.cov.12) & !is.null(cov.12)) {i.12<-rep(0, ncol(data.frame(cov.12)))}

if (!is.null(init.cov.13)) {i.13<-init.cov.13}
if (is.null(init.cov.13) & is.null(cov.13)) {i.13<-NULL}
if (is.null(init.cov.13) & !is.null(cov.13)) {i.13<-rep(0, ncol(data.frame(cov.13)))}

if (!is.null(init.cov.23)) {i.23<-init.cov.23}
if (is.null(init.cov.23) & is.null(cov.23)) {i.23<-NULL}
if (is.null(init.cov.23) & !is.null(cov.23)) {i.23<-rep(0, ncol(data.frame(cov.23)))}

ini <- c(i.12.dist, i.13.dist, i.23.dist, i.12, i.13, i.23)

if (is.null(weights)) {w <- rep(1, length(d1))} else {w <- weights}

#parameters for contributions associated to each transition
.w1 <- w[(sequences==1 & .na)]
.d1.1 <- d1[(sequences==1 & .na)]
.c1.12 <- cov.12.mat[(sequences==1 & .na),]
.c1.13 <- cov.13.mat[(sequences==1 & .na),]

.w13 <- w[(sequences==13 & .na)]
.d13.1 <- d1[(sequences==13 & .na)]
.c13.12 <- cov.12.mat[(sequences==13 & .na),]
.c13.13 <- cov.13.mat[(sequences==13 & .na),]

.w123 <- w[(sequences==123 & .na)]
.d123.1 <- d1[(sequences==123 & .na)]
.d123.2 <- d2[(sequences==123 & .na)]
.c123.12 <- cov.12.mat[(sequences==123 & .na),]
.c123.13 <- cov.13.mat[(sequences==123 & .na),]
.c123.23 <- cov.23.mat[(sequences==123 & .na),]

.w12 <- w[(sequences==12 & .na)]
.d12.1 <- d1[(sequences==12 & .na)]
.d12.2 <- d2[(sequences==12 & .na)]
.c12.12 <- cov.12.mat[(sequences==12 & .na),]
.c12.13 <- cov.13.mat[(sequences==12 & .na),]
.c12.23 <- cov.23.mat[(sequences==12 & .na),]

#log-likelihood
logV <- function(x)
{
if (dist[1]=="E" & is.null(cuts.12)) {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", 1, inherits = TRUE); assign("theta12", 1, inherits = TRUE); i<-1}
if (dist[1]=="W") {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", exp(x[2]), inherits = TRUE); assign("theta12", 1, inherits = TRUE); i<-2}
if (dist[1]=="WG") {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", exp(x[2]), inherits = TRUE); assign("theta12", exp(x[3]), inherits = TRUE); i<-3}
if (dist[1]=="PE" & !is.null(cuts.12)) {assign("sigma12", exp(x[1:(length(cuts.12)-1)]), inherits = TRUE); i<-(length(cuts.12)-1)}

if (dist[2]=="E" & is.null(cuts.13)) {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", 1, inherits = TRUE); assign("theta13", 1, inherits = TRUE); i<-i+1}
if (dist[2]=="W") {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", exp(x[i+2]), inherits = TRUE); assign("theta13", 1, inherits = TRUE); i<-i+2}
if (dist[2]=="WG") {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", exp(x[i+2]), inherits = TRUE); assign("theta13", exp(x[i+3]), inherits = TRUE); i<-i+3}
if (dist[2]=="PE" & !is.null(cuts.13)) {assign("sigma13", exp(x[(i+1):(i+length(cuts.13)-1)]), inherits = TRUE); i<-(i+length(cuts.13)-1)}

if (dist[3]=="E" & is.null(cuts.23)) {assign("sigma23", exp(x[i+1]), inherits = TRUE); assign("nu23", 1, inherits = TRUE); assign("theta23", 1, inherits = TRUE); i<-i+1}
if (dist[3]=="W") {assign("sigma23", exp(x[i+1]), inherits = TRUE); assign("nu23", exp(x[i+2]), inherits = TRUE); assign("theta23", 1, inherits = TRUE); i<-i+2}
if (dist[3]=="WG") {assign("sigma23", exp(x[i+1]), inherits = TRUE); assign("nu23", exp(x[i+2]), inherits = TRUE); assign("theta23", exp(x[i+3]), inherits = TRUE); i<-i+3}
if (dist[3]=="PE" & !is.null(cuts.23)) {assign("sigma23", exp(x[(i+1):(i+length(cuts.23)-1)]), inherits = TRUE); i<-(i+length(cuts.23)-1)}

if (is.null(cov.12)) {assign("coef12", 0, inherits = TRUE)}
 else {assign("coef12", x[(i+1):(i+ncol(data.frame(cov.12)))], inherits = TRUE); i <-i+ncol(data.frame(cov.12))}
if (is.null(cov.13)) {assign("coef13", 0, inherits = TRUE)}
 else {assign("coef13", x[(i+1):(i+ncol(data.frame(cov.13)))], inherits = TRUE); i <-i+ncol(data.frame(cov.13))}
 if (is.null(cov.23)) {assign("coef23", 0, inherits = TRUE)}
 else {assign("coef23", x[(i+1):(i+ncol(data.frame(cov.23)))], inherits = TRUE); i <-i+ncol(data.frame(cov.23))}
 
return( -1*(
 sum( .w1 * c1(.d1.1, .c1.12, .c1.13, cuts.12, cuts.13) ) +
 sum( .w13 * c13(.d13.1, .c13.12, .c13.13, cuts.12, cuts.13) ) +
 sum( .w123 * c123(.d123.1, .d123.2, .c123.12, .c123.13, .c123.23, cuts.12, cuts.13, cuts.23) ) + 
 sum( .w12 * c12(.d12.1, .d12.2, .c12.12, .c12.13, .c12.23, cuts.12, cuts.13, cuts.23) ) ) )
}


#initialization of the parameters
if (is.null(cov.12)) {cov.12.mat <- cbind(rep(0, length(d1))); n.12 <- NULL} else { cov.12.mat <- cbind(cov.12); n.12 <- paste("covariate(s) on trans. 12: num", 1:ncol(data.frame(cov.12))); if(!is.null(names.12)) {n.12 <- names.12} }

if (is.null(cov.13)) {cov.13.mat <- cbind(rep(0, length(d1))); n.13 <- NULL} else { cov.13.mat <- cbind(cov.13); n.13 <- paste("covariate(s) on trans. 13: num", 1:ncol(data.frame(cov.13))); if(!is.null(names.13)) {n.13 <- names.13} }

if (is.null(cov.23)) {cov.23.mat <- cbind(rep(0, length(d1))); n.23 <- NULL} else { cov.23.mat <- cbind(cov.23); n.23 <- paste("covariate(s) on trans. 23: num", 1:ncol(data.frame(cov.23))); if(!is.null(names.23)) {n.23 <- names.23} }

if (is.null(ini.dist.12)) {i.12.dist <- rep(0, 1*(dist[1]=="E" & is.null(cuts.12)) + 2*(dist[1]=="W") + 3*(dist[1]=="WG") + 1*(dist[1]=="PE" & !is.null(cuts.12))*(length(cuts.12)-1))}  else {i.12.dist<-ini.dist.12}

if (is.null(ini.dist.13)) {i.13.dist <- rep(0, 1*(dist[2]=="E" & is.null(cuts.13)) + 2*(dist[2]=="W") + 3*(dist[2]=="WG") + 1*(dist[2]=="PE" & !is.null(cuts.13))*(length(cuts.13)-1))}  else {i.13.dist<-ini.dist.13}
 
if (is.null(ini.dist.23)) {i.23.dist<-rep(0, 1*(dist[3]=="E" & is.null(cuts.23)) + 2*(dist[3]=="W") + 3*(dist[3]=="WG") + 1*(dist[3]=="PE" & !is.null(cuts.23))*(length(cuts.23)-1))}  else {i.23.dist<-ini.dist.23}
 
if (!is.null(init.cov.12)) {i.12<-init.cov.12}
if (is.null(init.cov.12) & is.null(cov.12)) {i.12<-NULL}
if (is.null(init.cov.12) & !is.null(cov.12)) {i.12<-rep(0, ncol(data.frame(cov.12)))}
 
if (!is.null(init.cov.13)) {i.13<-init.cov.13}
if (is.null(init.cov.13) & is.null(cov.13)) {i.13<-NULL}
if (is.null(init.cov.13) & !is.null(cov.13)) {i.13<-rep(0, ncol(data.frame(cov.13)))}
 
if (!is.null(init.cov.23)) {i.23<-init.cov.23}
if (is.null(init.cov.23) & is.null(cov.23)) {i.23<-NULL}
if (is.null(init.cov.23) & !is.null(cov.23)) {i.23<-rep(0, ncol(data.frame(cov.23)))}

ini <- c(i.12.dist, i.13.dist, i.23.dist, i.12, i.13, i.23)

#first maximum likelihood optimization
n<-1
res<-tryCatch(optim(ini, logV, hessian=conf.int, control=list(maxit=100000)))

if(inherits(res, "error"))  {
warning("Maximum likelihood optimization fails to converge", "\n") } else { if(silent==FALSE) {warning(-1*res$value, "\n") }  

#further maximum likelihood optimizations
if(is.null(precision)) {delta <- 10^(-6)} else {delta<-precision}


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
  names(table.res)<-c("Estimate", "Std.Error", "t.value", "Pr(>|t|)")
  table.covariance<-res$hessian
  }
  else {
  table.res <- data.frame(Estimate = round(res$par, 4) )
  table.covariance<-NULL
  warning("Hessian matrix not defined", "\n")
  } #end else for hessian matrix condition
}

if (conf.int==FALSE) {
table.res <- data.frame(Coef = round(res$par, 4) ) 
table.covariance<-NULL
}

if (dist[1]=="E" & is.null(cuts.12))  {  lab12<-c("log(sigma) on trans. 12")}
if (dist[1]=="W" & is.null(cuts.12))  {  lab12<-c("log(sigma) on trans. 12", "log(nu) on trans. 12")}
if (dist[1]=="WG" & is.null(cuts.12)) {  lab12<-c("log(sigma) on trans. 12", "log(nu) on trans. 12", "log(theta) on trans. 12")}
if (dist[1]=="PE" & !is.null(cuts.12)) { lab12<-rep("",length(cuts.12)-1)
 for (i in (1:(length(cuts.12)-1)))
  {
  lab12[i]<-paste("log(sigma) on trans. 12, interval [", cuts.12[i], ";", cuts.12[i+1], "[", sep="")
  }
 }

if (dist[2]=="E" & is.null(cuts.13))  { lab13<-c("log(sigma) on trans. 13")}
if (dist[2]=="W" & is.null(cuts.13))  { lab13<-c("log(sigma) on trans. 13", "log(nu) on trans. 13")}
if (dist[2]=="WG" & is.null(cuts.13)) { lab13<-c("log(sigma) on trans. 13", "log(nu) on trans. 13", "log(theta) on trans. 13")}
if (dist[2]=="PE" & !is.null(cuts.13)) { lab13<-rep("",length(cuts.13)-1)
 for (i in (1:(length(cuts.13)-1)))
  {
  lab13[i]<-paste("log(sigma) on trans. 13, interval [", cuts.13[i], ";", cuts.13[i+1], "[", sep="")
  }
 }

if (dist[3]=="E" & is.null(cuts.23))  {  lab23<-c("log(sigma) on trans. 23")}
if (dist[3]=="W" & is.null(cuts.23))  {  lab23<-c("log(sigma) on trans. 23", "log(nu) on trans. 23")}
if (dist[3]=="WG" & is.null(cuts.23)) {  lab23<-c("log(sigma) on trans. 23", "log(nu) on trans. 23", "log(theta) on trans. 23")}
if (dist[3]=="PE" & !is.null(cuts.23)) { lab23<-rep("",length(cuts.23)-1)
 for (i in (1(length(cuts.23)-1)))
  {
  lab23[i]<-paste("log(sigma) on trans. 23, interval [", cuts.23[i], ";", cuts.23[i+1], "[", sep="")
  }
 }

lab<-c(lab12, lab13, lab23, n.12, n.13, n.23)

rownames(table.res) <- paste(1:length(lab), lab)

warning("\n Number of data rows:",nrow(.D))
warning("Number of data rows with missing values (deleted):",nrow(.D)-sum(.na),"\n")

return(list(
object="m3 (3-state time-inhomogeneous markov model)",
dist=dist,
cuts.12=cuts.12,
cuts.13=cuts.13,
cuts.23=cuts.23,
covariates=c( max(0, length(n.12)), max(0, length(n.13)), max(0, length(n.23)) ),
table=table.res,
cov.matrix=table.covariance,
LogLik=(-1*res$value),
AIC=2*length(res$par)-2*(-1*res$value)))
  } 
  } 
}
