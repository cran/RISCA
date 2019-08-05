
mixture.2states <-function(times, sequences, weights=NULL, dist, cuts.12=NULL, cuts.13=NULL, ini.dist.12=NULL, ini.dist.13=NULL, cov.12=NULL, init.cov.12=NULL, names.12=NULL, cov.13=NULL, init.cov.13=NULL, names.13=NULL, cov.p=NULL, init.cov.p=NULL, names.p=NULL, init.intercept.p=NULL, conf.int=TRUE, silent=TRUE, precision=10^(-6))
{

#check conditions
if (missing(times)) 
        stop("Argument 'times' is missing with no default")
		
if (missing(sequences)) 
        stop("Argument 'sequences' is missing with no default")
		
if (missing(dist)) 
        stop("Argument 'dist' is missing with no default")
        
if (!is.vector(times) | !is.numeric(times))
        stop("Argument 'times' must be a numeric vector")
		
if (min(times,na.rm=T)<0)
		stop("Negative values for 'times' are not allowed")
		
if (is.na(min(times)))
		warning("individuals with missing values for 'times' will be removed from the analysis \n")
		
if (!is.vector(sequences) | !is.numeric(sequences) | (min(names(table(sequences)) %in% c(1,12,13))==0) )
        stop("Argument 'sequences' must be a numeric vector with values 1, 12, or 13")
		
if (min( c(1,12,13) %in% names(table(sequences)))==0)
        warning("all sequencess (1, 12, and 13) are not present \n ")
        
if (min(length(times), length(sequences)) != max(length(times), length(sequences)))
        stop("Arguments 'times' and 'sequences' need to have the same number of rows")
        
if(!is.null(weights))
{
	if (!is.vector(weights) | !is.numeric(weights))
        stop("Argument 'weights' must be a numeric vector")
		
	if (min(weights,na.rm=T)<0)
		stop("Negative values for 'weights' are not allowed")
		
	if (is.na(min(weights)))
		warning("individuals with missing values for 'weights' will be removed from the analysis \n")
}
 
if(length(dist)!=2)
 {stop("Argument 'dist' have to contain 2 values")} 
 
 if(!(dist[1] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 12")}
 
if(!(dist[2] %in% c("PE","E","W","WG")))
 {stop("Argument 'dist': incorrect distribution for transition 13")}
 
if(dist[1]!="PE" & (!is.null(cuts.12)))
 {stop("Arguments 'cuts.12' is only allowed for piecewise exponential distribution (PE for the first argument in 'dist')")}

if(dist[2]!="PE" & (!is.null(cuts.13)))
 {stop("Arguments 'cuts.13' is only allowed for piecewise exponential distribution (PE for the second argument in 'dist')")}

if(dist[1]=="PE" & !is.null(cuts.12))
 {
 if (!all(is.numeric(cuts.12)) | !all(!is.na(cuts.12)) | !all(cuts.12>0) | !all(is.finite(cuts.12)) | is.unsorted(cuts.12)) 
 {stop("Arguments 'cuts.12' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
 
if(dist[1]=="PE" & !is.null(cuts.12))
{
 if (max(cuts.12)>=max(times,na.rm=T)) 
 {stop("Arguments 'cuts.12': check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for 'times')")}
}

if(dist[2]=="PE" & !is.null(cuts.13))
 {
 if (!all(is.numeric(cuts.13)) | !all(!is.na(cuts.13)) | !all(cuts.13>0) | !all(is.finite(cuts.13)) | is.unsorted(cuts.13)) 
 {stop("Arguments 'cuts.13' must be a sorted vector with only positive and finite numeric values (internal timepoints)")}
 }
 
if(dist[2]=="PE" & !is.null(cuts.13))
{
 if (max(cuts.13)>=max(times,na.rm=T)) 
 {stop("Arguments 'cuts.13': check internal timepoints or time units (last internal timepoint is greater or equal to the maximum value for 'times')")}
}

if(!is.null(ini.dist.12) & !is.numeric(ini.dist.12))
 {stop("Argument 'ini.dist.12' must be a numeric vector (default is NULL)")} 
 
if(!is.null(ini.dist.13) & !is.numeric(ini.dist.13))
 {stop("Argument 'ini.dist.13' must be a numeric vector (default is NULL)")}  
 
if(dist[1]=="PE" & !is.null(ini.dist.12) & length(ini.dist.12)!=(length(cuts.12)+1))
 {stop("Incorrect number of parameters initialized for transition 12 (piecewise model)")}
 
if(dist[2]=="PE" & !is.null(ini.dist.13) & length(ini.dist.13)!=(length(cuts.13)+1))
 {stop("Incorrect number of parameters initialized for transition 13 (piecewise model)")}
 
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
 
if(!is.null(cov.12))
{
if ((!is.vector(cov.12) & !is.data.frame(cov.12) & !is.matrix(cov.12)) | !all(sapply(cov.12,is.numeric)))
 {stop("Argument 'cov.12' must be a numeric matrix or data.frame (default is NULL)")} 

if (nrow(data.frame(cov.12))!=length(times))
 {stop("Argument 'cov.12' needs to have the same number of rows than 'times'")}

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
 
if (nrow(data.frame(cov.13))!=length(times))
 {stop("Argument 'cov.13' needs to have the same number of rows than 'times'")}
 
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



if(!is.null(cov.p))
{
if ((!is.vector(cov.p) & !is.data.frame(cov.p) & !is.matrix(cov.p)) | !all(sapply(cov.p,is.numeric)))
 {stop("Argument 'cov.p' must be a numeric matrix or data.frame (default is NULL)")} 

if (nrow(data.frame(cov.p))!=length(times))
 {stop("Argument 'cov.p' needs to have the same number of rows than 'times'")}

if (sum(apply(sapply(data.frame(cov.p),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.p),is.na),1,sum))," individuals with missing values on 'cov.p' will be removed from the analysis\n")
		
if(!is.null(init.cov.p))
	{
	if (!is.numeric(init.cov.p))
	{stop("Argument 'init.cov.p' must be a numeric vector (default is NULL)")}
	
	if (ncol(data.frame(cov.p))!=length(init.cov.p))
	{stop("Argument 'init.cov.p' needs to have the same length than number of columns of 'cov.p'")}
	}
	
if (!is.null(names.p))
	{
	if (!is.character(names.p))
	{stop("Argument 'names.p' must be a character vector (default is NULL)")}
	
	if (ncol(data.frame(cov.p))!=length(names.p))
	{stop("Argument 'names.p' needs to have the same length than number of columns of 'cov.p'")}
	}
}


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
	

coef12<-NULL
sigma12<-NULL
nu12<-NULL
theta12<-NULL
coef13<-NULL
sigma13<-NULL
nu13<-NULL
theta13<-NULL
a0<-NULL
a1<-NULL

#sojourn time distributions

if(dist[1]=="WG" | dist[1]=="W" | (dist[1]=="E" & is.null(cuts.12)))
 {
 H12<-function(times,z,cuts) { exp(as.matrix(z) %*% coef12) * ((((1+(times/sigma12)^nu12))^(1/theta12))-1) }
 
 log.h12<-function(times,z,cuts) { (as.matrix(z) %*% coef12) - log(theta12) + ((1/theta12)-1) * log1p((times/sigma12)^nu12) + log(nu12) + (nu12-1)*log(times) - nu12*log(sigma12) }
 }

if(dist[1]=="PE" & !is.null(cuts.12))
 {
cuts.12 <- sort(cuts.12)
if ((cuts.12[1] <= 0) || (cuts.12[length(cuts.12)] == Inf)) 
   stop("'cuts.12' must be positive and finite.")
cuts.12 <- c(0, cuts.12, Inf)

H12<-function(times,z,cuts) {
 H<-rep(0,length(times))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(times>=cuts[i]))*exp(as.matrix(z) %*% coef12)*((pmin(cuts[i+1],times)-cuts[i])/sigma12[i])
  }
return(H)
rm(H)
 }
 
log.h12<-function(times,z,cuts) {
 log.h<-rep(0,length(times))
 for (i in (1:(length(cuts)-1)))
  {
  log.h<-log.h+(1*(times>=cuts[i])*(times<cuts[i+1]))*(as.matrix(z) %*% coef12-log(sigma12[i]))
  }
 return(log.h)
 rm(log.h)
}
}

if(dist[2]=="WG" | dist[2]=="W" | (dist[2]=="E" & is.null(cuts.13)))
 {
 H13<-function(times,z,cuts) { exp(as.matrix(z) %*% coef13) * ((((1+(times/sigma13)^nu13))^(1/theta13))-1) }
 
 log.h13<-function(times,z,cuts) { (as.matrix(z) %*% coef13) - log(theta13) + ((1/theta13)-1) * log1p((times/sigma13)^nu13) + log(nu13) + (nu13-1)*log(times) - nu13*log(sigma13) }
 }

if(dist[2]=="PE" & !is.null(cuts.13))
 {
cuts.13 <- sort(cuts.13)
if ((cuts.13[1] <= 0) || (cuts.13[length(cuts.13)] == Inf)) 
   stop("'cuts.13' must be positive and finite.")
cuts.13 <- c(0, cuts.13, Inf)

H13<-function(times,z,cuts) {
 H<-rep(0,length(times))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(times>=cuts[i]))*exp(as.matrix(z) %*% coef13)*((pmin(cuts[i+1],times)-cuts[i])/sigma13[i])
  }
return(H)
rm(H)
 }
 
log.h13<-function(times,z,cuts) {
 log.h<-rep(0,length(times))
 for (i in (1:(length(cuts)-1)))
  {
  log.h<-log.h+(1*(times>=cuts[i])*(times<cuts[i+1]))*(as.matrix(z) %*% coef13-log(sigma13[i]))
  }
 return(log.h)
 rm(log.h)
}
}

p2<-function(z)
	{exp(a0 + as.matrix(z) %*% as.vector(a1))/(exp(a0 + as.matrix(z) %*% as.vector(a1))+1)}

log.p2<-function(z)
	{a0 + as.matrix(z) %*% as.vector(a1) - log(exp(a0 + as.matrix(z) %*% as.vector(a1))+1)}

p3<-function(z)
	{1/(exp(a0 + as.matrix(z) %*% as.vector(a1))+1)}

log.p3<-function(z)
	{ -1 * log(exp(a0 + as.matrix(z) %*% as.vector(a1))+1)}


#contributions to the log-likelihood
c1<-function(times, z12, z13, zp, cut12, cut13)
 {return( log ( p2(zp) * exp(-H12(times, z12, cut12)) + p3(zp) * exp(-H13(times, z13, cut13)) ) ) }

c12<-function(times, z12, zp, cut12)
 {return(log.p2(zp) + log.h12(times, z12, cut12) - H12(times, z12, cut12)  )}

c13<-function(times, z13, zp, cut13)
 {return(log.p3(zp) + log.h13(times, z13, cut13) - H13(times, z13, cut13)  )}

#missing data
.D <- cbind(times, cov.12, cov.13, cov.p, weights)
.na <- (!is.na(apply(.D, FUN="sum", MARGIN=1)))

#initialization of the parameters
if (is.null(cov.12)) {cov.12.mat <- cbind(rep(0, length(times))); n.12 <- NULL} else { cov.12.mat <- cbind(cov.12); n.12 <- paste("covariate(s) on trans. 12:  num", 1:ncol(data.frame(cov.12))); if(!is.null(names.12)) {n.12 <- names.12} }

if (is.null(cov.13)) {cov.13.mat <- cbind(rep(0, length(times))); n.13 <- NULL} else { cov.13.mat <- cbind(cov.13); n.13 <- paste("covariate(s) on trans. 13:  num", 1:ncol(data.frame(cov.13))); if(!is.null(names.13)) {n.13 <- names.13} }

if (is.null(cov.p)) {cov.p.mat <- cbind(rep(0, length(times))); n.p <- NULL} else { cov.p.mat <- cbind(cov.p); n.p <- paste("covariate(s) associated with P(X=2):  num", 1:ncol(data.frame(cov.p))); if(!is.null(names.p)) {n.p <- names.p} }

if (is.null(ini.dist.12)) {i.12.dist<-rep(0, 1*(dist[1]=="E" & is.null(cuts.12)) + 2*(dist[1]=="W") + 3*(dist[1]=="WG") + 1*(dist[1]=="PE" & !is.null(cuts.12))*(length(cuts.12)-1))}
 else {i.12.dist<-ini.dist.12}

if (is.null(ini.dist.13)) {i.13.dist<-rep(0, 1*(dist[2]=="E" & is.null(cuts.13)) + 2*(dist[2]=="W") + 3*(dist[2]=="WG") + 1*(dist[2]=="PE" & !is.null(cuts.13))*(length(cuts.13)-1))}
 else {i.13.dist<-ini.dist.13}
 
if (!is.null(init.cov.12)) {i.12<-init.cov.12}
if (is.null(init.cov.12) & is.null(cov.12)) {i.12<-NULL}
if (is.null(init.cov.12) & !is.null(cov.12)) {i.12<-rep(0, ncol(data.frame(cov.12)))}
 
if (!is.null(init.cov.13)) {i.13<-init.cov.13}
if (is.null(init.cov.13) & is.null(cov.13)) {i.13<-NULL}
if (is.null(init.cov.13) & !is.null(cov.13)) {i.13<-rep(0, ncol(data.frame(cov.13)))}
 
if (!is.null(init.cov.p)) {i.p<-init.cov.p}
if (is.null(init.cov.p) & is.null(cov.p)) {i.p<-NULL}
if (is.null(init.cov.p) & !is.null(cov.p)) {i.p<-rep(0, ncol(data.frame(cov.p)))}

if (!is.null(init.intercept.p)) {int.p<-init.intercept.p} else {int.p<-0}

ini <- c(i.12.dist, i.13.dist, i.12, i.13, int.p, i.p)

if (is.null(weights)) {w <- rep(1, length(times))} else {w <- weights }

#parameters for contributions associated to each transition
.w1 <- w[(sequences==1 & .na)]
.t1 <- times[(sequences==1 & .na)]
.c1.12 <- cov.12.mat[(sequences==1 & .na),]
.c1.13 <- cov.13.mat[(sequences==1 & .na),]
.c1.p <- cov.p.mat[(sequences==1 & .na),]

.w12 <- w[(sequences==12 & .na)]
.t12 <- times[(sequences==12 & .na)]
.c12.12 <- cov.12.mat[(sequences==12 & .na),]
.c12.p <- cov.p.mat[(sequences==12 & .na),]

.w13 <- w[(sequences==13 & .na)]
.t13 <- times[(sequences==13 & .na)]
.c13.13 <- cov.13.mat[(sequences==13 & .na),]
.c13.p <- cov.p.mat[(sequences==13 & .na),]

#log-likelihood
logV<-function(x)
{
if (dist[1]=="E" & is.null(cuts.12)) {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", 1, inherits = TRUE); assign("theta12", 1, inherits = TRUE); i<-1}
if (dist[1]=="W") {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", exp(x[2]), inherits = TRUE); assign("theta12", 1, inherits = TRUE); i<-2}
if (dist[1]=="WG") {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", exp(x[2]), inherits = TRUE); assign("theta12", exp(x[3]), inherits = TRUE); i<-3}
if (dist[1]=="PE" & !is.null(cuts.12)) {assign("sigma12", exp(x[1:(length(cuts.12)-1)]), inherits = TRUE); i<-(length(cuts.12)-1)}

if (dist[2]=="E" & is.null(cuts.13)) {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", 1, inherits = TRUE); assign("theta13", 1, inherits = TRUE); i<-i+1}
if (dist[2]=="W") {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", exp(x[i+2]), inherits = TRUE); assign("theta13", 1, inherits = TRUE); i<-i+2}
if (dist[2]=="WG") {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", exp(x[i+2]), inherits = TRUE); assign("theta13", exp(x[i+3]), inherits = TRUE); i<-i+3}
if (dist[2]=="PE" & !is.null(cuts.13)) {assign("sigma13", exp(x[(i+1):(i+length(cuts.13)-1)]), inherits = TRUE); i<-(i+length(cuts.13)-1)}

if (is.null(cov.12)) {assign("coef12", 0, inherits = TRUE)}
 else {assign("coef12", x[(i+1):(i+ncol(data.frame(cov.12)))], inherits = TRUE); i <-i+ncol(data.frame(cov.12))}

if (is.null(cov.13)) {assign("coef13", 0, inherits = TRUE)}
 else {assign("coef13", x[(i+1):(i+ncol(data.frame(cov.13)))], inherits = TRUE); i <-i+ncol(data.frame(cov.13))}
 
assign("a0", x[i+1], inherits = TRUE); i<-i+1

if (is.null(cov.p)) {assign("a1", 0, inherits = TRUE)}
 else {assign("a1", x[(i+1):(i+ncol(data.frame(cov.p)))], inherits = TRUE)}
 
return( -1*(
 sum( .w1 * c1(.t1, .c1.12, .c1.13, .c1.p, cuts.12, cuts.13) ) +
 sum( .w12 * c12(.t12, .c12.12, .c12.p, cuts.12) ) +
 sum( .w13 * c13(.t13, .c13.13, .c13.p, cuts.13) ) ) )
}

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
  table.covariance<-solve(res$hessian)
  }
  else {
  table.res <- data.frame(Estimate = round(res$par, 4) )
  table.covariance<-NULL
  warning("\n Hessian matrix not defined", "\n")
  } #end else for hessian matrix condition
}

if (conf.int==FALSE) {
table.res <- data.frame(Coef = round(res$par, 4) ) 
table.covariance<-NULL
}

if (dist[1]=="E" & is.null(cuts.12))  { lab12<-c("log(sigma) on trans. 12")}
if (dist[1]=="W" & is.null(cuts.12))  { lab12<-c("log(sigma) on trans. 12", "log(nu) on trans. 12")}
if (dist[1]=="WG" & is.null(cuts.12)) { lab12<-c("log(sigma) on trans. 12", "log(nu) on trans. 12", "log(theta) on trans. 12")}
if (dist[1]=="PE" & !is.null(cuts.12)) {
 lab12<-rep("",length(cuts.12)-1)
 for (i in (1:(length(cuts.12)-1)))
  {
  lab12[i]<-paste("log(sigma) on trans. 12, interval [",cuts.12[i],";",cuts.12[i+1],"[",sep="")
  }
 }

if (dist[2]=="E" & is.null(cuts.13))  { lab13<-c("log(sigma) on trans. 13")}
if (dist[2]=="W" & is.null(cuts.13))  { lab13<-c("log(sigma) on trans. 13", "log(nu) on trans. 13")}
if (dist[2]=="WG" & is.null(cuts.13)) { lab13<-c("log(sigma) on trans. 13", "log(nu) on trans. 13", "log(theta) on trans. 13")}
if (dist[2]=="PE" & !is.null(cuts.13)) {
 lab13<-rep("",length(cuts.13)-1)
 for (i in (1:(length(cuts.13)-1)))
  {
  lab13[i]<-paste("log(sigma) on trans. 13, interval [",cuts.13[i],";",cuts.13[i+1],"[",sep="")
  }
 }

lab<-c(lab12, lab13, n.12, n.13, "intercept of the logit of P(X=2)", n.p)

rownames(table.res) <- paste(1:length(lab), lab)

warning("\n Number of data rows:",nrow(.D))
warning("Number of data rows with missing values (deleted):",nrow(.D)-sum(.na),"\n")

return(list(
object="mixture.2states",
dist=dist,
cuts.12=cuts.12,
cuts.13=cuts.13,
covariates=c( max(0, length(n.12)), max(0, length(n.13)), max(0, length(n.p)) ),
table=table.res,
cov.matrix=table.covariance,
LogLik=(-1*res$value),
AIC=2*length(res$par)-2*(-1*res$value)))
  } 
  } 
}

