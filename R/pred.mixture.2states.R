
pred.mixture.2states <- function(model, failure, times, cov.12=NULL, cov.13=NULL, cov.p=NULL)
{
#check conditions

if (missing(model)) 
        stop("Argument 'model' is missing with no default")
		
if (model$object!="mixture.2states") 
        stop("'model' have to be the list results of the function 'mixture.2states'")		

if (missing(failure)) 
        stop("Argument 'failure' is missing with no default")
		
if (!is.numeric(failure) | (failure!=2 & failure!=3)) 
        stop("Argument 'failure' must be a numeric value equals to 2 or 3")

if (missing(times)) 
        stop("Argument 'times' is missing with no default")
		
if (!is.numeric(times))
        stop("Argument 'times' must be a numeric vector")

if (min(times,na.rm=TRUE)<0)
		stop("Negative values for 'times' are not allowed")

if (model$covariates[1]!=max(0, dim(cbind(cov.12))[2]))
		stop("The number of covariates related to the time-to-failure #2 does not correspond with the 'model'")

if (model$covariates[2]!=max(0, dim(cbind(cov.13))[2]))
		stop("The number of covariates related to the time-to-failure #3 does not correspond with the 'model'")
		
if (model$covariates[3]!=max(0, dim(cbind(cov.p))[2]))
		stop("The number of covariates related to P(X=2) does not correspond with the 'model'")
		
if(!is.null(cov.12))
{
if ((!is.vector(cov.12) & !is.data.frame(cov.12) & !is.matrix(cov.12)) | !all(sapply(cov.12, is.numeric)))
 {stop("Argument 'cov.12' must be a numeric matrix or data.frame (default is NULL)")} 

if (sum(apply(sapply(data.frame(cov.12),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.12),is.na),1,sum))," individuals with missing values on 'cov.12' will be removed \n")
}

if(!is.null(cov.13))
{
if ((!is.vector(cov.13) & !is.data.frame(cov.13) & !is.matrix(cov.13)) | !all(sapply(cov.13,is.numeric)))
 {stop("Argument 'cov.13' must be a numeric matrix or data.frame (default is NULL)")}
 
if (sum(apply(sapply(data.frame(cov.13),is.na),1,sum))>0)
		warning(sum(apply(sapply(data.frame(cov.13),is.na),1,sum))," individuals with missing values on 'cov.13' will be removed \n")
}

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

dist <- model$dist
cuts.12 <- model$cuts.12
cuts.13 <- model$cuts.13

if(dist[1]=="WG" | dist[1]=="W" | (dist[1]=="E" & is.null(cuts.12)))
{
 H12<-function(t,z,cuts) { exp(as.matrix(z) %*% as.vector(coef12)) * ((((1+(t/sigma12)^nu12))^(1/theta12))-1) }
}

if(dist[1]=="PE" & !is.null(cuts.12))
{
cuts.12 <- sort(cuts.12)
if ((cuts.12[1] <= 0) || (cuts.12[length(cuts.12)] == Inf))
   stop("'cuts.12' must be positive and finite.")
cuts.12 <- c(0, cuts.12, Inf)

H12<-function(t,z,cuts) {
 H<-rep(0,length(t))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(t>=cuts[i]))*exp(as.matrix(z) %*% as.vector(coef12))*((pmin(cuts[i+1],t)-cuts[i])/sigma12[i])
  }
return(H)
rm(H)
 }
}

if(dist[2]=="WG" | dist[2]=="W" | (dist[2]=="E" & is.null(cuts.13)))
{
 H13<-function(t,z,cuts) { exp(as.matrix(z) %*% as.vector(coef13)) * ((((1+(t/sigma13)^nu13))^(1/theta13))-1) }
}

if(dist[2]=="PE" & !is.null(cuts.13))
{
cuts.13 <- sort(cuts.13)
if ((cuts.13[1] <= 0) || (cuts.13[length(cuts.13)] == Inf)) 
   stop("'cuts.13' must be positive and finite.")
cuts.13 <- c(0, cuts.13, Inf)

H13<-function(t,z,cuts) {
 H<-rep(0,length(t))
for (i in (1:(length(cuts)-1)))
  {
  H<-H+(1*(t>=cuts[i]))*exp(as.matrix(z) %*% as.vector(coef13))*((pmin(cuts[i+1],t)-cuts[i])/sigma13[i])
  }
return(H)
rm(H)
 }
}

p2<-function(z)
	{exp(a0 + as.matrix(z) %*% as.vector(a1))/(exp(a0 + as.matrix(z) %*% as.vector(a1))+1)}

#CIF definition
if(failure==2)
{
cif<-function(t, z12, z13, zp, cut12, cut13)
 {return( p2(zp) * (1-exp(-H12(t, z12, cut12))) ) }
}

if(failure==3)
{
cif<-function(t, z12, z13, zp, cut12, cut13)
 {return( (1-p2(zp)) * (1-exp(-H13(t, z13, cut13))) ) }
}

#missing data
.D <- cbind(cov.12, cov.13, cov.p)
.na <- (!is.na(apply(.D, FUN="sum", MARGIN=1)))

nb.cov.max <- max(1, dim(cbind(cov.12))[1], dim(cbind(cov.13))[1], dim(cbind(cov.p))[1])

#initialization of the parameters
if (is.null(cov.12)) { cov.12.mat <- cbind(rep(0, nb.cov.max)) } else { cov.12.mat <- cbind(cov.12) }

if (is.null(cov.13)) { cov.13.mat <- cbind(rep(0, nb.cov.max)) } else { cov.13.mat <- cbind(cov.13) }

if (is.null(cov.p))  { cov.p.mat  <- cbind(rep(0, nb.cov.max)) } else { cov.p.mat  <- cbind(cov.p)  }

x <- model$table[,1]

if (dist[1]=="E" & is.null(cuts.12)) {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", 1, inherits = TRUE); assign("theta12", 1, inherits = TRUE); i<-1}

if (dist[1]=="W") {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", exp(x[2]), inherits = TRUE); assign("theta12", 1, inherits = TRUE); i<-2}

if (dist[1]=="WG") {assign("sigma12", exp(x[1]), inherits = TRUE); assign("nu12", exp(x[2]), inherits = TRUE); assign("theta12", exp(x[3]), inherits = TRUE); i<-3}

if (dist[1]=="PE" & !is.null(cuts.12)) {assign("sigma12", exp(x[1:(length(cuts.12)-1)]), inherits = TRUE); i<-(length(cuts.12)-1)}

if (dist[2]=="E" & is.null(cuts.13)) {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", 1, inherits = TRUE); assign("theta13", 1, inherits = TRUE); i<-i+1}

if (dist[2]=="W") {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", exp(x[i+2]), inherits = TRUE); assign("theta13", 1, inherits = TRUE); i<-i+2}

if (dist[2]=="WG") {assign("sigma13", exp(x[i+1]), inherits = TRUE); assign("nu13", exp(x[i+2]), inherits = TRUE); assign("theta13", exp(x[i+3]), inherits = TRUE); i<-i+3}

if (dist[2]=="PE" & !is.null(cuts.13)) {assign("sigma13", exp(x[(i+1):(i+length(cuts.13)-1)]), inherits = TRUE); i<-(i+length(cuts.13)-1)}

if (is.null(cov.12)) {assign("coef12", 0, inherits = TRUE)} else {assign("coef12", x[(i+1):(i+ncol(data.frame(cov.12)))], inherits = TRUE); i <-i+ncol(data.frame(cov.12))}

if (is.null(cov.13)) {assign("coef13", 0, inherits = TRUE)} else {assign("coef13", x[(i+1):(i+ncol(data.frame(cov.13)))], inherits = TRUE); i <-i+ncol(data.frame(cov.13))}
 
assign("a0", x[i+1], inherits = TRUE); i<-i+1

if (is.null(cov.p)) {assign("a1", 0, inherits = TRUE)} else {assign("a1", x[(i+1):(i+ncol(data.frame(cov.p)))], inherits = TRUE); i <-i+ncol(data.frame(cov.p))}

k<-1
.cif<-matrix(-99, ncol=length(times), nrow= nb.cov.max)
for (j in sort(unique(times)))
{
.cif[,k] <- cif(j, cov.12.mat, cov.13.mat, cov.p.mat, cuts.12, cuts.13)
k<-k+1
}

return(list(times=times, cif=.cif))
}

