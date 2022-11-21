\name{pred.mixture.2states}
\alias{pred.mixture.2states}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative Incidence Function Form Horizontal Mixture Model With Two Competing Events
}
\description{
This function allows to estimate a cumulative incidence function (CIF) from an horizontal mixture model with two competing events, i.e. the results obtained from the function \code{mixture.2states}.
}
\usage{
pred.mixture.2states(model, failure, times, cov.12=NULL, cov.13=NULL, cov.p=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A list obtained by using the function \code{mixture.2states}.}
  \item{failure}{A numeric value for identifying the event for which the CIF has to be computed. Two possible values are allowed: 2 (for the CIF related to X=2) and 3 (for the CIF related to X=3).}
  \item{times}{A numeric vector with positive values related to the times for which the CIF has to be computed.}
 \item{cov.12}{A vector, matrix or data frame in which to look for variables related to the time from X=1 to X=2 with which to predict the CIF.}
 \item{cov.13}{A vector, matrix or data frame in which to look for variables related to the time from X=1 to X=3 with which to predict the CIF.}
 \item{cov.p}{A vector, matrix or data frame in which to look for variables related to the probability P(X=2).}
 }
 
\details{
The covariates has to be identical than the ones included in the mixture model declared in the argument \code{model}. More precisely, the columns of \code{cov.12}, \code{cov.13} and \code{cov.p} must correspond to the same variables.
}

\value{
\item{times}{A numeric vector with the times for which the CIF has to be computed.}
\item{cif}{A matrix with the predicted CIF for the \code{times} in columns and the individuals in rows.}
}

\references{
Trebern-Launay et al. Horizontal mixture model for competing risks: a method used in waitlisted renal transplant candidates. European Journal of Epidemiology. 33(3):275-286, 2018. <doi: 10.1007/s10654-017-0322-3>.
}

\author{
Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>
}

\examples{
# import the observed data
# X=1 corresponds to initial state with a functioning graft,
# X=2 to acute rejection episode (transient state), 
# X=3 to return to dialysis, X=4 to death with a functioning graft

data(dataDIVAT1)

# A subgroup analysis to reduce the time needed for this example

dataDIVAT1$id<-c(1:nrow(dataDIVAT1))
set.seed(2)
d2<-dataDIVAT1[dataDIVAT1$id \%in\% sample(dataDIVAT1$id, 300, replace = FALSE),]

# Data-management: two competing events
# the patient death is now X=2
# the return in dialysis is now X=3

d2$time<-NA
d2$time[d2$trajectory==1]<-d2$time1[d2$trajectory==1]

d2$time[d2$trajectory==12]<-d2$time2[d2$trajectory==12]
d2$trajectory[d2$trajectory==12]<-1

d2$time[d2$trajectory==13]<-d2$time1[d2$trajectory==13]

d2$time[d2$trajectory==123]<-d2$time2[d2$trajectory==123]
d2$trajectory[d2$trajectory==123]<-13

d2$time[d2$trajectory==14]<-d2$time1[d2$trajectory==14]

d2$time[d2$trajectory==124]<-d2$time2[d2$trajectory==124]
d2$trajectory[d2$trajectory==124]<-14

d2$trajectory[d2$trajectory==14]<-12

table(d2$trajectory)

# Univariable horizontal mixture model one binary explicative variable
# z is 1 if delayed graft function and 0 otherwise

mm2.model <- mixture.2states(times=d2$time, sequences=d2$trajectory,
  weights=NULL, dist=c("E","W"), cuts.12=NULL, cuts.13=NULL, 
  ini.dist.12=c(9.28), ini.dist.13=c(9.92, -0.23), 
  cov.12=d2$z, init.cov.12=0.84, names.12="beta_12",
  cov.13=d2$z, init.cov.13=0.76, names.13="beta_13",
  cov.p=NULL, init.cov.p=NULL, names.p=NULL, init.intercept.p=-0.75,
  conf.int=TRUE, silent=FALSE)

cif2.mm2 <- pred.mixture.2states(mm2.model, failure=2, times=seq(0, 4000, by=30),
 cov.12=c(0,1), cov.13=c(0,1), cov.p=NULL)

plot(cif2.mm2$times/365.25, cif2.mm2$cif[1,], col = 1, type="l", lty = 1,
 ylim=c(0,1), lwd =2, ylab="Cumulative Incidence Function",
 xlab="Times (years)", main="", xlim=c(0, 11), legend=FALSE)

lines(cif2.mm2$times/365.25, cif2.mm2$cif[2,], lwd=2, col=2)
}

\keyword{Competing events}
