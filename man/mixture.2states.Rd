\name{mixture.2states}
\alias{mixture.2states}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Horizontal Mixture Model for Two Competing Events
}
\description{
The 2-state mixture model which includes an initial state (X=1) and two absorbing states in competition (X=2 and X=3). Parameters are estimated by (weighted) Likelihood maximization.
}
\usage{
mixture.2states(times, sequences, weights=NULL, dist, cuts.12=NULL, cuts.13=NULL,
 ini.dist.12=NULL, ini.dist.13=NULL, cov.12=NULL, init.cov.12=NULL,
 names.12=NULL, cov.13=NULL, init.cov.13=NULL, names.13=NULL,
 cov.p=NULL, init.cov.p=NULL, names.p=NULL, init.intercept.p=NULL,
 conf.int=TRUE, silent=TRUE, precision=10^(-6))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{A numeric vector with the observed times in days from baseline to the last observation.}
  \item{sequences}{A numeric vector with the sequence of observed states. Three possible values are allowed: 1 (the individual is right-censored in X=1), 12 (the individual transits to X=2) and 13 (the individual transits to X=3).}
  \item{weights}{A numeric vector with the weights for correcting the contribution of each individual. When the vector is completed, the IPW estimator is implemented. Default is \code{NULL} which means that no weighting is applied.}
  \item{dist}{A character vector with two arguments describing respectively the distributions of duration time for transitions 1->2 and 1->3. Arguments allowed are \code{"E"} for Exponential distribution, \code{"PE"} for the piecewise exponential distribution, \code{"W"} for Weibull distribution or \code{"WG"} for Generalized Weibull distribution. When the user choose \code{"PE"}, the arguments \code{"cut.XX"} have also to be defined.}
  \item{cuts.12}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=1 to X=2. Only internal timepoints are allowed: timepoints cannot be \code{0} or \code{Inf}. Default is \code{NULL} which means that the distribution is not piecewise. Piecewise model is only allowed for exponential distribution.}
  \item{cuts.13}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=1 to X=3. Only internal timepoints are allowed: timepoints cannot be \code{0} or \code{Inf}. Default is \code{NULL} which means that the distribution is not piecewise. Piecewise model is only allowed for exponential distribution.}
  \item{ini.dist.12}{A numeric vector of initial values for the distribution from X=1 to X=2. The logarithm of the parameters have to be declared. Default value is 1.}
  \item{ini.dist.13}{A numeric vector of initial values for the distribution from X=1 to X=3. The logarithm of the parameters have to be declared. Default value is 1.}
 \item{cov.12}{A matrix (or data frame) with the explicative time-fixed variable(s) related to the time from X=1 to X=2.}
 \item{init.cov.12}{A numeric vector of initial values for regression coefficients (logarithm of the cause-specific hazards ratios) associated to \code{cov.12}. Default initial value is 0.}
 \item{names.12}{An optional character vector with name of explicative variables associated to \code{cov.12}.}
 \item{cov.13}{A numeric matrix (or data frame) with the explicative time-fixed variable(s) related to the time from X=1 to X=3.}
 \item{init.cov.13}{A numeric vector of initial values for regression coefficients (logarithm of the cause-specific hazards ratios) associated to \code{cov.13}. Default initial value is 0.}
 \item{names.13}{An optional character vector with name of explicative variables associated to \code{cov.13}.}
  \item{cov.p}{A matrix (or data frame) with the explicative time-fixed variable(s) related to the probability P(X=2), which is regressing according to a logistic function.}
 \item{init.cov.p}{A numeric vector of initial values for regression coefficients (logarithm of the cause-specific hazards ratios) associated to \code{cov.p}. Default initial value is 0.}
 \item{names.p}{An optional character vector with name of explicative variables associated to \code{cov.p}.}
 \item{init.intercept.p}{A numeric value to iniate the intercept of the logit of P(X=2). Default value is 0.}
 \item{conf.int}{A logical value specifying if the pointwise confidence intervals for parameters and the variance-covariance matrix should be returned. Default is \code{TRUE}.}
 \item{silent}{A logical value specifying if the log-likelihood value should be returned at each iteration. Default is \code{TRUE}, which corresponds to silent mode (no display).}
 \item{precision}{A numeric positive value indicating the required precision for the log-likelihood maximization between each iteration. Default is \eqn{10^{-6}}.}
 }
 
\details{
Hazard functions available are:
\tabular{rl}{
Exponential distribution \tab \eqn{\lambda(t)=1/\sigma} \cr											
Weibull distribution \tab \eqn{\lambda(t)=\nu(\frac{1}{\sigma})^{\nu}t^{\nu-1}} \cr
Generalized Weibull distribution \tab \eqn{\lambda(t)=\frac{1}{\theta}\left(1+\left(\frac{t}{\sigma}\right)^{\nu}\right)^{\frac{1}{\theta}-1} \nu\left(\frac{1}{\sigma}\right)^{\nu} t^{\nu-1}}
}
with \eqn{\sigma}, \eqn{\nu},and \eqn{\theta>0}. The parameter \eqn{\sigma} varies for each interval when the distribution is piecewise Exponential. We advise to initialize the logarithm of these parameters in \code{ini.dist.12}, \code{ini.dist.13} and \code{ini.dist.23}.

To estimate the marginal effect of a binary exposure, the \code{weights} may be equal to \code{1/p}, where \code{p} is the estimated probability that the individual belongs to his or her own observed group of exposure. The probabilities \code{p} are often estimated by a logistic regression in which the dependent binary variable is the exposure. The possible confounding factors are the explanatory variables of this logistic model.
}

\value{
\item{object}{The character string indicating the estimated model: "mixture.2states (mixture model with two competing events)".}
\item{dist}{A character vector with two arguments describing respectively the distributions of duration time for transitions 1->2 and 1->3.}
\item{cuts.12}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=1 to X=2.}
\item{cuts.13}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=1 to X=3.}
\item{covariates}{A numeric vector indicating the numbers of covariates respectively related to the time to the event X=2, the time to the event X=3, the long-term probability P(X=2).}
\item{table}{A data frame containing the estimated parameters of the model (\code{Estimate}). When the option \code{conf.int=TRUE} is specified, this data frame includes three additional columns: the Standard Errors of parameters (\code{Std.Error}), the value of the Wald statistic (\code{t.value}), and the p-value for the Wald test (\code{Pr(>|t|)}).}
\item{cov.matrix}{A data frame corresponding to variance-covariance matrix of the parameters.}
\item{LogLik}{A numeric value corresponding to the (weighted) log-likelihood of the model.}
\item{AIC}{A numeric value corresponding to the Akaike Information Criterion of the model.}
}

\references{
Trebern-Launay K, Kessler M, Bayat-Makoei S, Querard AH, Briancon S, Giral M, Foucher Y. Horizontal mixture model for competing risks: a method used in waitlisted renal transplant candidates. European Journal of Epidemiology. 33(3):275-286, 2018. <doi: 10.1007/s10654-017-0322-3>.
}
\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>
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

mm2.test <- mixture.2states(times=d2$time, sequences=d2$trajectory, weights=NULL,
  dist=c("E","W"), cuts.12=NULL, cuts.13=NULL, 
  ini.dist.12=c(9.28), ini.dist.13=c(9.92, -0.23), 
  cov.12=d2$z, init.cov.12=0.84, names.12="beta_12",
  cov.13=d2$z, init.cov.13=0.76, names.13="beta_13",
  cov.p=NULL, init.cov.p=NULL, names.p=NULL, init.intercept.p=-0.75,
  conf.int=TRUE, silent=FALSE)
  
mm2.test$table
}

\keyword{Competing events}
