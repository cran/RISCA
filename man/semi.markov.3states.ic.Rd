\name{semi.markov.3states.ic}
\alias{semi.markov.3states.ic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3-State Semi-Markov Model With Interval-Censored Data
}
\description{
The 3-state SM model includes an initial state (X=1), a transient state (X=2) and an absorbing state (X=3). Usually, X=1 corresponds to disease-free or remission, X=2 to relapse, and X=3 to death. In this illness-death model, the possible transitions are: 1->2, 1->3 and 2->3. The time from X=1 to X=2 is interval-censored. Parameters are estimated by (weighted) Likelihood maximization.
}
\usage{
semi.markov.3states.ic(times0, times1, times2, sequences, weights=NULL,
dist, cuts.12=NULL, cuts.13=NULL, cuts.23=NULL, 
ini.dist.12=NULL, ini.dist.13=NULL, ini.dist.23=NULL,
cov.12=NULL, init.cov.12=NULL, names.12=NULL,
cov.13=NULL, init.cov.13=NULL, names.13=NULL,
cov.23=NULL, init.cov.23=NULL, names.23=NULL, 
conf.int=TRUE, silent=TRUE, precision=10^(-6), 
legendre=30, homogeneous=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times0}{A numeric vector with the observed times in days from baseline to the last observation time in X=1.}
  \item{times1}{A numeric vector with the observed times in days from baseline to the first observation time in X=2. \code{NA} for individuals right-censored in X=1 or individuals who are directly in X=3 after X=1 (without any observation in X=2).}
  \item{times2}{A numeric vector with the observed times in days from baseline to the last follow-up.}
  \item{sequences}{A numeric vector with the sequences of observed states. Four possible values are allowed: 1 (individual right-censored in X=1), 12 (individual right-censored in X=2), 13 (individual who directly observed in X=3 after X=3, without any observation of X=2), 123 (individual who transited from X=1 to X=3 through X=2).}
  \item{weights}{A numeric vector with the weights for correcting the contribution of each individual. When the vector is completed, the IPW estimator is implemented. Default is \code{NULL} which means that no weighting is applied.}
  \item{dist}{A character vector with three arguments describing respectively the distributions of duration time for transitions 1->2, 1->3 and 2->3. Arguments allowed are \code{"E"} for Exponential distribution, \code{"PE"} for the piecewise exponential distribution, \code{"W"} for Weibull distribution or \code{"WG"} for Generalized Weibull distribution. When the user choose \code{"PE"}, the arguments \code{"cut.XX"} have also to be defined.}
  \item{cuts.12}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=1 to X=2. Only internal timepoints are allowed: timepoints cannot be \code{0} or \code{Inf}. Default is \code{NULL} which means that the distribution is not piecewise. Piecewise model is only allowed for exponential distribution.}
  \item{cuts.13}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=1 to X=3. Only internal timepoints are allowed: timepoints cannot be \code{0} or \code{Inf}. Default is \code{NULL} which means that the distribution is not piecewise. Piecewise model is only allowed for exponential distribution.}
  \item{cuts.23}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=2 to X=3. Only internal timepoints are allowed: timepoints cannot be \code{0} or \code{Inf}. Default is \code{NULL} which means that the distribution is not piecewise. Piecewise model is only allowed for exponential distribution.}
  \item{ini.dist.12}{A numeric vector of initial values for the distribution from X=1 to X=2. The logarithm of the parameters have to be declared. Default value is 1.}
  \item{ini.dist.13}{A numeric vector of initial values for the distribution from X=1 to X=3. The logarithm of the parameters have to be declared. Default value is 1.}
  \item{ini.dist.23}{A numeric vector of initial values for the distribution from X=2 to X=3. The logarithm of the parameters have to be declared. Default value is 1.}
 \item{cov.12}{A matrix (or data frame) with the explicative time-fixed variable(s) related to the time from X=1 to X=2.}
 \item{init.cov.12}{A numeric vector of initial values for regression coefficients (logarithm of the cause-specific hazards ratios) associated to \code{cov.12}. Default initial value is 0.}
 \item{names.12}{An optional character vector with name of explicative variables associated to \code{cov.12}.}
 \item{cov.13}{A numeric matrix (or data frame) with the explicative time-fixed variable(s) related to the time from X=1 to X=3.}
 \item{init.cov.13}{A numeric vector of initial values for regression coefficients (logarithm of the cause-specific hazards ratios) associated to \code{cov.13}. Default initial value is 0.}
 \item{names.13}{An optional character vector with name of explicative variables associated to \code{cov.13}.}
  \item{cov.23}{A numeric matrix (or data frame) with the explicative time-fixed variable(s) related to the time from X=2 to X=3.}
 \item{init.cov.23}{A numeric vector of initial values for regression coefficients (logarithm of the cause-specific hazards ratios) associated to \code{cov.23}. Default initial value is 0.}
 \item{names.23}{An optional character vector with name of explicative variables associated to \code{cov.23}.}
 \item{conf.int}{A logical value specifying if the pointwise confidence intervals for parameters and the variance-covariance matrix should be returned. Default is \code{TRUE}.}
 \item{silent}{A logical value specifying if the log-likelihood value should be returned at each iteration. Default is \code{TRUE}, which corresponds to silent mode (no display).}
 \item{precision}{A numeric positive value indicating the required precision for the log-likelihood maximization between each iteration. Default is \eqn{10^{-6}}.}
 \item{legendre}{A numeric value indicating the number of knots and weights for Gaussian quadrature used in convolution products. Default is 30.}
 \item{homogeneous}{A logical value specifying if the time spent in the state X=1 is considered as non-associated with the distribution of the time from the entry in the state X=2 to the transition in the state X=3. Default is \code{TRUE}, assuming no association.}
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

Two kinds of model can be estimated: homogeneous and non-homogeneous semi-Markov model. In the first one, the hazard functions only depend on the times spent in the corresponding state. Note that for the transitions from the state X=1, the time spent in the state corresponds to the chronological time from the baseline of the study, as for Markov models. In the second one, the hazard function of the transition from the state X=2 to X=3 depends on two time scales: the time spent in the state 2 which is the random variable of interest, and the time spend in the state X=1 as a covariate.
}

\value{
\item{object}{The character string indicating the model: "semi.markov.3states.ic (3-state semi-markov model with interval-censored data)".}
\item{dist}{A character vector with two arguments describing respectively the distributions of duration time for transitions 1->2, 1->3 and 2->3.}
\item{cuts.12}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=1 to X=2.}
\item{cuts.13}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=1 to X=3.}
\item{cuts.23}{A numeric vector indicating the timepoints in days for the piecewise exponential distribution related to the time from X=2 to X=3.}
\item{covariates}{A numeric vector indicating the numbers of covariates respectively related to the transition 1->2, 1->3 and 2->3.}
\item{table}{A data frame containing the estimated parameters of the model (\code{Estimate}). When the option \code{conf.int=TRUE} is specified, this data frame includes three additional columns: the Standard Errors of parameters (\code{Std.Error}), the value of the Wald statistic (\code{t.value}), and the related p-value for the Wald test (\code{Pr(>|t|)}).}
\item{cov.matrix}{A data frame corresponding to variance-covariance matrix of the parameters.}
\item{LogLik}{A numeric value corresponding to the (weighted) log-likelihood of the model.}
\item{AIC}{A numeric value corresponding to the Akaike Information Criterion of the model.}
}

\references{
Gillaizeau et al. Inverse Probability Weighting to control confounding in an illness-death model for interval-censored data. Stat Med. 37(8):1245-1258, 2018. <doi: 10.1002/sim.7550>.
}
\author{
Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>

Florence Gillaizeau <Florence.Gillaizeau@univ-nantes.fr>
}

\examples{
# The example is too long to compute for a submission on the CRAN
# Remove the characters '#'

# import the observed data (read the application in Gillaizeau et al. for more details)
# X=1 corresponds to initial state with a functioning graft, X=2 to acute rejection episode, 
# X=3 to return to dialysis, X=4 to death with a functioning graft

# data(dataDIVAT1)

# A subgroup analysis to reduce the time needed for this example

# dataDIVAT1$id<-c(1:nrow(dataDIVAT1))
# set.seed(2)
# d3<-dataDIVAT1[dataDIVAT1$id \%in\% sample(dataDIVAT1$id, 100, replace = FALSE),]

# To illustrate the use of a 3-state model, the return in dialysis are right-censored

# d3$trajectory[d3$trajectory==13]<-1
# d3$trajectory[d3$trajectory==123]<-12
# d3$trajectory[d3$trajectory==14]<-13
# d3$trajectory[d3$trajectory==124]<-123
# table(d3$trajectory)

# X=2 is supposed to be interval-censored between 'times0' and 'times1' because
# health examinations take place each year after inclusion

# d3$times0<-NA
# d3$times1<-NA
# d3$time2_<-NA

# i<-d3$trajectory==1
# d3$times0[i]<-trunc(d3$time1[i]/365.24)*365.24+1
# d3$times1[i]<-NA
# d3$times2[i]<- d3$time1[i]+1

# i<-d3$trajectory==12
# d3$times0[i]<-trunc(d3$time1[i]/365.24)*365.24+1
# d3$times1[i]<-(trunc(d3$time1[d3$trajectory==12]/365.24)+1)*365.24
# d3$times2[i]<-pmax(d3$time2[i], (trunc(d3$time1[i]/365.24)+2)*365.24)

# i<-d3$trajectory==13
# d3$times0[i]<-trunc(d3$time1[i]/365.24)*365.24+1
# d3$times1[i]<-NA
# d3$times2[i]<-d3$time1[i]

# i<-d3$trajectory==123
# d3$times0[i]<-trunc(d3$time1[i]/365.24)*365.24+1
# d3$times1[i]<-(trunc(d3$time1[i]/365.24)+1)*365.24 
# d3$times2[i]<- pmax(d3$time2[i], (trunc(d3$time1[i]/365.24)+2)*365.24)

# 3-state homogeneous semi-Markov model with interval-censored data 
# including one binary explicative variable (z is 1 if delayed graft function and
# 0 otherwise).
# Estimation of the marginal effect of z on the transition from X=1 to X=2 
# by adjusting for 2 possible confounding factors (age and gender)
# We only reduced the precision and the number of iteration to save time in this example,
# prefer the default values.

# propensity.score <- glm(z ~ ageR + sexR, family=binomial(link="logit"),data=d3)  
# d3$fit<-propensity.score$fitted.values
  
# p1<-mean(d3$z)
# d3$w <- p1/d3$fit
# d3$w[d3$z==0]<-(1-p1)/(1-d3$fit[d3$z==0])
 
# semi.markov.3states.ic(times0=d3$times0, times1=d3$times1,
# times2=d3$times2, sequences=d3$trajectory,
# weights=d3$w, dist=c("E","E","E"), cuts.12=NULL, cuts.13=NULL, cuts.23=NULL,
# ini.dist.12=c(8.23), ini.dist.13=c(10.92), ini.dist.23=c(10.67),
# cov.12=d3$z, init.cov.12=c(0.02), names.12=c("beta12_z"),
# conf.int=TRUE, silent=FALSE, precision=0.001, legendre=20)$table
}

\keyword{Semi-Markov}

