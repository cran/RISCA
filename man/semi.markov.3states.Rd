\name{semi.markov.3states}
\alias{semi.markov.3states}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3-State Semi-Markov Model
}
\description{
The 3-state SM model includes an initial state (X=1), a transient state (X=2) and an absorbing state (X=3). Usually, X=1 corresponds to disease-free or remission, X=2 to relapse, and X=3 to death. In this illness-death model, the possible transitions are: 1->2, 1->3 and 2->3.
}
\usage{
semi.markov.3states(times1, times2, sequences, weights=NULL,
dist, cuts.12=NULL, cuts.13=NULL, cuts.23=NULL, 
ini.dist.12=NULL, ini.dist.13=NULL, ini.dist.23=NULL,
cov.12=NULL, init.cov.12=NULL, names.12=NULL,
cov.13=NULL, init.cov.13=NULL, names.13=NULL,
cov.23=NULL, init.cov.23=NULL, names.23=NULL, 
conf.int=TRUE, silent=TRUE, precision=10^(-6))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times1}{A numeric vector with the observed times in days from baseline to the first transition (X=2 or X=3) or to the right-censoring (in X=1 at the last follow-up).}
  \item{times2}{A numeric vector with the observed times in days from baseline to the second transition or to the right censoring (in X=2 at the last follow-up). \code{NA} for individuals right-censored in X=1 or individuals who directly transited from X=1 to X=3.}
  \item{sequences}{A numeric vector with the sequences of observed states. Four possible values are allowed: 1 (individual right-censored in X=1), 12 (individual right-censored in X=2), 13 (individual who directly transited from X=1 to X=3), 123 (individual who transited from X=1 to X=3 through X=2).}
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
\item{object}{The character string indicating the estimated model: "semi.markov.3states (3-state semi-markov model)".}
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

Gillaizeau et al. A multistate additive relative survival semi-Markov model. Statistical methods in medical research. 26(4):1700-1711, 2017. <doi: 10.1177/ 0962280215586456>.
}

\author{
Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>

Florence Gillaizeau <Florence.Gillaizeau@univ-nantes.fr>
}

\examples{
# import the observed data
# (X=1 corresponds to initial state with a functioning graft, X=2 to acute rejection episode, 
# X=3 to return to dialysis, X=4 to death with a functioning graft)

data(dataDIVAT1)

# A subgroup analysis to reduce the time needed for this example

dataDIVAT1$id<-c(1:nrow(dataDIVAT1))
set.seed(2)
d3<-dataDIVAT1[dataDIVAT1$id \%in\% sample(dataDIVAT1$id, 250, replace = FALSE),]

# To illustrate the use of a 3-state model, individuals with trajectory 13 and 123 are 
# censored at the time of transition into state X=3

d3$trajectory[d3$trajectory==13]<-1
d3$trajectory[d3$trajectory==123]<-12
d3$trajectory[d3$trajectory==14]<-13
d3$trajectory[d3$trajectory==124]<-123

# 3-state parametric semi-Markov model including one explicative variable 
# on the transition 1->2 (z is 1 if delayed graft function and 0 otherwise).
# We only reduced the precision and the number of iteration to save time in this example,
# prefer the default values.

semi.markov.3states(times1=d3$time1, times2=d3$time2, sequences=d3$trajectory,
 dist=c("E","E","E"), ini.dist.12=c(9.93), ini.dist.13=c(11.54), ini.dist.23=c(10.21),
  cov.12=d3$z, init.cov.12=c(-0.13), names.12=c("beta12_z"),
  conf.int=TRUE, silent=FALSE, precision=0.001)$table
}

\keyword{Semi-Markov}

