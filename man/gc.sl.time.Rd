\name{gc.sl.time}
\alias{gc.sl.time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marginal Effect for Censored Outcome by G-computation with a Super Learner for the Outcome Model.
}
\description{
This function allows to estimate the (population-average) marginal effect of an exposure or a treatment by G-computation for a censored times-to-event, the Q-model being specified by a Super Learner.
}
\usage{
gc.sl.time(object, max.time, effect, estim.tune, estim.weights,
  methods, conditional,  iterations, length.out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{An object returned by the function \code{sl.time}. It shall be included the exposition/treatment of interest in the argument \code{group} of the \code{sl.time} function.}
\item{max.time}{The maximum time of follow-up to estimate of restricted mean survival time (RMST). The default value is \code{NULL}, which corresponds to \code{object$pro.time} if this argument is not NULL, or the time at which half of the sample is still at-risk otherwise.}
\item{effect}{The type of marginal effect to be estimated. Three types are possible (see details): "ATE" (by default), "ATT" and "ATU".}
\item{estim.tune}{A logical value. If \code{TRUE}, the tuning parameters of the algorithms of the Super Learner are estimated at each boostraap iteration.  The default is \code{TRUE}.}
\item{estim.weights}{A logical value. If \code{TRUE}, the weights of the algorithms of the Super Learner are estimated at each bootstrap iteration. The default is \code{TRUE}.}
\item{methods}{A vector of character with the methods for which the marginal effect are computing. The methods included in \code{object$methods} or "sl" are allowed. You can use "all" if you want to estimate the marginal effect of the Super Learner and all the included algorithmsor models. You can use "sl" if you want only the estimation for the Super Learner. The default value is "sl".}
\item{conditional}{A logical value. If \code{TRUE}, the conditional (subject-specific) Hazard Ratio (HR) is also computing. See the help related to the \code{hr.sl.time} function for more details on this option. The default is \code{FALSE}.}
\item{iterations}{The number of bootstrap resamples to estimate of the variances and the confidence intervals.}
\item{length.out}{An non-negative number argument related to the desired length of time intervals to compute the marginal estimands. The default is 100.}

}
 
\details{
The ATE corresponds to Average Treatment effect on the Entire population, i.e. the marginal effect if all the sample is treated versus all the sample is untreated. The ATT corresponds to Average Treatment effect on the Treated, i.e. the marginal effect if the treated patients (\code{group} = 1) would have been untreated. The ATU corresponds to Average Treatment effect on the Untreated , i.e. the marginal effect if the untreated patients (\code{group} = 0) would have been treated. The RMST is the mean survival time of all subjects in the study population followed up to \code{max.time}.
}

\value{
\item{effect}{A character string with the type of selected effect.}
\item{max.time}{A scalar related to the maximum time of follow-up for the RMST computation.}
\item{RMST0}{A table related to the RMST in the unexposed/untreated sample: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{RMST1}{A table related to the RMST in the exposed/treated sample: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{delta}{A table related to the difference between the RMST in the exposed/treated sample minus in the unexposed/untreated one: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval, and \code{p.value} is the p-value of the bilateral test of the null hypothesis \code{RMST0} = \code{RMST1}.}
\item{logHR}{A table related to the logarithm of the average HR: \code{estimate} is the estimated value, \code{std.error} is the corresponding standard error, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval, and \code{p.value} is the p-value of the bilateral test of the null hypothesis HR = 1.}
\item{table.surv}{This data frame presents the survival probabilities (\code{survival}) in each group (\code{variable}) according to the \code{times}. The number of individuals at risk (\code{n.risk}) is also provided.}
\item{logHR.conditional.values}{A vector with the estimated individual HRs (only if \code{conditional} = \code{TRUE}).}
\item{logHR.conditional}{A table related to the logarithm of the conditional HR  (only if \code{conditional} = \code{TRUE}): \code{estimate} is the estimated value, \code{std.error} is the corresponding standard error, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval, and \code{p.value} is the p-value of the bilateral test of the null hypothesis HR = 1.}
\item{learner}{A list with the same arguments (\code{RMST0}, \code{RMST1}, etc.) for the other algorithms/model included in the Super Learner ans declared in the argument \code{methods} (only if \code{methods=NULL}).}
}

\references{
Sabathe and Foucher. Super Learner for survival prediction from censored data: Extension of the R package RISCA. Manuscript submitted. 2022. 
}

\author{
Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>
}

\examples{

data(dataDIVAT2)

#The outcome model base on a Super Learner and the first 150 individuals of the data base
sl1<-sl.time( methods=c("aft.gamma", "ph.gompertz"),  metric="ibs",
  data=dataDIVAT2[1:150,],  times="times", failures="failures", group="ecd",
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant"), cv=3)
  
#Marginal effect of the treatment (ATE): use 1000 iterations instead of 2
gc.ate <- gc.sl.time(sl1, max.time=12, effect="ATE", iterations=2, 
  estim.tune=FALSE, estim.weights=FALSE)

#Plot the survival curves
plot(gc.ate, ylab="Confounder-adjusted survival",
 xlab="Time post-transplantation (years)", col=c(1,2))
}

\keyword{G-computation}
\keyword{Super Learner}