\name{gc.survival}
\alias{gc.survival}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marginal Effect for Censored Outcome by G-computation.
}
\description{
This function allows to estimate the marginal effect of an exposure or a treatment by G-computation for a censored times-to-event.
}
\usage{
gc.survival(coxph.obj, data, group, times, failures, max.time, effect, iterations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coxph.obj}{A coxph object obtained by using the function \code{coxph} with the argument \code{x=TRUE} must be specified. It shall be included the exposition/treatment of interest variable and at least one covariate.}
  \item{data}{A data frame in which to look for the variables related to the status of the event (observed or censored), the follow-up time, the treatment/exposure and the covariables included in the previous model \code{coxph.obj}. The covariates' names have to be identical than the ones included in \code{coxph.obj}.}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones.}
  \item{times}{The name of the variable related the numeric vector with the follow-up times.}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event).}
  \item{max.time}{The maximum time of follow-up to estimate of restricted mean survival time (RMST).}
 \item{effect}{The type of marginal effect to be estimated. Three types are possible (see details): "ATE" (by default), "ATT" and "ATU".}
 \item{iterations}{The number of bootstrap resamples to estimate of the variances and the confidence intervals.}
 }
 
\details{
The ATE corresponds to Average Treatment effect on the Entire population, i.e. the marginal effect if all the sample is treated versus all the sample is untreated. The ATT corresponds to Average Treatment effect on the Treated, i.e. the marginal effect if the treated patients (\code{group} = 1) would have been untreated. The ATU corresponds to Average Treatment effect on the Untreated , i.e. the marginal effect if the untreated patients (\code{group} = 0) would have been treated. The RMST is the mean survival time of all subjects in the study population followed up to \code{max.time}, and is simply the area under the survival curve up to \code{max.time} obtained by the trapezoidal rule.
}

\value{
\item{effect}{A character string with the type of selected effect.}
\item{max.time}{A scalar related to the maximum time of follow-up.}
\item{RMST0}{A table related to the RMST in the unexposed/untreated sample: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{RMST1}{A table related to the RMST in the exposed/treated sample: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{delta}{A table related to the difference between the RMST in the exposed/treated sample minus in the unexposed/untreated one: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval, and \code{p.value} is the p-value of the bilateral test of the null hypothesis \code{RMST0} = \code{RMST1}.}
\item{logHR}{A table related to the logarithm of the average Hazard Ratio (HR): \code{estimate} is the estimated value, \code{std.error} is the corresponding standard error, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval, and \code{p.value} is the p-value of the bilateral test of the null hypothesis HR = 1.}
}

\references{
Chatton A, Le Borgne F, Leyrat C, Foucher Y. G-computation versus inverse probability weighting for time-to-event analyses. Manuscript in progress. 2019. 
}

\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>

Arthur Chatton <Arthur.Chatton@etu.univ-nantes.fr>
}

\examples{
data(dataDIVAT2)

#Raw effect of the treatment
cox.raw <- coxph(Surv(times,failures) ~ ecd, data=dataDIVAT2, x=TRUE)
summary(cox.raw)

#Conditional effect of the treatment 
cox.cdt <- coxph(Surv(times,failures) ~ ecd + age + retransplant, data=dataDIVAT2, x=TRUE)
summary(cox.cdt)

#Marginal effect of the treatment (ATE): use 1000 iterations instead of 10
#We restricted to 10 to respect the CRAN policy in terms of time for computation 
gc.ate <- gc.survival(coxph.obj=cox.cdt, data=dataDIVAT2, group="ecd", times="times",
 failures="failures", max.time=max(dataDIVAT2$times), iterations=10, effect="ATE")
gc.ate 

#Sum-up of the 3 HRs
data.frame( raw=exp(cox.raw$coefficients),
conditional=exp(cox.cdt$coefficients[1]),
marginal.ate=exp(gc.ate$logHR[,1]) )

}

\keyword{G-computation}
