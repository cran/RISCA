\name{rmst}
\alias{rmst}

\title{Restricted Mean Survival Times.}

\description{
This function allows to estimate the Restricted Mean Survival Times (RMST).}

\usage{rmst(times, surv.rates, max.time, type) }

\arguments{
  \item{times}{A numeric vector with the times.}
  \item{surv.rates}{A numeric vector with the survival rates.}
  \item{max.time}{The maximum follow-up time.}
  \item{type}{The type of input survival cirves. Two possible arguments: "s" for a step function or "l" for a continous function.} }
 
\details{
RMST represents an interesting alternative to the hazard ratio in order to estimate the effect of an exposure. The RMST is the mean survival time in the population followed up to \code{max.time}. It corresponds to the area under the survival curve up to \code{max.time}.
}

\references{
Royston R, Parmar M. Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. BMC Medical Research Methodology 2013;13:152. <DOI: 10.1186/ 1471-2288-13-152>
}

\examples{
data(dataDIVAT2)

#Survival according to the donor status (ECD versus SCD)
res <- summary(survfit(Surv(times,failures) ~ ecd, data=dataDIVAT2))

#The mean survival time in ECD recipients followed-up to 10 years
rmst(times = res$time[as.character(res$strata)=="ecd=1"],
 surv.rates = res$surv[as.character(res$strata)=="ecd=1"],
 max.time = 10, type = "s")
 
#The mean survival time in SCD recipients followed-up to 10 years
rmst(times=res$time[as.character(res$strata)=="ecd=0"],
 surv.rates=res$surv[as.character(res$strata)=="ecd=0"],
 max.time=10, type = "s")
}