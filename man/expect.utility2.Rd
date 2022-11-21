\name{expect.utility2}
\alias{expect.utility2}
\title{
Cut-Off Estimation Of A Prognostic Marker (Two Groups Are observed).
}
\description{
This function allows the estimation of a cut-off of a prognostic marker for medical decision making between two treatments A and B by maximizing the expected utility in a time-dependent context.
}
\usage{
expect.utility2(times, failures, variable, treatment, pro.time,
 u.A0, u.A1, u.B0, u.B1,  n.boot)
}
\arguments{
  \item{times}{A numeric vector with the follow up times.}
  \item{failures}{A numeric vector with the event indicator (0=right censoring, 1=event).}
  \item{variable}{A numeric vector with the observed values of the marker/variable under interest \code{X}. This variable is collected at the baseline (\code{times=0}). By convention, we assume that patients with \code{X>k} will preferentially receive A, \code{k} being the optimal cut-off. In contrast, patients with \code{X<k} will preferentially receive B.}
  \item{treatment}{A character vector with the observed treatment. Only character strings "A" and "B" are allowed.}
  \item{pro.time}{The prognostic time for which the capacities of the marker and the patient outcomes are considered in the same unit than the one used in the argument \code{times}.}
  \item{u.A0}{A value of the utility of a patient receiving the treatment A before the event occurrence. This value should respect the 0-1 scale (from death to perfect health).}
  \item{u.A1}{A value of the utility of a patient receiving the treatment A after the event occurrence. This value should respect the 0-1 scale.}
  \item{u.B0}{A value of the utility of a patient receiving the treatment B before the event occurrence. This value should respect the 0-1 scale.}
  \item{u.B1}{A value of the utility of a patient receiving the treatment B after the event occurrence. This value should respect the 0-1 scale.}
  \item{n.boot}{Number of bootstrap iterations to compute the 95\% confidence interval of the optimal cut-off. The default value is NULL: no confidence interval is estimated.}
}
\details{
The treatments A and B are allocated independently to the value of the marker \code{X}. We assume that the treatment A will be more convenient for patients with high value of the marker \code{X}. The aim is to compute the optimal cut-off value for a future stratified medical decision rule: treatment A for patients with \code{X>k} and treatment B for patients with \code{X<k}.
}
\value{
\item{estimation}{This is a single value if \code{n.boot=NULL}, which corresponds to the estimated cut-off that maximizes the time-dependent expected utility of the medical decision. If this value corresponds to the minimum of the marker, all the patients should be treated with A. If this value corresponds to the maximum of the marker, all the patients should conserve the treatment B. When \code{n.boot} is not null, two additional values are returned: \code{CIinf} is the lower bound of the 95\% confidence interval and  \code{CIsup} is the upper bound of the 95\% confidence interval.}
\item{max.eu}{This value corresponds to the maximum expected utility associated with the \code{estimation}.}
\item{table}{This data frame is composed by 8 columns representing respectively the cut-off values (\code{cut.off}), the time-dependent expected utilities (\code{utility}), the proportions of patients with a marker value higher (\code{pA}) and lower (\code{pB}) than the cut-off value, the numbers of QALYs for patients with a marker value higher (\code{qA}) and lower (\code{qB}) than the cut-off value, the RMST for patients with a marker value higher (\code{eA}) and lower (\code{eB}) than the cut-off value.}
\item{delta.rmst}{This is a vector with two values. The first value, entitled \code{high.level}, represents the  RMST for patients with a marker higher than the \code{estimation} and treated with A minus the RMST for patients with a marker higher than the \code{estimation} and treated with B. The second value, entitled \code{low.level}, represents the RMST for patients with a marker lower than or equal to the \code{estimation} and treated by B minus the RMST for patients with a marker lower than or equal to the \code{estimation} and treated with A.}
\item{delta.qaly}{This is a vector with two values. The first value, entitled \code{high.level}, represents the number of QALYs for patients with a marker higher than the \code{estimation} and treated by A minus the number of QALYs for patients with a marker higher than the \code{estimation} and treated with B. The second value, entitled \code{low.level}, represents the number of QALYs for patients with a marker lower than or equal to the \code{estimation} and treated with B minus the number of QALYs for patients with a marker lower than or equal to the \code{estimation} and treated with A.}
\item{missing}{Number of deleted observations due to missing data.}
}

\references{
Dantan et al. Optimal threshold estimator of a prognostic marker by maximizing a time-dependent expected utility function for a patient-centered stratified medicine. Statistical Methods in Medical Research, 27(6) :1847-1859. 2016. <doi:10.1177/  0962280216671161>
}

\author{
Etienne Dantan <Etienne.Dantan@univ-nantes.fr>

Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>
}

\examples{

# import and attach the data example

data(dataCSL)
dataCSL <- dataCSL[order(dataCSL$id, dataCSL$time),]
dataCSL$ordre <-0
for (i in unique(dataCSL$id))
 {dataCSL$ordre[dataCSL$id==i] <- 1:sum(dataCSL$id==i)}

dataCSL$ttt[dataCSL$treat==0]<-"A"
dataCSL$ttt[dataCSL$treat==1]<-"B"

dataCSL0 <- dataCSL[dataCSL$ordre==1,]
dataCSL0<-dataCSL0[,c(1,4,5,14,9)]

# the expected utility function for a prognostic up to 8 years

EUt.obj <- expect.utility2(dataCSL0$eventT, dataCSL0$dc,
 dataCSL0$prot.base, treatment= dataCSL0$ttt,
 pro.time=8, u.A0=0.75*0.95, u.A1=0, u.B0=0.75, u.B1=0)

plot(EUt.obj$table$cut.off, EUt.obj$table$utility, type="l",
 xlab="Cut-off values", ylab="Expected utility",col=1, lty=1)

segments(EUt.obj$estimation, 0, EUt.obj$estimation, EUt.obj$max.eu, lty=3)
segments(0, EUt.obj$max.eu, EUt.obj$estimation, EUt.obj$max.eu, lty=3)

text(EUt.obj$estimation-2, 3.38,
 paste("Optimal cut-off=",round(EUt.obj$estimation,2)), srt=90, cex=0.8)
 text(min(dataCSL0$prot.base)+40, EUt.obj$max.eu-0.005,
 paste("Expected utility=",round(EUt.obj$max.eu,2)), cex=0.8) 

# the optimal cut-off: patients with an higher value should receive the treatment A

EUt.obj$estimation
}

\keyword{Utility}
