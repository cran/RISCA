\name{expect.utility1}
\alias{expect.utility1}
\title{
Cut-Off Estimation Of A Prognostic Marker (Only One Observed Group).
}
\description{
This function allows the estimation of a cut-off for medical decision making between two treatments A and B from a prognostic marker by maximizing the expected utility in a time-dependent context. Only the observations of one group are available.
}
\usage{
expect.utility1(times, failures, variable, pro.time, u.A0, u.A1, u.B0, u.B1,
 n.boot, rmst.change)
}
\arguments{
  \item{times}{A numeric vector with the follow up times for the patients receiving the treatment B.}
  \item{failures}{A numeric vector with the event indicator for the patients receiving the treatment B (0=right censoring, 1=event).}
  \item{variable}{A numeric vector with the observed values of the marker under interest \code{X} for the patients receiving the treatment B. This variable is collected at the baseline (\code{times=0}). By convention, we assume that patients with \code{X>k} will preferentially receive A, \code{k} being the optimal cut-off. In contrast, patients with \code{X<k} will receive preferentially B.}
  \item{pro.time}{The prognostic time for which the prognostic capacities of the marker and the patient outcomes are considered in the same unit than the one used in the argument \code{times}.}
  \item{u.A0}{A value of the utility of a patient receiving the treatment A before the event occurrence. This value should respect the 0-1 scale (from death to perfect health).}
  \item{u.A1}{A value of the utility of a patient receiving the treatment A after the event occurrence. This value should respect the 0-1 scale.}
  \item{u.B0}{A value of the utility of a patient receiving the treatment B before the event occurrence. This value should respect the 0-1 scale.}
  \item{u.B1}{A value of the utility of a patient receiving the treatment B after the event occurrence. This value should respect the 0-1 scale.}
  \item{n.boot}{Number of bootstrap iterations to compute the 95\% confidence interval of the optimal cut-off. The default value is NULL: no confidence interval is estimated.}
  \item{rmst.change}{A numeric vector with the expected relative change in the Restricted Mean Survival Time (RMST) by using the treatment A instead of the treatment B among patients with \code{X>k}.}
}
\details{
The user observes a cohort of patients receiving the treatment B. She(he) assumes that an alternative treatment A would be more convenient for patients with high-values of the marker \code{X}. She(he) aims to compute the optimal cut-off value for a future stratified medical decision rule: treatment A for patients with \code{X>k} and treatment B for patients with \code{X<k}. The user has to enter the observed cohort of patients with the treatment B. Additional to the assumptions related to health-state utilities, the user have to specify in \code{rmst.change}, i.e. the expected relative change in terms of RMST between the two treatments. For instance, if the observed life expectancy of a patient with treatment B over the next 8 years (value entered in \code{pro.time}) is 6.70 years, and assuming that the treatment A increases this life expectancy during the next 8 years by 1.33 years, the expected relative change in RMST is 0.20 (=1.33/6.7).
}
\value{
\item{estimation}{This is a single value if \code{n.boot=NULL}, which corresponds to the estimated cut-off that maximizes the time-dependent expected utility of the medical decision. If this value corresponds to the minimum of the marker, all the patients should be treated with A. If this value corresponds to the maximum of the marker, all the patients should conserve the treatment B. When \code{n.boot} is not null, two additional values are returned: \code{CIinf} is the lower bound of the 95\% confidence interval and  \code{CIsup} is the upper bound of the 95\% confidence interval.}
\item{max.eu}{This value corresponds to the maximum expected utility associated with the \code{estimation}.}
\item{table}{This data frame is composed by 8 columns representing respectively the cut-off values, the time-dependent expected utilities (\code{utility}), the proportions of patients with a marker value higher (\code{pA}) and lower (\code{pB}) than the cut-off value, the numbers of QALYs for patients with a marker value higher (\code{qA}) and lower (\code{qB}) than the cut-off value, the RMST for patients with a marker value higher (\code{eA}) and lower (\code{eB}) than the cut-off value.}
\item{delta.rmst}{This value represents the expected RMST for patients with a marker higher than the \code{estimation} (treated with A) minus the observed RMST for patients with a marker higher than the \code{estimation}(treated with B).}
\item{delta.qaly}{This value represents the number of QALYs for patients with a marker higher than the \code{estimation} (treated with A) minus the observed number of QALYs for patients with a marker higher than the \code{estimation} (treated with B).}
\item{missing}{Number of deleted observations due to missing data.}
}

\references{
Dantan et al. Optimal threshold estimator of a prognostic marker by maximizing a time-dependent expected utility function for a patient-centered stratified medicine. Statistical Methods in Medical Research, 27(6) :1847-1859. 2016. <doi:10.1177/  0962280216671161>
}
\author{
Etienne Dantan <Etienne.Dantan@univ-nantes.fr>

Yohann Foucher <Yohann.Foucher@univ-nantes.fr>
}

\examples{

data(dataKTFS)

# to respect the CRAN policy (run times < 5s), we reduced the database
# to the first 1500 patients. Removed the first line to use the antire sample.

dataKTFS <- dataKTFS[1:1500,]
dataKTFS$score <- round(dataKTFS$score, 1) 

# the expected utility function for a prognostic up to 8 years

EUt.obj <- expect.utility1(dataKTFS$time, dataKTFS$failure, dataKTFS$score,
 pro.time=8, u.A0=0.81*0.95, u.A1=0.53, u.B0=0.81, u.B1=0.53, rmst.change=0.2)

plot(EUt.obj$table$cut.off, EUt.obj$table$utility, type="l",
 xlab="Cut-off values", ylab="Expected utility", col=1, lty=1) 
 
segments(EUt.obj$estimation, 0, EUt.obj$estimation, EUt.obj$max.eu, lty=3)
segments(0, EUt.obj$max.eu, EUt.obj$estimation, EUt.obj$max.eu, lty=3)

text(EUt.obj$estimation-0.2, 6.22,
 paste("Optimal cut-off=", round(EUt.obj$estimation,2)), srt=90, cex=0.8)
 text(min(dataKTFS$score)+1.4, EUt.obj$max.eu-0.006,
 paste("Expected utility=", round(EUt.obj$max.eu, 2)), cex=0.8) 

# the optimal cut-off: patients with an higher value should receive the treatment A

EUt.obj$estimation
}


\keyword{Utility}
