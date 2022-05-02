\name{survival.summary}
\alias{survival.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary Survival Curve From Aggregated Data
}
\description{
Estimation of the summary survival curve from the survival rates and the numbers of at-risk individuals extracted from studies of a meta-analysis.
}
\usage{
survival.summary(study, time, n.risk, surv.rate, confidence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{study}{A numeric vector with the numbering of the studies included in the meta-analysis. The numbering of a study is repeated for each survival probabilities extracted from this study.}
  \item{time}{A numeric vector with the time at which the survival probabilities are collected.}
  \item{n.risk}{A numeric vector with the number of at-risk patients in the study for each value of thr \code{time}.}
  \item{surv.rate}{A numeric vector with the survival rates collected per study for each value of \code{time}.}
  \item{confidence}{A text argument indicating the method to calculate the 95\% confidence interval of the summary survival probabilities: "Greenwood" or "MonteCarlo".}
  }
\details{
The survival probabilities have to be extracted at the same set of points in time for all studies. Missing data are not allowed. The studies included in the meta-analysis can have different length of follow-up. For a study ending after the time t, all survival probabilities until t have to be entered in data. The data are sorted by study and by time. The conditional survival probabilities are arc-sine transformed and thus pooled assuming fixed effects or random effects. A correction of 0.25 is applied to the arc-sine transformation. For random effects, the multivariate methodology of DerSimonian and Laird is applied and the between-study covariances are accounted. The summary survival probabilities are obtained by the product of the pooled conditional survival probabilities. The mean and median survival times are derived from the summary survival curve assuming a linear interpolation of the survival between the points.
}
\value{
\item{verif.data}{A data frame in which the first column (\code{study}) correspond to the number of the study and the second column (\code{check}) equals 1 if the time of collection for this study respects the other times for the other studies and 0 otherwise. Remember that the times of survival rates have to be identical between studies. The end of each study can be different. If at least one study did not respect this format, the other arguments values are non-attributed (NA).}
\item{summary.fixed}{A matrix containing the summarized survival probabilities assuming fixed effects. The first column contains the time at which the summary survivals are computed. The second column contains the estimations of the summary survival probabilities. The third and fourth columns contain the lower and the upper bound of the 95\% confidence interval, computed by either the Greenwood or the Monte Carlo approach as specified by the user.}
\item{median.fixed}{A numerical vector containing the estimated median survival time computed from the summary survival curve assuming fixed effects and the lower and upper bounds of the 95\% confidence interval computed by a Monte Carlo approach.}
\item{mean.fixed}{A numerical vector containing the estimated mean survival time computed from the summary survival curve assuming fixed effects and the lower and upper bounds of the 95\% confidence interval computed by a Monte Carlo approach.}
\item{heterogeneity}{A numerical vector containing the value of the Q statistic for the heterogeneity, the H index and the I-squared index (in percentage).}
\item{summary.random}{A matrix containing the summarized survival probabilities assuming random effects. The first column contains the time at which the summary survivals are computed. The second column contains the estimations of the summary survival probabilities. The third and fourth columns contain the lower and the upper bound of the 95\% confidence interval around the summary survival probabilities, computed by either the Greenwood or the Monte Carlo approach as specified by the user.}
\item{median.random}{A numerical vector containing the estimated median survival time computed from the summary survival curve assuming random effects and the lower and upper bounds of the 95\% confidence interval computed by a Monte Carlo approach.}
\item{mean.random}{A numerical vector containing the estimated mean survival time computed from the summary survival curve assuming random effects and the lower and upper bounds of the 95\% confidence interval computed by a Monte Carlo approach.}
}
\references{
Combescure et al. The multivariate DerSimonian and Laird's methodology applied to meta-analysis of survival curves. 10;33(15):2521-37, 2014. Statistics in Medicine. <doi:10.1002/sim.6111>
}
\author{
Y. Foucher <Yohann.Foucher@univ-nantes.fr>

D. Jackson <daniel.jackson@mrc-bsu.cam.ac.uk>

C. Combescure <Christophe.Combescure@hcuge.ch>
}

\examples{

# import and attach the data example
data(dataHepatology)
attach(dataHepatology)

# computation of the summary survivals

results<-survival.summary(study, time, n.risk, survival, confidence="Greenwood")
results

# plot the estimated summary survival curve against the extracted ones

RandomEffectSummary<- results$summary.random

plot(time, survival, type="n", col="grey", ylim=c(0,1),xlab="Time",
 ylab="Survival")
 
for (i in unique(sort(study)))
{
lines(time[study==i], survival[study==i], type="l", col="grey")
points(max(time[study==i]),
 survival[study==i & time==max(time[study==i])], pch=15)
}

lines(RandomEffectSummary[,1], RandomEffectSummary[,2], type="l",
 col="red", lwd=3)
points(RandomEffectSummary[,1], RandomEffectSummary[,3], type="l",
 col="red", lty=3, lwd=3)
points(RandomEffectSummary[,1], RandomEffectSummary[,4], type="l",
 col="red", lty=3, lwd=3)
}

\keyword{Survival}

