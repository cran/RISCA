\name{survival.summary.strata}
\alias{survival.summary.strata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary Survival Curve And Comparison Between Strata.
}
\description{
Estimation of the summary survival curve from the survival rates and the numbers of at-risk individuals extracted from studies of a meta-analysis and comparisons between strata of studies.
}
\usage{
survival.summary.strata(study, time, n.risk, surv.rate, confidence, strata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{study}{A numeric vector with the numbering of the studies included in the meta-analysis. The numbering of a study is repeated for each survival probabilities extracted from this study.}
  \item{time}{A numeric vector with the time at which the survival probabilities are collected. }
  \item{n.risk}{A numeric vector with the number of at-risk patients in the study for each value of thr \code{time}.}
  \item{surv.rate}{A numeric vector with the survival rates collected per study for each value of \code{time}.}
  \item{confidence}{A text argument indicating the method to calculate the 95\% confidence interval of the summary survival probabilities: "Greenwood" or "MonteCarlo".}
  \item{strata}{A factor designing the strata. Each stratum has to contain at least two studies.}
  }
\details{
The survival probabilities have to be extracted at the same set of points in time for all studies. Missing data are not allowed. The studies included in the meta-analysis can have different length of follow-up. For a study ending after the time t, all survival probabilities until t have to be entered in data. The data are sorted by study and by time. The conditional survival probabilities are arc-sine transformed and thus pooled assuming fixed effects or random effects. A correction of 0.25 is applied to the arc-sine transformation. For random effects, the multivariate  methodology of DerSimonian and Laird is applied and the between-study covariances are accounted. The summary survival probabilities are obtained by the product of the pooled conditional survival probabilities. The mean and median survival times are derived from the summary survival curve assuming a linear interpolation of the survival between the points. The summary survival curve is assessed in each stratum. The duration of follow-up is the greatest duration for which each stratum contains at least two studies reporting the survival at this duration. The between-strata is assessed and tested.
}
\value{
\item{verif.data}{A data frame in which the first column (\code{study}) correspond to the number of the study and the second column (\code{check}) equals 1 if the time of collection for this study respects the other times for the other studies and 0 otherwise. Remember that the times of survival rates have to be identical between studies. The end of each study can be different. If at least one study did not respect this format, the other arguments values are non-attributed (NA).}
\item{summary.fixed}{A list list of matrix. Each matrix contains the summarized survival probabilities assuming fixed effects. Each matrix provides the results for one stratum. The first column contains the time at which the summary survivals are computed. The second column contains the estimations of the summary survival probabilities. The third and fourth columns contain the lower and the upper bound of the 95\% confidence interval, computed by either the Greenwood or the Monte Carlo approach as specified by the user. The last element of the list is the summary survival when all strata are pooled.}
\item{summary.random}{A list object containing the summarized survival probabilities in each stratum assuming random effects. The results are presented similarly as \code{summary.fixed}.}
\item{median.fixed}{A numerical vector containing the estimated median survival time computed from the summary survival curve assuming fixed effects and the lower and upper bounds of the 95\% confidence interval computed by a Monte Carlo approach.}
\item{mean.fixed}{A numerical vector containing the estimated mean survival time computed from the summary survival curve assuming fixed effects and the lower and upper bounds of the 95\% confidence interval computed by a Monte Carlo approach.}
\item{heterogeneity}{A numerical vector containing the value of the Q statistic for the heterogeneity, the H index and the I-squared index (in percentage).}
\item{p.value}{The p-value of the test for the null hypothesis that the between-strata heterogeneity is null.}
}
\references{
Combescure et al. The multivariate DerSimonian and Laird's methodology applied to meta-analysis of survival curves. 10;33(15):2521-37, 2014. Statistics in Medicine. <doi:10.1002/sim.6111>.
}
\author{
Y. Foucher <Yohann.Foucher@univ-poitiers.fr>

D. Jackson <daniel.jackson@mrc-bsu.cam.ac.uk>

C. Combescure <Christophe.Combescure@hcuge.ch>
}

\examples{
# import and attach the data example
data(dataHepatology)

attach(dataHepatology)

# computation of the summary survivals

results<-survival.summary.strata(study = study, time = time, n.risk = n.risk,
 surv.rate = survival, confidence="Greenwood", strata = location)

results
}

\keyword{Survival}
