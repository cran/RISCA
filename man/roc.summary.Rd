\name{roc.summary}
\alias{roc.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary ROC Curve For Aggregated Data.
}
\description{
This function computes summary ROC curve (Combescure et al., 2016).
}
\usage{
roc.summary(study.num, classe, n, year, surv, nrisk, proba, marker.min,
 marker.max, init.nlme1, precision, pro.time, time.cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{study.num}{A numeric vector (1,2,3,...) with the study identification.}
  \item{classe}{A numeric vector with integers (1,2,3,...) for identifying the groups  defined using the studied marker. 1 is the first group with the lowest values of the marker.}
  \item{n}{A numeric vector with the number of subjects at the baseline (date of marker collection).}
  \item{year}{A numeric vector with the survival times.}
  \item{surv}{A numeric vector with the survival probabilities corresponding to the previous times (often obtained graphically using the published survival curves).}
  \item{nrisk}{A numeric vector with the number of subjects at-risk of the event at the corresponding \code{year}.}
  \item{proba}{This numeric vector represents the proportion of the patients in a center which belong to the corresponding group.}
  \item{marker.min}{A numeric vector with the minimum values of the marker interval corresponding to the previous class.}
  \item{marker.max}{A numeric vector with the maximum values of the marker interval corresponding to the previous class.}
  \item{init.nlme1}{A numeric vector with the initiate values (mean, sd) of the maker distribution which is assumed to be Gaussian. Default is (0,1).}
  \item{precision}{A numeric vector with the initiate values (mean, sd) of the maker distribution which is assumed to be Gaussian. Default is \eqn{10^{-6}}.}
  \item{pro.time}{The value of prognostic time is the maximum delay for which the capacity of the variable is evaluated. The same unit than the one used in the argument \code{time}.}
  \item{time.cutoff}{The value of internal threasholds for the definition of the piecewise hazard function (3 values for a 4-piece constant function and 4 values for a 5-piece constant function).}
  }
\details{
This function computes summary ROC curve. The hazard function associated with the time-to-event was defined as a 4-piece or a 5-piece constant function with a specific association with the marker at each interval. The maker distribution is assumed to be Gaussian distributed.
}
\value{
\item{nlme1}{An object of class \code{nlme} representing the nonlinear mixed-effects model of the marker distribution. The marker is assumed Gaussian distributed. \code{mu} and \code{sigma} represent the mean and the standard deviation. The inter-study variability is modeled with a random effect on the mean. See nlmeObject for the components of the fit.}
\item{nlme2}{An object of class \code{nlme} representing the nonlinear mixed-effects model of the time distribution. The hazard function is a stepwise function with 5 intervals. \code{exp(beta0.1)} and \code{exp(beta0.2)} represent the baseline hazard and the hazard ratio in the first interval. \code{exp(beta1.1)} and \code{exp(beta1.2)} represent the corrections of these parameters for the second interval... The inter-study variability is modeled with a random effect on the baseline parameter \code{beta0.1}. See nlmeObject for the components of the fit.}
\item{table}{This data frame presents the sensitivities (\code{se}) and specificities (\code{sp}) associated with the cut-off values (\code{cut.off}).}
\item{auc}{The area under the SROC curve for a prognostic up to prognostic time.}
}
\references{
Combescure et al. A literature-based approach to evaluate the predictive capacity of a marker using time-dependent Summary Receiver Operating Characteristics. Stat Methods Med Res, 25(2):674-85, 2016. <doi: 10.1177/ 0962280212464542>
}

\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>

Christophe Combescure <christophe.combescure@hcuge.ch>
}

\examples{

# The example is too long to compute for a submission on the CRAN
# Remove the characters '#'

### import and attach the data example
data(dataKi67)

### Compute the SROC curve for a prognostic up to 9 years
roc9y<-roc.summary(dataKi67$study.num, dataKi67$classe, dataKi67$n,
 dataKi67$year, dataKi67$surv, dataKi67$nrisk, dataKi67$proba,
 dataKi67$log.marker.min, dataKi67$log.marker.max,
 init.nlme1=c(2.55, -0.29), precision=50, pro.time=9,
 time.cutoff=c(2, 4, 8))

### The ROC graph associated to these to SROC curves
plot(roc9y, col=1, lty=1, lwd=2, type="l", xlab="1-specificity", ylab="sensibility")

### Check of the goodness-of-fit: the observed proportions of
### patients in the $g$th interval of the study $k$ versus the
### fitted proportions (equation 3).

plot(roc9y$data.marker$proba, roc9y$data.marker$fitted,
 xlab="Observed probabilities", ylab="Fitted probabilities",
 ylim=c(0,1), xlim=c(0,1))
abline(0,1)

### Check of the goodness-of-fit: the observed bivariate
### probabilities versus the fitted bivariate
### probabilities (equation 4).

plot(roc9y$data.surv$p.joint, roc9y$data.surv$fitted,
 xlab="Observed probabilities", ylab="Fitted probabilities",
 ylim=c(0,1), xlim=c(0,1))
abline(0,1)

### Check of the goodness-of-fit: the residuals of the bivariate
### probabilities (equation 4) versus the times.

plot(roc9y$data.surv$year, roc9y$data.surv$resid,
 xlab="Survival time (years)", ylab="Residuals")
lines(lowess(roc9y$data.surv$year,
 I(roc9y$data.surv$resid), iter=0))
}

\keyword{ROC}

