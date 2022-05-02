\name{roc.prognostic.individual}
\alias{roc.prognostic.individual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prognostic ROC Curve based on Individual Data
}
\description{
The user enters individual survival data. The area under the prognostic ROC curve is assessed by using the trapezoidal rules. The extrapolated areas (when survival curves do not reach 0) are performed by assuming pessimist, optimist and non-informative situation. The confidence intervals are obtained by non-parametric bootstrapping.
}
\usage{
roc.prognostic.individual(times, failures, variable, iterations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{A numeric vector with the follow up times.}
  \item{failures}{A numeric vector with the event indicator (0=right censored, 1=event).}
  \item{variable}{A numeric vector with the result of the binary test (only two groups). The variable equals 1 for the high risk group and 0 for the low risk group.}
  \item{iterations}{The number of bootstrap samples to compute the confidence intervals. The default value is 0, which corresponds to no confidence interval computation.}
  }
\details{
The maximum prognostic time is the minimum between the two last observed times of failure in both groups. Prognostic ROC curve is an alternative graphical approach to represent the discriminative capacity of the marker: a receiver operating characteristic (ROC) curve by plotting 1 minus the survival in the high-risk group against 1 minus the survival in the low-risk group. The area under the curve (AUC) corresponds to the probability that a patient in the low-risk group has a longer lifetime than a patient in the high-risk group. The prognostic ROC curve provides complementary information compared to survival curves. The AUC is assessed by using the trapezoidal rules. When survival curves do not reach 0, the prognostic ROC curve is incomplete and the extrapolations of the AUC are performed by assuming pessimist, optimist and non-informative situations.
}
\value{
\item{max.time}{This is the maximum prognostic time used for the analysis}
\item{table}{This data frame presents the different time cut-offs associated with the coordinates of the ROC curves.}
\item{auc}{This data frame presents the different estimations of the area under the prognostic ROC curve: the lower bound, the pessimist, the non-informative, the optimist and the upper bound.}
\item{CI.95}{This data frame presents the 95 percentage confidence intervals of the area under the prognostic ROC curve: the lower bound, the pessimist, the non-informative, the optimist and the upper bound.}
\item{auc.boot}{This data frame presents the different estimations of the area under the prognostic ROC curve for each of the \code{B} bootstrap samples: the lower bound, the pessimist, the non-informative, the optimist and the upper bound.}
}
\references{
Combescure C, Perneger TV, Weber DC, Daures JP and Foucher Y. Prognostic ROC curves: a method for representing the overall discriminative capacity of binary markers with right-censored time-to-event endpoints. Epidemiology 2014 Jan;25(1):103-9. <doi: 10.1097/EDE.0000000000000004>
}
\author{
Y. Foucher <Yohann.Foucher@univ-nantes.fr> and

C. Combescure <Christophe.Combescure@hcuge.ch>
}

\examples{

###################################################################
# example of two samples with different exponential distributions #
###################################################################

n1 <- 200
n2 <- 200
grp  <- c(rep(1, n1), rep(0, n2))
time.evt <- c(rexp(n1, rate = 1.2), rexp(n2, rate = 0.5))
time.cen <- rexp(n1+n2, rate = 0.2)
time <- pmin(time.evt, time.cen)
evt  <- 1*(time.evt < time.cen)

# Illustration of both survival curves
surv.temp <- survfit(Surv(time, evt) ~ grp) 
plot(surv.temp, lty = 2:3)

# Computation of the prognostic ROC curve
proc.result <- roc.prognostic.individual(time, evt, grp, iterations=50)
  # Use more than 50 bootstrap samples for real applications
  
# Representation of the prognostic ROC curve
plot(proc.result$table$x, proc.result$table$y, type="l",
 lwd=2, xlim=c(0,1), ylim=c(0,1),
 xlab="1-Survival in the low risk group",
 ylab="1-Survival in the high risk group")
abline(c(0,0), c(1,1), lty=2)

# The corresponding 95% CI of the pessimist value
proc.result$CI.95$pessimist

}

\keyword{ROC}
\keyword{Survival}