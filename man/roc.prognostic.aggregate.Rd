\name{roc.prognostic.aggregate}
\alias{roc.prognostic.aggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prognostic ROC Curve Based on Survival Probabilities
}
\description{
Prognostic ROC curve is an alternative graphical approach to represent the discriminative capacity of the marker: a receiver operating characteristic (ROC) curve by plotting 1 minus the survival in the high-risk group against 1 minus the survival in the low-risk group. The area under the curve (AUC) corresponds to the probability that a patient in the low-risk group has a longer lifetime than a patient in the high-risk group. The prognostic ROC curve provides complementary information compared to survival curves. The AUC is assessed by using the trapezoidal rules. When survival curves do not reach 0, the prognostic ROC curve is incomplete and the extrapolations of the AUC are performed by assuming pessimist, optimist and non-informative situations. The user enters the survival according to the model she/he chooses. The area under the prognostic ROC curve is assessed by using the trapezoidal rules. The extrapolated areas (when survival curves do not reach 0) are performed by assuming pessimist, optimist and non-informative situation.
}
\usage{
roc.prognostic.aggregate(time.lr, surv.lr, time.hr, surv.hr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time.lr}{A numeric vector with the survival times in the low risk group.}
  \item{surv.lr}{A numeric vector with the survival probabilities corresponding to \code{time.lr}}
  \item{time.hr}{A numeric vector with the survival times in the high risk group.}
  \item{surv.hr}{A numeric vector with the survival probabilities corresponding to \code{time.hr}.}
  }
\details{
The maximum prognostic time is the minimum between the maximum of \code{time.lr} and the maximum of \code{time.hr}.
}
\value{
\item{max.time}{This is the maximum prognostic time used for the analysis}
\item{table}{This data frame presents the different time cut-offs associated with the coordinates of the ROC curves.}
\item{auc}{This data frame presents the different estimations of the area under the prognostic ROC curve: the lower bound, the pessimist, the non-informative, the optimist and the upper bound.}
}
\references{
Combescure C, Perneger TV, Weber DC, Daures JP and Foucher Y. Prognostic ROC curves: a method for representing the overall discriminative capacity of binary markers with right-censored time-to-event endpoints. Epidemiology 2014 Jan;25(1):103-9. <doi: 10.1097/EDE.0000000000000004>
}
\author{
Y. Foucher <Yohann.Foucher@univ-nantes.fr> and

C. Combescure <Christophe.Combescure@hcuge.ch>
}

\examples{

# example of two survival curves using exponential distributions
time.hr <- seq(0, 600, by=5)
time.lr <- seq(0, 500, by=2)
surv.hr <- exp(-0.005*time.hr)
surv.lr <- exp(-0.003*time.lr)

# Illustration of both survival curves
plot(time.hr, surv.hr, xlab="Time (in days)",
 ylab="Patient survival", lwd=2, type="l")
lines(time.lr, surv.lr, lty=2, col=2, lwd=2)
legend("topright", c("High-Risk Group", "Low-Risk Group"), lwd=2,
 col=1:2, lty=1:2)

# Computation of the prognostic ROC curve
proc.result <- roc.prognostic.aggregate(time.lr, surv.lr, time.hr, surv.hr)

# Representation of the prognostic ROC curve
plot(proc.result$table$x, proc.result$table$y, type="l",
 lwd=2, xlim=c(0,1), ylim=c(0,1),
 xlab="1-Survival in the low risk group",
 ylab="1-Survival in the high risk group")
abline(c(0,0), c(1,1), lty=2)

# The pessimist value of the area under the curve
proc.result$auc$pessimist
}

\keyword{ROC}
\keyword{Survival}
