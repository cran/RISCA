\name{plot.sl.time}
\alias{plot.sl.time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Caliration Plot for Super Learner
}
\description{
A calibration plot of a Super Learner obtained by the function \code{sl.time}.
}

\usage{
\method{plot}{sl.time}(x, ..., col=1, lty=1, lwd=1, type="b", pch = 16, ylab=NULL, xlab=NULL, method="sl",
 n.groups=5, pro.time=NULL, newdata=NULL, times=NULL, failures=NULL)
}

\arguments{
  \item{x}{An object returned by the function \code{sl.time}.}
  \item{...}{Additional arguments affecting the plot.}
  \item{col}{A numeric value with the color of the survival curves. The default is 1 for black.}
  \item{lty}{A numeric value with the type of the survival curves. The default is 1.}
  \item{lwd}{A numeric value with the type of the survival curves. The default is 1.}
  \item{type}{A character string giving the type of plot desired. The default is "b" for both the lines and points}
  \item{pch}{A numeric value for plotting characters or symbols. The default is 16.}
  \item{xlab}{The label for the x axis. The default is "Predicted survival".}
  \item{ylab}{The label for the y axis. The default is "Observed survival".}
  \item{method}{A character string with the name of the algorithm included in the SL for wich the calibration plot is performed. The defaut is "sl" for the Super Learner.}
  \item{n.groups}{A numeric value with the number of groups by their class probabilities. The default is 5.}
  \item{pro.time}{The prognostic time at wich the calibration plot of the survival probabilities.}
    \item{newdata}{An optional data frame containing the new sample for validation with covariate values, folow-up times, and event status. The default value is\code{NULL}, the calibration plot is performed from the same subjects of the training sample.}
  \item{times}{The name of the variable related the numeric vector with the follow-up times in \code{newdata} (optional argument only necessary when newdata is not NULL).}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event) in \code{newdata} (optional argument only necessary when newdata is not NULL).}
}


\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>

Camille Sabathe <camille.sabathe@univ-nantes.fr>
}

\examples{

data(dataDIVAT2)

dataDIVAT2$train <- 1*rbinom(n=dim(dataDIVAT2)[1], size = 1, prob=1/3)

# The training of the super learner with 2 models
sl<-sl.time(method=c("aft.gamma", "ph.exponential"),
  metric="ribs",  data=dataDIVAT2[dataDIVAT2$train==1,],
  times="times", failures="failures", pro.time = 12,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"),  cv=3)

# The calibration plot from the validation sample
plot(sl, method="sl", n.groups=5, pro.time=12, col=2,
     xlab="Predicted 12-year survival", ylab="Observed 12-year survival",
     newdata=dataDIVAT2[dataDIVAT2$train==0,], times="times", failures="failures")

}

\keyword{Plotting}
