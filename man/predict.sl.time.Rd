\name{predict.sl.time}
\alias{predict.sl.time}

\title{
Prediction from an Super Learner (SL) for Censored Outcomes
}
\description{
Predict the survival of new observations based on an SL by unsing the \code{sl.time} function.
}

\usage{
\method{predict}{sl.time}(object, ..., newdata, newtimes)
}

\arguments{
  \item{object}{An object returned by the function \code{sl.time}.}
  \item{...}{Further arguments passed.}
  \item{newdata}{An optional data frame containing covariate values at which to produce predicted values. There must be a column for every covariate included in \code{cov.quanti} and \code{cov.quali} included in the training sample. The default value is\code{NULL}, the predicted values are computed for the subjects of the training sample.}
  \item{newtimes}{The times at which to produce predicted values. The default value is\code{NULL}, the predicted values are computed for the observed times in the training data frame.}
}

\value{
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed times (columns).}
}


\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>

Camille Sabathe <camille.sabathe@univ-nantes.fr>
}

\examples{

data(dataDIVAT2)

dataDIVAT2$train <- 1*rbinom(n=dim(dataDIVAT2)[1], size = 1, prob=1/3)

# The training of the super learner with 2 algorithms
sl<-sl.time(method=c("cox.ridge", "aft.ggamma"),  metric="ribs",
  data=dataDIVAT2[dataDIVAT2$train==1,],  times="times", failures="failures", pro.time = 12,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"), cv=3)

# Individual prediction for 2 new subjects
pred <- predict(sl,
  newdata=data.frame(age=c(52,52), hla=c(0,1), retransplant=c(1,1), ecd=c(0,1)))

plot(y=pred$predictions$sl[1,], x=pred$times, xlab="Time (years)", ylab="Predicted survival",
     col=1, type="l", lty=1, lwd=2, ylim=c(0,1))

lines(y=pred$predictions$sl[2,], x=pred$times, col=2, type="l", lty=1, lwd=2)

legend("bottomright", col=c(1,2), lty=1, lwd=2, c("Subject #1", "Subject #2"))

}

\keyword{Predicting}
\keyword{Super Learner}