\name{print.sl.time}
\alias{print.sl.time}

\title{
S3 method for Printing an 'sl.time' object
}
\description{
Print the contribution of learners included in the \code{sl.time} function.
}

\usage{
\method{print}{sl.time}(x,  ..., digits)
}

\arguments{
  \item{x}{An object returned by the function \code{sl.time}.}
  \item{...}{Further arguments passed.}
  \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
}

\author{
Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>

Camille Sabathe <camille.sabathe@univ-nantes.fr>
}

\examples{

data(dataDIVAT2)

sl1<-sl.time(method=c("cox.ridge", "aft.ggamma"),  metric="ribs",
  data=dataDIVAT2[1:150,],  times="times", failures="failures", pro.time = 12,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"), cv=3)

print(sl1, digits=4)
}

\keyword{Printing}
\keyword{Super Learner}