\name{plot.rocrisca}
\alias{plot.rocrisca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Method for 'rocrisca' Objects
}
\description{
A plot of ROC curves is produced. In the RISCA package, it concerns the functions \code{roc.binary}, \code{roc.net}, \code{roc.summary}, and \code{roc.time}.
}
\usage{
\method{plot}{rocrisca}(x, ..., information=TRUE)
}
\arguments{
  \item{x}{An object of class \code{rocrisca}, returned by the functions \code{roc.binary}, \code{roc.net}, \code{roc.summary}, and \code{roc.time}.}
  \item{...}{Additional arguments affecting the plot.}
  \item{information}{A logical value indicating whether the non-information line is plotted. The default values is TRUE.}
}
\author{
Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>
}

\examples{

# import and attach the data example

data(dataDIVAT3)

# A subgroup analysis to reduce the time needed for this exemple

dataDIVAT3 <- dataDIVAT3[1:400,]

# The ROC curve to evaluate the crude capacities of the recipient age for the
# prognosis of post kidney transplant mortality (we ignore the censoring process)

roc1 <- roc.binary(status="death", variable="ageR", confounders=~1,
data=dataDIVAT3, precision=seq(0.1,0.9, by=0.1) )

# The corresponding ROC graph with color and symbols

plot(roc1, type="b", xlab="1-specificity", ylab="sensibility", col=2, pch=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Plotting}
