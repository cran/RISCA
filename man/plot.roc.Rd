\name{plot.roc}
\alias{plot.roc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Method for 'roc' Objects
}
\description{
A plot of ROC curves is produced. In the RISCA package, it concerns the functions \code{roc.binary}, \code{roc.net}, \code{roc.summary}, and \code{roc.time}.
}
\usage{
\method{plot}{roc}(x, ..., information=TRUE)
}
\arguments{
  \item{x}{An object of class \code{roc}, returned by the functions \code{roc.binary}, \code{roc.net}, \code{roc.summary}, and \code{roc.time}.}
  \item{...}{Additional arguments affecting the plot.}
  \item{information}{A logical value indicating whether the non-information line is plotted. The default values is TRUE.}
}
\references{
Le Borgne et al. Standardized and weighted time-dependent ROC curves to evaluate the intrinsic prognostic capacities of a marker by taking into account confounding factors. Manuscript submitted. Stat Methods Med Res. 27(11):3397-3410, 2018. <doi: 10.1177/ 0962280217702416.>
}
\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>
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

plot(roc1, type="b", col=2, pch=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ROC}
