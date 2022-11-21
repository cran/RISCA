\name{lrs.multistate}
\alias{lrs.multistate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood Ratio Statistic to Compare Embedded Multistate Models
}
\description{
This function computes a Likelihood Ratio Statistic to compare two embedded multistate models.
}
\usage{
lrs.multistate(model1, model0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model1}{A list containing the results after using a function included in the present \code{RISCA} package.}
  \item{model0}{A list containing the results after using a function included in the present \code{RISCA} package. The function used to obtained the \code{model0} have to be the same than the one used to obtain the \code{model1}. The \code{model0} have to be embedded in the \code{model1}.}
 }

\value{
\item{statistic}{The value of the statistic.}
\item{ddl}{The degrees of freedom.}
\item{pvalue}{The p-value.}
}

\author{
Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>

Florence Gillaizeau <Florence.Gillaizeau@univ-nantes.fr>
}

\examples{
data(dataDIVAT1)

# A subgroup analysis to reduce the time needed for this example

dataDIVAT1$id<-c(1:nrow(dataDIVAT1))
set.seed(2)
d3<-dataDIVAT1[dataDIVAT1$id \%in\% sample(dataDIVAT1$id, 250, replace = FALSE),]

# To illustrate the use of a 3-state model, individuals with trajectory 13 and 123 are 
# censored at the time of transition into state X=3

d3$trajectory[d3$trajectory==13]<-1
d3$trajectory[d3$trajectory==123]<-12
d3$trajectory[d3$trajectory==14]<-13
d3$trajectory[d3$trajectory==124]<-123

# 3-state parametric semi-Markov model : does 'z' influence both the
# transition 1->3 ? We only reduced the precision and the number of iteration
# to save time in this example, prefere the default values.

m1 <- semi.markov.3states(times1=d3$time1, times2=d3$time2,
  sequences=d3$trajectory, dist=c("E","E","E"),
  ini.dist.12=c(9.93), ini.dist.13=c(11.54), ini.dist.23=c(10.21),
  cov.12=d3$z, init.cov.12=c(-0.13), names.12=c("beta12_z"),
  cov.13=d3$z, init.cov.13=c(1.61),  names.13=c("beta13_z"),
  conf.int=TRUE, silent=FALSE, precision=0.001)
  
m1

m0 <- semi.markov.3states(times1=d3$time1, times2=d3$time2,
  sequences=d3$trajectory, dist=c("E","E","E"),
  ini.dist.12=c(9.93), ini.dist.13=c(11.54), ini.dist.23=c(10.21),
  cov.12=d3$z, init.cov.12=c(-0.13), names.12=c("beta12_z"),
  conf.int=TRUE, silent=FALSE, precision=0.001)

m0

lrs.multistate(model1=m1, model0=m0)
}

\keyword{Markov}
\keyword{Semi-Markov}