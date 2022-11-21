\name{ipw.log.rank}
\alias{ipw.log.rank}
\title{
Log-Rank Test for Adjusted Survival Curves.
}
\description{
The user enters individual survival data and the weights previously calculated (by using logistic regression for instance). The usual log-rank test is adapted to the corresponding adjusted survival curves.
}
\usage{
ipw.log.rank(times, failures, variable, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{A numeric vector with the follow up times.}
  \item{failures}{A numeric vector with the event indicator (0=right censored, 1=event).}
  \item{variable}{A numeric vector with the binary variable under interest (only two groups).}
  \item{weights}{The weights for correcting the contribution of each individual. By default, the weights are all equaled to 1 and the survival curves correspond to the usual Kaplan-Meier estimator.}
 }
\details{
For instance, the \code{weights} may be equal to 1/\code{p}, where \code{p} is the estimated probability of the individual to be in its group. The probabilities \code{p} are often estimated by a logistic regression in which the dependent binary variable is the group. The possible confounding factors are the explanatory variables of this model. 
}
\value{
\item{statistic}{The value of the log-rank statistic.}
\item{p.value}{The p-value associated to the previous log-rank statistic.}
}
\references{
Le Borgne et al. Comparisons of the performances of different statistical tests for time-to-event analysis with confounding factors: practical illustrations in kidney transplantation. Statistics in medicine. 30;35(7):1103-16, 2016. <doi:10.1002/ sim.6777>

Jun Xie and Chaofeng Liu. Adjusted Kaplan-Meier estimator and log-rank test with inverse probability of treatment weighting for survival data. Statistics in medicine, 24(20):3089-3110, 2005. <doi:10.1002/sim.2174>
}

\author{
Yohann Foucher <Yohann.Foucher@univ-poitiers.fr>

Jun Xie <junxie@purdue.edu>

Florant Le Borgne <fleborgne@idbc.fr>
}

\examples{

data(dataDIVAT2)

# adjusted log-rank test
Pr0 <- glm(ecd ~ 1, family = binomial(link="logit"), data=dataDIVAT2)$fitted.values[1]
Pr1 <- glm(ecd ~ age + hla + retransplant, data=dataDIVAT2,
 family=binomial(link = "logit"))$fitted.values
W <- (dataDIVAT2$ecd==1) * (1/Pr1) + (dataDIVAT2$ecd==0) * (1)/(1-Pr1)

ipw.log.rank(dataDIVAT2$times, dataDIVAT2$failures, dataDIVAT2$ecd, W)

}

\keyword{Inverse Probability Weighting}
