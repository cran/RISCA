\name{roc.net}
\alias{roc.net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Net Time-Dependent ROC Curves With Right Censored Data.
}
\description{
This function performs the characteristics of a net time-dependent ROC curve (Lorent, 2013) based on k-nearest neighbor's (knn) estimator or only based on the Pohar-Perme estimator (Pohar, 2012).
}
\usage{
roc.net(times, failures, variable, p.age, p.sex, p.year,
 rate.table, pro.time, cut.off, knn=FALSE, prop=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{A numeric vector with the follow up times.}
  \item{failures}{A numeric vector with the event indicator (0=right censored, 1=event).}
  \item{variable}{A numeric vector with the prognostic variable. This variable is collected at the baseline.}
  \item{p.age}{A numeric vector with the age of the individuals at the baseline (in days).}
  \item{p.sex}{A character vector with the gender the individuals ('male' or 'female').}
  \item{p.year}{A numeric vector with the calendar year at the baseline (number of days from the January 1, 1960).}
  \item{rate.table}{A rate-table object with the expected mortality rates by age, sex, and cohort year. The same units used in \code{p.age}, \code{p.sex}, \code{p.year}}
  \item{pro.time}{The value of prognostic time represents the maximum delay for which the capacity of the variable is evaluated. The same unit than the one used in the argument \code{time}.}
  \item{cut.off}{The cut-off values of the variable used to define the possible binary tests.}
  \item{knn}{A logical value indicating whether k-nearest neighbor's estimator should be used.}
  \item{prop}{This is the proportion of the nearest neighbors. The estimation will be based on 2*prop (both right and left proportions) of the total sample size.}
}
\details{
This function computes net time-dependent ROC curve with right-censored data using estimator defined by Pohar-Perm et al. (2012) and the k-nearest neighbor's (knn) estimator. The aim is to evaluate the capacity of a variable (measured at the baseline) to predict the excess of mortality of a studied population compared to the general population mortality. Using the knn estimator ensures a monotone and increasing ROC curve, but the computation time may be long. This approach may thus be avoided if the sample size is large because of computing time.
}
\value{
\item{table}{This data frame presents the sensitivities and specificities associated with the cut-off values. One can observe NA if the value cannot be computed. \code{J} represents the Youden index.}
\item{auc}{The area under the time-dependent ROC curve for a prognostic up to prognostic time.}
\item{missing}{Number of deleted observations due to missing data.}
\item{warning}{This message indicates possible difficulties in the computation of the net ROC curve, for instance if the net survival was not lower or equal1 to 1 for particular cut-off values or times.}
}
\references{
Lorent et al. Net time-dependent ROC curves: a solution for evaluating the accuracy of a marker to predict disease-related mortality. Statistics in Medicine, 33, 2379-89. 2013. <doi:10.1002/sim.6079>
}
\author{
Y. Foucher <Yohann.Foucher@univ-poitiers.fr>
}

\examples{

# import the observed data

data(dataDIVAT3)

# A subgroup analysis to reduce the time needed for this example

dataDIVAT3 <- dataDIVAT3[1:400,]

# import the expected mortality rates

data(fr.ratetable)

# the values of recipient age used for computing the sensibilities and
# specificities (choose more values in practice)

age.cut <- quantile(dataDIVAT3$ageR, probs=seq(0.1, 0.9, by=0.1))

# recoding of the variables for matching with the ratetable

dataDIVAT3$sex <- "male"
dataDIVAT3$sex[dataDIVAT3$sexeR==0] <- "female"
dataDIVAT3$year <- date::mdy.date(month=01, day=01, year=dataDIVAT3$year.tx,
 nineteen = TRUE, fillday = FALSE, fillmonth = FALSE)
dataDIVAT3$age <- dataDIVAT3$ageR*365

# the ROC curve (without correction by the knn estimator) to 
# reduce the time for computing this example. In practice, the
# correction should by used in case of non-montone results.

roc1 <- roc.net(times=dataDIVAT3$death.time,
 failures=dataDIVAT3$death,  variable=dataDIVAT3$ageR,
 p.age=dataDIVAT3$age, p.sex=dataDIVAT3$sex, p.year=dataDIVAT3$year,
 rate.table=fr.ratetable, pro.time=3000, cut.off=age.cut, knn=FALSE)
 
# the sensibilities and specificities associated with the cut off values

roc1$table

# the traditional ROC graph

plot(roc1, col=2, pch=2, lty=2, type="b", xlab="1-specificity", ylab="sensibility")

legend("bottomright", paste("Without knn, (AUC=",
 round(roc1$auc, 2), ")", sep=""),lty=1, lwd=2, col=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ROC}
\keyword{Survival}