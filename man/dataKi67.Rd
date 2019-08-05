\name{dataKi67}
\alias{dataKi67}
\docType{data}
\title{
The Aggregated Data Published By de Azambuja et al. (2007).
}
\description{
The aggregated data from the meta-analysis proposed by Azambuja et al. (2007).
}
\usage{data(dataKi67)}
\format{
A data frame with 406 observations (rows) with the 10 following variables (columns).
  \describe{
    \item{\code{classe}}{This numeric vector represents the groups of patients defined using KI-67. 1 is the first group which is defined by the lowest KI-67 values.}
	\item{\code{n}}{This numeric vector represents the number of recipients at the baseline (date of KI-67 collection) in each group.}
    \item{\code{year}}{This numeric vector represents the survival time (in years).}
    \item{\code{surv}}{This numeric vector represents the survival probabilities at each year (obtained using the Kaplan and Meier estimator from the published papers).}
    \item{\code{nrisk}}{This numeric vector represents the number of subjects at-risk of the event at the corresponding \code{year}.}
    \item{\code{proba}}{This numeric vector represents the proportion of the patients for a given paper which belong to the corresponding group.}
    \item{\code{log.marker.min}}{This numeric vector represents the logarithm of the minimum value of the KI-67 interval.}
    \item{\code{log.marker.max}}{This numeric vector represents the logarithm of the maximum value of the KI-67 interval.}
    \item{\code{study.num}}{This numeric vector identifies the studies.}
    \item{\code{author}}{This character vector identifies the first author of the paper.}
    \item{\code{year.paper}}{This numeric vector identifies the year of publication.}

  }
}
\details{
KI-67 is a marker of the proliferative activity of breast cancer, but its prognostic capacity is still unclear. In their meta-analysis, de Azambuja et al. (2007) concluded that KI-67 positivity conferred a worse survival. This work focused on the 35 evaluable studies of the relationship between KI-67 and the overall survival. 23 studies described survival curves according to the level of KI-67. Survival probabilities were measured every year.
}
\references{
de Azambuja et al. Ki-67 as prognostic marker in early breast cancer: a meta-analysis of published studies involving 12 155 patients. British Journal of Cancer. 96:1504-1513, 2007. <doi: 10.1038/ sj.bjc.6603756>
}

\examples{
data(dataKi67)

# Kaplan Meier estimations of graft survivals in Wintzer et al. (1991)
plot(dataKi67$year[dataKi67$study.num==1],
 dataKi67$surv[dataKi67$study.num==1],
 xlab="Post transplantation time (years)",
 ylab="Graft survival", ylim=c(0.6,1), xlim=c(0, 4), type="n")

# Goup 1
lines(c(0, dataKi67$year[dataKi67$study.num==1 & dataKi67$classe==1]),
 c(1, dataKi67$surv[dataKi67$study.num==1 & dataKi67$classe==1]),
 type="b", col=1, lty=1, lwd=2)

# Goup 2
lines(c(0, dataKi67$year[dataKi67$study.num==1 & dataKi67$classe==2]),
 c(1, dataKi67$surv[dataKi67$study.num==1 & dataKi67$classe==2]),
 type="b", col=2, lty=2, lwd=2)

# legend
legend("bottomleft", c("group #1 (log Ki67 < 2.49)",
 "group #2 (log Ki67 > 2.49)"), col=c(1, 2), lty=c(1, 2), lwd=c(2, 2))
}

\keyword{Dataset}

