\name{dataDIVAT5}
\alias{dataDIVAT5}
\docType{data}
\title{
The Aggregated Kidney Graft Survival Stratified By The 1-year Serum Creatinine.
}
\description{
A data frame that presents the aggregated outcomes per center of 5943 French kidney transplant recipients from the DIVAT cohort.  It was used by Combescure et al. (2016).
}
\usage{data(dataDIVAT5)}
\format{
A data frame with 106 observations for the 8 following variables:
  \describe{
    \item{\code{classe}}{This numeric vector represents the groups of recipients defined using the 1-year serum creatinine. 1 is the first group with the lowest values of 1-year serum creatinine.}
	\item{\code{n}}{This numeric vector represents the number of recipients at the baseline (date of the transplantation) in each group.}
    \item{\code{year}}{This numeric vector represents the post transplant time (in years).}
    \item{\code{surv}}{This numeric vector represents the survival probabilities at each year (obtained using the Kaplan and Meier estimator from the individual data).}
    \item{\code{n.risk}}{This numeric vector represents the number of subjects at-risk of the event at the corresponding \code{year}.}
    \item{\code{proba}}{This numeric vector represents the proportion of the patients in a center which belong to the corresponding group.}
    \item{\code{marker.min}}{This numeric vector represents the minimum value of the interval of the 1-year serum creatinine (in \eqn{\mu mol/l}).}
    \item{\code{marker.max}}{This numeric vector represents the maximum value of the interval of the 1-year serum creatinine (in \eqn{\mu mol/l}).}
	\item{\code{centre.num}}{This numeric vector represents the centers.}
  }
}
\details{
The immunology and nephrology department of the Nantes University hospital constituted a data bank with the monitoring of medical records for kidney and/or pancreas transplant recipients. Here, we considered a subpopulation of 4195 adult patients and who had received a first kidney graft between January 1996 and Jun 2008. Five centers participated. A total of 511 graft failures were observed (346 returns to dialysis and 165 deaths with a functional kidney). Based on this database, we constructed an aggregated dataset to perform a meta-analysis on 5 published monocentric studies. The medical objective was to evaluate whether 1-year serum creatinine (Cr) is a good predictive marker of graft failure. Cr is a breakdown product and is removed from the body by the kidneys. If kidney function is abnormal, blood Cr levels increase.
}
\source{
URL: www.divat.fr
}
\references{
Combescure et al. A literature-based approach to evaluate the predictive capacity of a marker using time-dependent Summary Receiver Operating Characteristics. Stat Methods Med Res, 25(2):674-85, 2016. <doi: 10.1177/ 0962280212464542>
}
\examples{
data(dataDIVAT5)

# Kaplan Meier estimations of the graft survival in the first center
plot(dataDIVAT5$year[dataDIVAT5$centre.num==1],
 dataDIVAT5$surv[dataDIVAT5$centre.num==1],
 xlab="Post transplantation time (years)", ylab="Graft survival",
 ylim=c(0.7,1), xlim=c(0, 9), type="n")

# Goup 1
lines(c(0, dataDIVAT5$year[dataDIVAT5$centre.num==1 &
 dataDIVAT5$classe==1]),
 c(1, dataDIVAT5$surv[dataDIVAT5$centre.num==1 &
 dataDIVAT5$classe==1]),
 type="b", col=1, lty=1, lwd=2)

# Goup 2
lines(c(0, dataDIVAT5$year[dataDIVAT5$centre.num==1 &
 dataDIVAT5$classe==2]),
 c(1, dataDIVAT5$surv[dataDIVAT5$centre.num==1 &
 dataDIVAT5$classe==2]),
 type="b", col=2, lty=2, lwd=2)

# legend
legend("bottomleft", c("group #1 (1-year Cr<4.57)",
 "group #2 (1-year Cr>4.57)"), col=c(1, 2),
 lty=c(1, 2), lwd=c(2, 2))
}

\keyword{Dataset}


