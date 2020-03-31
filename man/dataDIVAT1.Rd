\name{dataDIVAT1}
\alias{dataDIVAT1}
\docType{data}
\title{
A First Sample From The DIVAT Data Bank.
}
\description{
A data frame with 5943 French kidney transplant recipients from the DIVAT cohort. 
}
\usage{data(dataDIVAT1)}
\format{
 A data frame with 5943 observations for the 7 following variables:
  \describe{
\item{\code{trajectory}}{A numeric vector with the sequences of observed states. The patient evolution can 
 be described according to a 4-state structure: X=1 represents the healthy state, X=2 represents the acute rejection episode, X=3 the definitive return to dialysis and X=4 the death. These times can be right-censored. A vector of covariates is also collected at the transplantation, i.e. the baseline of the cohort.}
\item{\code{time1}}{A numeric vector with the times (in days) between the transplantation and the first clinical event (acute rejection episode, return to dialysis, or death with a functioning graft), or the times to censoring if \code{trajectory=1}.}
\item{\code{time2}}{A numeric vector with the time between the transplantation and the second clinical event (return to dialysis or death with a functioning graft), or the the time to censoring if \code{trajectory=12}.}
\item{\code{ageR}}{A numeric vector with the recipient age (in years) at the transplantation.}
\item{\code{sexR}}{A character vector with the recipient gender.}
\item{\code{year.tx}}{A numeric vector with the calendar year of the transplantation.}
\item{\code{z}}{A numeric vector represents the explicative variable under interest, i.e. the delayed graft function (1=yes, 0=no).}
}
}

\source{
URL: www.divat.fr
}
\examples{
data(dataDIVAT1)

### a description of transitions
table(dataDIVAT1$trajectory)

### patient-graft survival (first event between the return to dialysis and the patient 
### death with a functioning graft)

dataDIVAT1$failure<-1*(dataDIVAT1$trajectory!=1 & dataDIVAT1$trajectory!=12)

dataDIVAT1$time<-NA
dataDIVAT1$time<-ifelse(dataDIVAT1$trajectory \%in\% c(1,12,13,14),
dataDIVAT1$time1,dataDIVAT1$time1+dataDIVAT1$time2)

plot(survfit(Surv(time/365.24, failure) ~ 1 , data=dataDIVAT1), mark.time=FALSE, 
      xlim=c(0,12), ylim=c(0,1), cex=1.5, col=1, lwd=2, lty=1, 
      xlab="Times after the transplantation (years)", 
      ylab="Patient-graft survival")

}

\keyword{Dataset}
