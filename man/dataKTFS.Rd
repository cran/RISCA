\name{dataKTFS}
\alias{dataKTFS}
\docType{data}
\title{
A Fifth Sample Of The DIVAT Cohort.
}
\description{
A data frame with 2169 French kidney transplant recipients for who the Kidney Transplant Failure Score (KTFS) was collected. The KTFS is a score proposed by Foucher et al. (2010) to stratify the recipients according to their risk of return in dialysis.
}
\usage{data(dataKTFS)}
\format{
  A data frame with 2169 observations (raws) with the 3 following variables (columns).
  \describe{
    \item{\code{time}}{This numeric vector represents the follow up time in years (until return in dialyis or censoring)}
    \item{\code{failure}}{This numeric vector represents the graft failure indicator at the follow-up end (1=return, 0=censoring)}
    \item{\code{score}}{This numeric vector represents the KTFS values.}
  }
}

\references{
Foucher Y. al. A clinical scoring system highly predictive of long-term kidney graft survival. Kidney International, 78:1288-94, 2010. <doi:10.1038/ki.2010.232>
}
\source{ URL: www.divat.fr }
\examples{
data(dataKTFS)

### a short summary of the recipient age at the transplantation
summary(dataKTFS$score)

### Kaplan and Meier estimation of the recipient survival
plot(survfit(Surv(time, failure) ~ I(score>4.17), data = dataKTFS),
 xlab="Post transplantation time (in years)", ylab="Patient survival",
 mark.time=FALSE, col=c(2,1), lty=c(2,1))
 
legend("bottomleft", c("Recipients in the high-risk group",
 "Recipients in the low-risk group"), col=1:2, lty=1:2)
}

\keyword{Dataset}

