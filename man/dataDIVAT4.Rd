\name{dataDIVAT4}
\alias{dataDIVAT4}
\docType{data}
\title{
A Fourth Sample From the DIVAT Data Bank.
}

\description{
A data frame with 6648 French kidney transplant recipients from the DIVAT cohort. According to this data set, patient and graft survival can be computed for each hospital. This database was used by Combsecure et al. (2014).
}

\usage{data(dataDIVAT4)}

\format{
A data frame with the 3 following variables:
  \describe{
    \item{\code{hospital}}{This numeric vector represents the hospital in which the kidney transplantation was performed. Six hospitals were included.}
	\item{\code{time}}{For right censored data, this numeric vector represents the follow up time (in days). When the event is observed, this numeric vector represents the exact time-to-event.}
    \item{\code{status}}{This indicator describes the end of the follow-up. The \code{status} equals 0 for right censored data and 1 if the event is observed.}
  }
}

\details{
The data were extracted from the prospective multicentric DIVAT cohort (www.divat.fr). Patients transplanted between 2000 and 2012 in 6 French centers were included. Only adult patients receiving a single kidney transplant and treated by Calcineurin inhibitors  and Mycophenolate Mofetil for maintenance therapy after transplantation were considered. In parallel, patients who received multi-organ transplantation were excluded from the study.  The time-to-event under interest is the time between the kidney transplantation and the graft failure, i.e. the first event between the return in dialysis or the death of the patient with functional kidney.
}

\references{
Combescure et al. The multivariate DerSimonian and Laird's methodology applied to meta-analysis of survival curves. 10;33(15):2521-37, 2014. Statistics in Medicine. <doi:10.1002/sim.6111>
}


\source{URL: www.divat.fr}

\examples{
data(dataDIVAT4)

divat.surv <- survfit(Surv(time/365.24, status) ~ hospital, data = dataDIVAT4) 

plot(divat.surv, lty = 1:6, col=1:6, lwd=2, mark.time=FALSE,
 xlab="Post transplantation time (days)", ylab="Patient and graft survival") 
}

\keyword{Dataset}

