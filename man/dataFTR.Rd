\name{dataFTR}
\alias{dataFTR}
\docType{data}
\title{
Data for First Kidney Transplant Recipients.
}
\description{
Data were extracted from the DIVAT cohort. It corresponds to the reference sample constituted by first transplant recipients (FTR).  
}
\usage{data(dataFTR)}
\format{
A data frame with the 4 following variables:
  \describe{
    \item{\code{Tps.Evt}}{This numeric vector provides the post-transplantation time (in days).}
	\item{\code{Evt}}{This numeric vector provides the indicator of graft failure at the end of the follow-up (1 for failure and 0 for right censoring).}
    \item{\code{ageR2cl}}{This numeric vector provides the recipient age at transplantation (1 for older than 55 years and 0 otherwise).}
    \item{\code{sexeR}}{This numeric vector provides the recipient gender (1 for men and 0 for women).}
  }
}
\details{
First transplant recipients (FTR) constituted the reference group. Recipients older than 18 years at the date of transplantation between 1996 and 2010 were selected from the French DIVAT (www.divat.fr/en) multicentric prospective cohort. Only recipients with a maintenance therapy with calcineurin inhibitors, mammalian target of rapamycin inhibitors or belatacept, in addition to mycophenolic acid and steroid were included. Simultaneous transplantations were excluded. Two explicative variables are proposed: recipient age at transplantation and recipient gender.  
}
\references{
K. Trebern-Launay, M. Giral, J. Dantal and Y. Foucher. Comparison of the risk factors effects between two populations: two alternative approaches illustrated by the analysis of first and second kidney transplant recipients. BMC Med Res Methodol. 2013 Aug 6;13:102.  <doi: 10.1186/1471-2288-13-102>
}

\examples{

data(dataFTR)

# Compute a Cox PH model with both explicative variables
summary(coxph(Surv(Tps.Evt/365.24, Evt) ~ ageR2cl + sexeR, data=dataFTR))
}

\keyword{Dataset}

