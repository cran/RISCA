\name{dataCSL}
\alias{dataCSL}
\docType{data}
\title{
CSL Liver Chirrosis Data.
}
\description{
Survival status for the liver chirrosis patients of Schlichting et al.
}
\usage{data(dataCSL)}
\format{
This data frame contains the following columns:
\tabular{rl}{
\code{id} \tab a numeric vector related to the subject ID. \cr									
\code{time} \tab A numeric vector with the time of measurement. \cr
\code{prot} \tab A numeric vector with the prothrombin level at measurement time. \cr
\code{dc} \tab A numeric vector code. 0: censored observation, 1: died at eventT. \cr
\code{eventT} \tab A numeric vector with the time of event (death). \cr
\code{treat} \tab A numeric vector code. 0: active treatment of prednisone, 1: placebo treatment.  \cr
\code{sex} \tab A numeric vector code. 0: female, 1: male.  \cr
\code{age} \tab A numeric vector with the age of subject at inclusion time subtracted 60.  \cr
\code{prot.base} \tab  A numeric vector with the prothrombin base level before entering the study \cr
\code{prot.prev} \tab A numeric vector with the level of prothrombin at previous measurement time.  \cr
\code{lt} \tab  A numeric vector with the starting time for the time-interval. \cr
\code{rt} \tab  A numeric vector with the stopping time for the time-interval.}
}

\source{
The \code{timreg} package
}
\references{
Schlichting et al. Prognostic factors in cirrhosis identified by Cox's regression model. Hepatology. 1983 Nov-Dec;3(6):889-95. <doi https://doi.org/ 10.1002/ hep.1840030601>
}
\examples{
data(dataCSL)
names(dataCSL)
}

\keyword{Dataset}
