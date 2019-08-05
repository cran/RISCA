\name{rein.ratetable}
\alias{rein.ratetable}
\docType{data}
\title{
Expected Mortality Of French Patients With ESKD.
}
\description{
An object of the class rate-table which contains the expected mortality of French patients with End Stage Kidney Disease (ESKD) in dialysis and registered previously on waiting list for renal transplantation. It is an array with two dimensions: the time since the registry on waiting list and the patient profile at this baseline.
}
\usage{data(rein.ratetable)}
\format{
 The format is "ratetable". The attributes are:
 \describe{
 \item{\code{dim}}{A numeric vector with the length of each dimension.}
 \item{\code{dimnames}}{A list of vectors with the names of each variable of both dimensions (see details for modalities of \code{profile}).}
 \item{\code{dimid}}{A character vector with the identification of the dimensions: \code{time}, \code{profile}.}
 \item{\code{factor}}{A vector of indicators=1 if the corresponding dimension do not vary according to the \code{time}. The \code{profile} is associated to 1.}
 \item{\code{cutpoints}}{A list of the thresholds to identify the mortality rates according to the time-varying dimensions (NULL for \code{profile}).}
 \item{\code{class}}{The class of the object: \code{ratetable}.}
 }
}
\details{
The organization of a rate-table object is described in details by Therneau (1999) and Pohar (2006). We have estimated these data by using the renal epidemiology and information network (REIN). A competing risk model was performed in order to deal with the loss of follow-up to the transplantation. The \code{profile} represents the concatenation of the significant binary variables associated with the mortality on waiting list. In order to obtain the same variable in her/his database, the user have to concatenate the following values in the same order. These values are collected at the registration date on the waiting list:

\tabular{rl}{
a. \tab 1 if the patient is male and 0 if female. \cr
b. \tab 1 if the time between the first dialysis and the registration is higher than 1 year and 0 otherwise. \cr
c. \tab 1 if hemodialysis and 0 if peritoneal dialysis. \cr
d. \tab 1 if the patient is diabetic and 0 otherwise. \cr
e. \tab 1 if the patient had history of cardiac failure diabetic and 0 otherwise. \cr
f. \tab 1 if the patient had history of vascular failure diabetic and 0 otherwise. \cr
g. \tab 1 if the patient had history of cancer and 0 otherwise. \cr
h. \tab 1 if the age of the patients is in between 40 and 60 years and 0 otherwise. \cr
i. \tab 1 if the age of the patients is higher than 60 years and 0 otherwise. \cr
j. \tab 1 if the calendar year at registration is strictly higher than 2007 and 0 \cr
}
}
\source{
URL: www.agence-biomedecine.fr/Le-programme-REIN
}
\references{
Pohar and Stare Relative survival analysis in R. Computer methods and programs in biomedicine, 81, 272-278, 2006.  <doi:10.1016/j.cmpb.2006.01.004>
}
\examples{
data(rein.ratetable)

is.ratetable(rein.ratetable)

}

\keyword{Ratetable}

