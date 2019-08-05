\name{dataHepatology}
\alias{dataHepatology}
\docType{data}
\title{
The Data Extracted From The Meta-Analysis By Cabibbo et al. (2010).
}
\description{
Data were extracted from the studies included in the meta-analysis by Cabibbo et al. which aimed to assess the survival rate in untreated patients with hepatocellular carcinoma.
}
\usage{data(dataHepatology)}
\format{
A data frame with with the 8 following variables (columns).
  \describe{
    \item{\code{study}}{This numeric vector represents number of the study.}
	\item{\code{first.author}}{This vector represents the name of the first author.}
    \item{\code{year.pub}}{This numeric vector represents the publication year.}
    \item{\code{time}}{This numeric vector represents the times for which the survival rates are collected in years.}
    \item{\code{survival}}{This numeric vector represents the survival rates for each value of \code{time}}
    \item{\code{n.risk}}{This numeric vector represents the number of at-risk patients for each value of \code{time}}
    \item{\code{location}}{This factor indicates the location of the study (Asia, North Amercia or Europe)}
    \item{\code{design}}{This factor indicates if the study is monocentric ou multicentric.}
  }
}
\details{
The survival probabilities were extracted from the published survival curves each month during the first six months and then by step of three months. The pictures of the curves were digitalized using the R package ReadImage and the probabilities were extracted using the package digitize proposed by Poisot. The numbers of at-risk patients for each interval of time were derived from the numbers of at-risk patients reported in the studies, and using the methods of Parmar or Williamson to account for censorship. Studies have different length of follow-up. For each study, survival probabilities and the numbers of at-risk patients were collected at all points in time before the end of follow-up.
}
\references{
Cabibbo, G., et al., A meta-analysis of survival rates of untreated patients in randomized clinical trials of hepatocellular carcinoma. Hepatology, 2010. 51(4): p. 1274-83.

Poisot, T., The digitize Package: Extracting Numerical Data from Scatter-plots. The R Journal, 2011. 3(1): p. 25-26.

Parmar, M.K., V. Torri, and L. Stewart, Extracting summary statistics to perform meta-analyses of the published literature for survival endpoints. Stat Med, 1998. 17(24): p. 2815-34
}
\examples{

data(dataHepatology)
times <- dataHepatology$time
survival <- dataHepatology$survival
study  <- dataHepatology$study

plot(times, survival, type="n", 
 ylim=c(0,1), xlab="Time",ylab="Survival")

for (i in unique(sort(study)))
{
lines(times[study==i], survival[study==i], type="l", col="grey")
points(max(times[study==i]),
 survival[study==i & times == max( times[study==i])], pch=15)
}
}

\keyword{Dataset}

