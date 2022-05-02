\name{survival.mr}
\alias{survival.mr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiplicative-Regression Model to Compare the Risk Factors Between Two Reference and Relative Populations
}
\description{
Compute a multiplicative-regression model to compare the risk factors between a reference and a relative population.
}
\usage{
survival.mr(times, failures, cov.relative, data,
cox.reference, cov.reference, ini, iterations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{The column name in \code{data}, in which the time of follow-up of each individual is collected.}
  \item{failures}{The column name in \code{data}, in which the indicator of event at the end of follow-up is collected (1 if the event is observed and 0 if right censoring).}
  \item{cov.relative}{The column(s) name(s) in \code{data} in order to declare the explicative variable included in the multiplicative relative model.}
  \item{data}{A data frame with the variables (columns) of the individuals (raw) of the relative sample.}
  \item{cox.reference}{The results of the Cox model performed in the reference sample, i.e an object obtained by the \code{coxph} function.}
  \item{cov.reference}{The column(s) name(s) in \code{data} in order to declare the explicative variable corresponding to those included in the Cox model \code{cox.reference}. Please, note that the order of these variables is important and have to be similar with the order in \code{cox.reference}.}
  \item{ini}{A vector with the same length than \code{cov.relative} with the initial values for the parameters to be optimized.}
  \item{iterations}{The number of iterations of the bootstrap resampling.}
  }
\details{
We proposed here an adaptation of a multiplicative-regression model for relative survival to study the heterogeneity of risk factors between two groups of patients. Estimation of parameters is based on partial likelihood maximization and Monte-Carlo simulations associated with bootstrap re-sampling yields to obtain the corresponding standard deviations.  The expected hazard ratios are obtained by using a PH Cox model. 
}
\value{
\item{matrix.coef}{A matrix containing the parameters estimations at each of the B iterations.}
\item{estim.coef}{A numerical vector containing the mean of the previous estimation}
\item{lower95.coef}{A numerical vector containing the lower bounds of the 95\% confidence intervals.}
\item{upper95.coef}{A numerical vector containing the upper bounds of the 95\% confidence intervals.}
}
\references{
K. Trebern-Launay, M. Giral, J. Dantal and Y. Foucher. Comparison of the risk factors effects between two populations: two alternative approaches illustrated by the analysis of first and second kidney transplant recipients. BMC Med Res Methodol. 2013 Aug 6;13:102. <doi: 10.1186/1471-2288-13-102>
}
\author{
Y. Foucher <Yohann.Foucher@univ-nantes.fr>

K. Trebern-Launay <katygre@yahoo.fr>
}

\examples{

# import and attach both samples
data(dataFTR)
data(dataSTR)

# We reduce the dimension to save time for this example (CRAN policies)
# Compute the Cox model in the First Kidney Transplantations (FTR)
cox.FTR<-coxph(Surv(Tps.Evt, Evt)~ ageR2cl + sexeR, data=dataFTR[1:100,])
summary(cox.FTR)

# Compute the multiplicative relative model
# for Second Kidney Transplantations  (STR)
# Choose iterations>>5 for real applications
mrs.STR <- survival.mr(times="Tps.Evt", failures="Evt",
 cov.relative=c("ageR2cl", "Tattente2cl"), data=dataSTR[1:100,],
 cox.reference=cox.FTR, cov.reference=c("ageR2cl", "sexeR"),
 ini=c(0,0), iterations=5)

  
# The parameters estimations (mean of the values)
mrs.STR$estim.coef 

# The 95 percent. confidence intervals
cbind(mrs.STR$lower95.coef, mrs.STR$upper95.coef) 
}

\keyword{Survival}

