\name{roc.binary}
\alias{roc.binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ROC Curves For Binary Outcomes.
}
\description{
This function allows for the estimation of ROC curve by taking into account possible confounding factors. Two methods are implemented: i) the standardized and weighted ROC based on an IPW estimator, and ii) the placement values ROC.
}
\usage{
roc.binary(status, variable, confounders,  data, precision, estimator)
}
\arguments{
  \item{status}{A character string with the name of the variable in \code{data} which represents the disease status indicator (for instance: 0=healthy, 1=diseased).}
  \item{variable}{A character string with the name of the variable in \code{data} which represents the diagnostic/prognostic variable under interest. For the standardized and weighted ROC method, the variable must be previously standardized according to the covariates among the controls as proposed by Le Borgne et al. (2017).}
  \item{confounders}{An object of class "formula". More precisely only the right part with an expression of the form \code{~ model}, where \code{model} is the linear predictor of the logistic regressions performed for each cut-off value. The user can use \code{~1} to obtain the crude estimation.}
  \item{data}{An object of the class \code{data.frame} containing the variables previously detailed.}
  \item{precision}{A numeric vector with values between 0 and 1. The values represent the x-axis (1-specificity) of the ROC graph for which the user want to obtain the corresponding sensitivities. 0 and 1 are not allowed.}
  \item{estimator}{Two possible estimators can be used: "ipw" and "pv". IPW is based on the Inverse Probability Weigthing theory as proposed by Le Borgne et al. (2017). This estimator applied on a variable  standardized according to the covariates among the controls allows to obtain a standardized and weighted ROC. The IPW estimator is selected by default. The user can also use the placement values (pv) estimator as proposed by Pepe and Cai (Biometrics, 2004).}
}
\details{
This function computes confounder-adjusted ROC curve for uncensored data. We adapted the usual estimator by considering the individual probabilities to be diseased, given the possible confounding factors. The standardized and weighted ROC is obtained by both providing a variable under interest standardized according to the possible confounding factors and using "ipw" in the option \code{estimator}. The user can also use the estimator first proposed by Pepe and Cai (2004) which is based on placement values.
}
\value{
\item{table}{This data frame presents the sensitivities and specificities. \code{J} represents the Youden index.}
\item{auc}{The area under the ROC curve.}
}
\references{
Blanche et al. (2013) Review and comparison of roc curve estimators for a time-dependent outcome with marker-dependent censoring. Biometrical Journal, 55, 687-704. <doi:10.1002/bimj.201200045>

Pepe and Cai. (2004) The analysis of placement values for evaluating discriminatory measures. Biometrics, 60(2), 528-35. <doi:10.1111/j.0006-341X.2004.00200.x>

Le Borgne et al. Standardized and weighted time-dependent ROC curves to evaluate the intrinsic prognostic capacities of a marker by taking into account confounding factors. Stat Methods Med Res. 27(11):3397-3410, 2018. <doi: 10.1177/ 0962280217702416>.
}
\author{
Y. Foucher <Yohann.Foucher@univ-poitiers.fr>
}

\examples{

# import and attach the data example

data(dataDIVAT3)

# A subgroup analysis to reduce the time needed for this exemple

dataDIVAT3 <- dataDIVAT3[1:400,]

# The ROC curve to evaluate the crude capacities of the recipient age for the
# prognosis of post kidney transplant mortality (we ignore the censoring process)

roc1 <- roc.binary(status="death", variable="ageR", confounders=~1,
data=dataDIVAT3, precision=seq(0.1,0.9, by=0.1) )

# The corresponding ROC graph with basic options

plot(roc1, xlab="1-specificity", ylab="sensibility")

# The corresponding ROC graph with color and symbols

plot(roc1, col=2, pch=2, type="b", xlab="1-specificity", ylab="sensibility")

# The standardized and weighted ROC curve to evaluate the capacities
# of the recipient age for the prognosis of post kidney transplant
# mortality by taking into account the donor age and the recipient
# gender (we ignore the censoring process).

# 1. Standardize the marker according to the covariates among the controls
lm1 <- lm(ageR ~ ageD + sexeR, data=dataDIVAT3[dataDIVAT3$death == 0,])
dataDIVAT3$ageR_std <- (dataDIVAT3$ageR - (lm1$coef[1] + lm1$coef[2] *
 dataDIVAT3$ageD + lm1$coef[3] * dataDIVAT3$sexeR)) / sd(lm1$residuals)

# 2. Compute the sensitivity and specificity from the proposed IPW estimators
roc2 <- roc.binary(status="death", variable="ageR_std",
 confounders=~bs(ageD, df=3) + sexeR, data=dataDIVAT3, precision=seq(0.1,0.9, by=0.1))

# The corresponding ROC graph

plot(roc2, col=2, pch=2, lty=2, type="b", xlab="1-specificity", ylab="sensibility")

lines(roc1, col=1, pch=1, type="b")

legend("bottomright", lty=1:2, lwd=1, pch=1:2, col=1:2,
 c(paste("Crude estimation, (AUC=", round(roc1$auc, 2), ")", sep=""),
 paste("Adjusted estimation, (AUC=", round(roc2$auc, 2), ")", sep="") ) ) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ROC}
