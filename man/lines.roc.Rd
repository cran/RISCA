\name{lines.roc}
\alias{lines.roc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add Lines to a ROC Plot
}
\description{
Used to add an additionnal ROC curve to ROC plot generated with \code{plot.roc}.
}
\usage{
\method{lines}{roc}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{roc}, returned by the functions \code{roc.binary},
  \code{roc.net}, \code{roc.summary}, and \code{roc.time}.}
  \item{...}{Additional arguments affecting the plot line.}
}
\references{
Le Borgne et al. Standardized and weighted time-dependent ROC curves to evaluate the intrinsic prognostic capacities of a marker by taking into account confounding factors. Manuscript submitted. Stat Methods Med Res. 27(11):3397-3410, 2018. <doi: 10.1177/ 0962280217702416.>
}
\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>
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

plot(roc2, type="b", col=2, pch=2, lty=2)

lines(roc1, type="b", col=1, pch=1)

legend("bottomright", lty=1:2, lwd=1, pch=1:2, col=1:2,
 c(paste("Crude estimation, (AUC=", round(roc1$auc, 2), ")", sep=""),
 paste("Adjusted estimation, (AUC=", round(roc2$auc, 2), ")", sep="") ) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Plotting}
