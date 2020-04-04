\name{gc.sl.binary}
\alias{gc.sl.binary}
\title{
Marginal Effect for Binary Outcome by Super Learned G-computation.
}
\description{
This function allows to estimate the marginal effect of an exposure or a treatment by G-computation for binary outcomes, the outcome prediction (Q-model) is obtained by using a super learner.
}
\usage{
gc.sl.binary(outcome, group, cov.quanti, cov.quali, data, effect,
 tuneLength, cv, iterations, n.cluster)
}

\arguments{
  \item{outcome}{The name of the variable related to the binary outcome variable. This numeric variable must have only two levels: 0 and 1 (for instance: 0=healthy, 1=diseased).}
  \item{group}{The name of the variable related to the exposure/treatment variable. This numeric variable must have only two levels: 0 for the untreated/unexposed patients and 1 for the treated/exposed patients.}
  \item{cov.quanti}{The name(s) of the variable(s) related to the possible quantitative covariates. These variables must be numeric.}
    \item{cov.quali}{The name(s) of the variable(s) related to the possible qualitative covariates. These variables must be numeric with two levels: 0 and 1. A complete disjunctive form have two be used for covariates with more levels.}
  \item{data}{A data frame in which to look for the variable related to the outcome, the treatment/exposure and the covariables.}
 \item{effect}{The type of the marginal effect to be estimated. Three types are possible (see details): "ATE" (by default), "ATT" and "ATU".}
 \item{tuneLength}{It defines the total number of parameter combinations that will be evaluated in the machine learning techniques.} 
 \item{cv}{The number of splits for cross-validation.}
 \item{iterations}{The number of bootstrap resamples to estimate of the variances and confidence intervals.}
  \item{n.cluster}{The number of cores to use, i.e., the upper limit for the number of child processes that run simultaneously (1 by default).}
 }
 
\details{
The ATE corresponds to Average Treatment effect on the Entire population, i.e. the marginal effect if all the sample is treated versus all the sample is untreated. The ATT corresponds to Average Treatment effect on the Treated, i.e. the marginal effect if the treated patients (\code{group = 1}) would have been untreated. The ATU corresponds to Average Treatment effect on the Untreated , i.e. the marginal effect if the untreated patients (\code{group = 0}) would have been treated. The Super Learner includes the following machine learning techniques: logistic regression with Lasso penalization, logistic regression with Elasticnet penalization, neural network with one hidden layer, and support vector machine with radial basis, as explained in details by Chatton et al. (2020).
}

\value{
\item{misssing}{Number of deleted observations due to missing data.}
\item{effect}{A character string with the type of the marginal effect.}
\item{p0}{A table related to the average proportion of events in the unexposed/untreated sample: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{p1}{A table related to the average proportion of events in the exposed/treated sample: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{delta}{A table related to the difference between the average proportions of events in the exposed/treated sample minus in the unexposed/untreated sample: \code{estimate} is the estimated value, \code{std.error} is the corresponding standard error, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{logOR}{A table related to the logarithm of the average Odds Ratio (OR): \code{estimate} is the estimated value, \code{std.error} is the corresponding standard error, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{p.value}{The p-value of the bilateral test of the null hypothesis \code{p0} = \code{p1}, i.e. OR = 1.}
}

\references{
Chatton A et al. Improving Causal Inference Through the Super Learned G-computation. Manuscript submitted. 2020.
}

\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>
}

\examples{
#data simulation

#treatment = 1 if the patients have been the exposure or treatment of interest and 0 otherwise
#treatment <- rbinom(200, 1, prob=0.5)

#one quantitative covariate

#covariate1 <- rnorm(200, 0, 1)
#covariate1[treatment==1] <- rnorm(sum(treatment==1), 0.3, 1)

#one qualitative covariate

#covariate2 <- rbinom(200, 1, prob=0.5)
#covariate2[treatment==1] <- rbinom(sum(treatment==1), 1, prob=0.4)

#outcome <- rbinom(200, 1, prob = exp(-2+0.26*treatment+0.4*covariate1-

#0.4*covariate2)/(1+exp( -2+0.26*treatment+0.4*covariate1-0.4*covariate2)))

#tab <- data.frame(outcome, treatment, covariate1, covariate2)

#Raw effect of the treatment

#glm.raw <- glm(outcome ~ treatment, data=tab, family = binomial(link=logit))
#summary(glm.raw)

#Conditional effect of the treatment

#glm.multi <- glm(outcome ~ treatment + covariate1 + covariate2,
# data=tab, family = binomial(link=logit))
#summary(glm.multi)

#Marginal effects of the treatment (ATE) by using logistic regression as the Q-model

#gc.ate1 <- gc.logistic(glm.obj=glm.multi, data=tab, group="treatment", effect="ATE",
# var.method="bootstrap", iterations=1000, n.cluster=1)

#Marginal effects of the treatment (ATE) by using a super learner as the Q-model

#gc.ate2 <- gc.sl.binary(outcome="outcome", group="treatment", cov.quanti="covariate1",
# cov.quali="covariate2", data=tab, effect="ATE", tuneLength=10, cv=3,
# iterations=1000, n.cluster=1)

#Sum-up of the 3 ORs

#data.frame( raw=exp(glm.raw$coefficients[2]),
#conditional=exp(glm.multi$coefficients[2]),
#marginal.ate.logistic=exp(gc.ate1$logOR[,1]),
#marginal.ate.sl=exp(gc.ate2$logOR[,1]) )

}

\keyword{G-computation}
