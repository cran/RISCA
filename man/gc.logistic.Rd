\name{gc.logistic}
\alias{gc.logistic}
\title{
Marginal Effect for Binary Outcome by G-computation.
}
\description{
This function allows to estimate the marginal effect of an exposure or a treatment by G-computation for binary outcomes.
}
\usage{
gc.logistic(glm.obj, data, group, effect, var.method, iterations, n.cluster)
}

\arguments{
  \item{glm.obj}{A glm object obtained by using the function glm with the argument \code{family = binomial(link=logit)} to obtain a multivariate logistic regression. It shall be included the exposition/treatment of interest variable and at least one covariate.}
  \item{data}{A data frame in which to look for the variable related to the outcome, the treatment/exposure and the covariables included in the previous logistic regression  \code{glm.obj}. The covariates' names have to be identical than the ones included in the logistic model.}
  \item{group}{The name or the variable related to the exposure/treatment variable. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed patients.}
 \item{effect}{The type of the marginal effect to be estimated. Three types are possible (see details): "ATE" (by default), "ATT" and "ATU".}
 \item{var.method}{The method to estimate the variances and the confidence intervals. Two methods are possible: "simulations" (by default) which consists in parametric simulation based on the maximum likelihood estimates of the multivariate logistic regression \code{glm.obj} and "bootstrap" which consists in bootstrap resampling of the database \code{data}.} 
 \item{iterations}{The number of iterations (simulations or resamples depending on the argument in \code{var.method}) to estimate of the variances and confidence intervals.}
  \item{n.cluster}{The number of cores to use, i.e., the upper limit for the number of child processes that run simultaneously (1 by default).}
 }
 
\details{
The ATE corresponds to Average Treatment effect on the Entire population, i.e. the marginal effect if all the sample is treated versus all the sample is untreated. The ATT corresponds to Average Treatment effect on the Treated, i.e. the marginal effect if the treated patients (\code{group = 1}) would have been untreated. The ATU corresponds to Average Treatment effect on the Untreated , i.e. the marginal effect if the untreated patients (\code{group = 0}) would have been treated. Simulation method for variance estimation has a shorter computing time than the boostrap method, but bootstrap is more accurate.
}

\value{
\item{effect}{A character string with the type of the marginal effect.}
\item{p0}{A table related to the average proportion of events in the unexposed/untreated sample: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{p1}{A table related to the average proportion of events in the exposed/treated sample: \code{estimate} is the estimated value, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{delta}{A table related to the difference between the average proportions of events in the exposed/treated sample minus in the unexposed/untreated sample: \code{estimate} is the estimated value, \code{std.error} is the corresponding standard error, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{logOR}{A table related to the logarithm of the average Odds Ratio (OR): \code{estimate} is the estimated value, \code{std.error} is the corresponding standard error, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval.}
\item{p.value}{The p-value of the bilateral test of the null hypothesis \code{p0} = \code{p1}, i.e. OR = 1.}
}

\references{
Chatton A, Le Borgne F, Gillaizeau F, Leyrat C, Barbin L, Laplaud D, Leger M, Giraudeau B, Foucher Y. G-computation and propensity score-based methods with binary outcome: a comparative study to evaluate their performances for causal inference. Manuscript submitted. 2019.
}

\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>

Arthur Chatton <Arthur.Chatton@etu.univ-nantes.fr>
}

\examples{
#data simulation
#treatment = 1 if the patients have been the exposure or treatment of interest and 0 otherwise
treatment <- rbinom(600, 1, prob=0.5)

covariate <- rnorm(600, 0, 1)
covariate[treatment==1] <- rnorm(sum(treatment==1), 0.3, 1)

outcome <- rbinom(600, 1, prob=
exp(-2+0.26*treatment+0.7*covariate)/(1+exp(-2+0.26*treatment+0.7*covariate)))

tab <- data.frame(outcome, treatment, covariate)

#Raw effect of the treatment
glm.raw <- glm(outcome ~ treatment, data=tab, family = binomial(link=logit))
summary(glm.raw)

#Conditional effect of the treatment
glm.multi <- glm(outcome ~ treatment + covariate, data=tab, family = binomial(link=logit))
summary(glm.multi)

#Marginal effects of the treatment (ATE)
gc.ate <- gc.logistic(glm.obj=glm.multi, data=tab, group="treatment", effect="ATE",
 var.method="simulations", iterations=1000, n.cluster=1)

#Sum-up of the 3 ORs
data.frame( raw=exp(glm.raw$coefficients[2]),
conditional=exp(glm.multi$coefficients[2]),
marginal.ate=exp(gc.ate$logOR[,1]) )
}


\keyword{G-computation}
