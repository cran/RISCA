\name{hr.sl.time}
\alias{hr.sl.time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditionnal Effect for Censored Outcome with a Super Learner for the Outcome Model.
}
\description{
This function allows to estimate the conditional (subject-specific) effect of an exposure or a treatment for a censored times-to-event by using a Super Learner (SL).
}
\usage{
hr.sl.time(object, pro.time, iterations, n.cluster, cluster.type, conda.env)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned by the function \code{sl.time}. It shall be included the exposition/treatment of interest in the argument \code{group} of the \code{sl.time} function.}
  \item{pro.time}{The maximum time of follow-up to estimate of restricted mean survival time (RMST). The default value is \code{NULL}, which corresponds to \code{object$pro.time} if this argument is not NULL, or the time at which half of the sample is still at-risk otherwise.}
 \item{iterations}{The number of bootstrap resamples to estimate of the variances and the confidence intervals.}
  \item{n.cluster}{The number of cores to use, i.e., the upper limit for the number of child processes that run simultaneously (1 by default).}
 \item{cluster.type}{A character string with the type of parallelization. The default type is "PSOCK" (it calls makePSOCKcluster, faster on MacOS or Linux platforms). An alternative is "FORK" (it calls makeForkCluster, it does not work on Windows platforms).}
 \item{conda.env}{The name of the Conda environment to use if Python is used via reticulate. The default is \code{NULL}, which corresponds th the default environment by default located at ~/.virtualenvs.}
 }
 
\details{
The conditional effect is estimated by a conditional hazard ratio (HR), which corresponds to average of the individual HR on the studied population, i.e. the subject-specific effect if the individual is treated (\code{group} = 1) versus untreated (\code{group} = 0).
}

\value{
\item{logHR.conditional.values}{A vector with the estimated individual HRs.}
\item{logHR.conditional}{A table related to the logarithm of the conditional HR: \code{estimate} is the estimated value, \code{std.error} is the corresponding standard error, \code{ci.lower} and \code{ci.upper} represent the 95\% confidence interval, and \code{p.value} is the p-value of the bilateral test of the null hypothesis HR = 1.}
}

\references{
Foucher et al. Manuscript in progress. 2021. 
}

\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>
}

\examples{

data(dataDIVAT2)

#The outcome model base on a Super Learner from the first 150 individuals
sl1<-sl.time(  methods=c("aft.gamma", "ph.gompertz"),  metric="bs",
  data=dataDIVAT2[1:150,],  times="times", failures="failures",  group="ecd",
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant"),  cv=3)

#Set the prognostic time, for instance:
pro.time <- min(c(
  max(dataDIVAT2$times[dataDIVAT2$ecd==0]),
  max(dataDIVAT2$times[dataDIVAT2$ecd==1]) ))

#Conditional effect of the ECD by using a SL: use 1000 iterations instead of 2
hr.sl <- hr.sl.time(sl1, pro.time=pro.time, iterations=2, n.cluster=1)

#Conditional effect of the ECD by using a Cox model
cox.cdt <- coxph(Surv(times,failures) ~ ecd + age + retransplant,
  data=dataDIVAT2)

data.frame( rbind(
summary(cox.cdt)$coef[1,c(1,2,5)],
hr.sl$logHR.conditional[,c(1,2,5)] ))
}

\keyword{Super Learner}
