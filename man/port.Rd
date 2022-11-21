\name{port}
\alias{port}
\title{
POsitivity-Regression Tree (PoRT) Algorithm to Identify Positivity Violations.
}
\description{
This function allows to identify potential posivity violations by using the PoRT algorithm.
}
\usage{
port(group, cov.quanti, cov.quali, data, alpha, beta, gamma, pruning,
  minbucket, minsplit, maxdepth)
}

\arguments{
  \item{group}{A character string with the name of the exposure in \code{data}: 0 for the untreated/unexposed patients and 1 for the treated/exposed patients.}
  \item{cov.quanti}{A character string with the names of the quantitative predictors in \code{data}.}
  \item{cov.quali}{A character string with the names of the qualitative predictors in \code{data}.}
  \item{data}{A data frame in which to look for the variables related to the treatment/exposure and the predictors.}
  \item{alpha}{The minimal proportion of the whole sample size to consider a problematic subgroup. The default value is 0.05.}
 \item{beta}{The exposed or unexposed proportion under which one can consider a positivity violation. The default value is 0.05.}
 \item{gamma}{The maximum number of predictors used to define the subgroup.  The default value is 2. See 'Details'.} 
 \item{pruning}{If \code{TRUE}, provide only the violations contained between two values for quantitative predictors. The default value is \code{FALSE}.}
 \item{minbucket}{An \code{rpart} parameter: minimum number of observations in any leaf. The default value is 6.}
 \item{minsplit}{An \code{rpart} parameter: minimum number of observations that must exist in a node in order for a
split to be attempted. If only one of \code{minbucket} or \code{minsplit} is specified, the code either sets \code{minsplit} to
\code{minbucket*3} or \code{minbucket} to \code{minsplit/3}, as appropriate. The default value is 20.}
 \item{maxdepth}{An \code{rpart} parameter. Set the maximum depth of any node of the final tree, with the root node counted as depth 0. Values greater than 30 rpart will give nonsense results on 32-bit machines. The default value is 30.}
 }
 
\details{
In a first step, the PoRT algorithm estimates one tree for each predictor and memorises the leaves corresponding to problematic subgroups according to the hyperparameters \code{alpha} and \code{beta} (i.e., the subgroup must at least include \code{alpha*100} percent of the whole sample, and the exposure prevalence in the subgroup must be superior to \code{1-beta} or inferior to \code{beta}). If \code{gamma=1}, the algorithm stops. Otherwise, if at least one problematic subgroup is identified in this first step, the corresponding predictor(s) is(are) not considered in the second step, which estimates one tree for all possible couples of remaining predictors and memorizes the leaves corresponding to problematic subgroups. If \code{gamma=2}, the algorithm stops; otherwise, the third step consists of building one tree for all possible trios of remaining covariates not involved in the previously identified subgroups, etc.
}

\value{
The \code{port} function returns a characters string summarising all the subgroups identified as violating the positivity assumption, and provides for each of these subgroups the exposure prevalence, the subgroup size and the relative subgroup size (with respect to the sample size).  
}

\references{
Danelian et al. Identification of positivity violations' using regression trees: The PoRT algorithm. Manuscript submitted. 2022.
}

\author{
Arthur Chatton <Arthur.Chatton@univ-nantes.fr>
}

\examples{
data("dataDIVAT2")

# PoRT with default hyperparameters
port(group="ecd", cov.quanti="age", cov.quali=c("hla", "retransplant"),
data=dataDIVAT2)

# Illustration of the 'pruning' argument
port(group="ecd", cov.quanti="age", cov.quali=c("hla", "retransplant"),
    data=dataDIVAT2, beta=0.01)
    
port(group="ecd", cov.quanti="age", cov.quali=c("hla", "retransplant"),
    data=dataDIVAT2, beta=0.01, pruning=TRUE)
}

\keyword{Positivity}