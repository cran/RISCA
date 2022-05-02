
gc.sl.binary<-function(outcome, group, cov.quanti=NULL, cov.quali=NULL, data, effect="ATE",
                tuneLength=20, cv=10, iterations=1000, n.cluster=1, cluster.type="PSOCK")
{

if(is.null(cov.quanti) & is.null(cov.quali)) {
	stop("This function is for causal inference with at least one covariate in \"cov.quanti\" or in \"cov.quali\" ") }

if(!is.character(outcome)){
	stop("The argument \"outcome\" must be a character string") }
	
if(length(outcome)!=1){
	stop("The argument \"outcome\" must be have a length equaled to 1") }
	
if(!is.character(group)){
	stop("The argument \"group\" must be a character string") }
	
if(length(group)!=1){
	stop("The argument \"group\" must be have a length equaled to 1") }

if(!is.data.frame(data)){
	stop("The argument \"data\" need to be a data.frame") }
	
if(!is.character(effect)){
	stop("The argument \"effect\" must be a vector of character strings") }

if(is.na(match(effect,c("ATT","ATU","ATE"))) | length(effect)!=1){
	stop("Incorrect modality specified in the argument \"effect\": only one option among \"ATE\", \"ATT\" or \"ATU\" is possible") }

if(length(cv)!=1) { stop("The argument \"cv\" must be have a length equaled to 1") }
	
if(!is.numeric(cv)) { stop("The argument \"cv\" must be a numeric") }
	
if(cv<1) {stop("The argument \"cv\" must higher than 1") }
	
if(length(n.cluster)!=1) { stop("The argument \"n.cluster\" must be have a length equaled to 1") }
	
if(!is.numeric(n.cluster)) { stop("The argument \"n.cluster\" must be a numeric") }
	
if(n.cluster<1) {stop("The argument \"n.cluster\" must higher than 1") }
	
if(length(tuneLength)!=1) { stop("The argument \"tuneLength\" must be have a length equaled to 1") }
	
if(!is.numeric(tuneLength)) { stop("The argument \"tuneLength\" must be a numeric") }
	
if(tuneLength<1) {stop("The argument \"tuneLength\" must higher than 1") }
	
mod <- unique(data[,group])

if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
	stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the argument \"group\" ") }

mod <- unique(data[,outcome])

if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
	stop("Two modalities encoded 1 (for event) and 0 (for absence of event) are required in the argument \"outcome\" ") }

if(!is.null(cov.quanti)){
	
if(!is.character(cov.quanti)){
	stop("The argument \"cov.quanti\" must be a vector of character strings") }
	
if(sum(apply(data.frame(data[,cov.quanti]), MARGIN=2, FUN="is.numeric"))!=length(cov.quanti)){
	stop("The argument \"cov.quanti\" must correspond to numeric variables in the dataframe") }  }
	
if(!is.null(cov.quali)){
	
if(!is.character(cov.quali)){
	stop("The argument \"cov.quali\" must be a vector of character strings") }
	
if(sum(apply(data.frame(data[,cov.quali]), MARGIN=2, FUN="is.numeric"))!=length(cov.quali)){
	stop("The argument \"cov.quali\" must correspond to numeric variables in the dataframe") }
	
if(sum(apply(data.frame(data[,cov.quali]), MARGIN=2, FUN="unique"))!=length(cov.quali)){
	stop("Two modalities encoded 1 and 0 are required in the variables declared in the argument \"cov.quali\" ") } }
	
.na <- is.na(apply(data[,c(outcome, group, cov.quali, cov.quanti)], MARGIN=1, FUN="sum"))

if(sum(.na)>0) {warning("Individuals with missing values will be removed from the analysis \n")}

data <- data[.na==FALSE,]

data.obs.caret <- data
data.obs.caret[data[,outcome]==1, outcome] <- "yes"
data.obs.caret[data[,outcome]==0, outcome] <- "no"
data.obs.caret[ , outcome] <- as.factor(data.obs.caret[ , outcome])

if(cv>=dim(data.obs.caret)[1]){ stop("The argument \"cv\" must be inferior than the number of subjects with no missing data") }

if(n.cluster==1) {allowParallel <- FALSE}  else { 
  allowParallel <- TRUE
  cl <- makeCluster(n.cluster, type=cluster.type)
  registerDoParallel(cl)
  clusterEvalQ(cl, {library(splines); library(SuperLearner)}) }

control <- trainControl(allowParallel = TRUE,  verboseIter = FALSE, classProbs = TRUE,
                        summaryFunction = twoClassSummary,  method = "cv", number = cv)

.f.nnet <- as.formula(paste(outcome, "~ .", sep = " "))

nnet <-  train(.f.nnet, data = data.obs.caret, method = 'nnet', tuneLength = tuneLength,
               metric = "ROC", trControl = control, trace = FALSE)


if(!is.null(cov.quanti) & !is.null(cov.quali)) {
.f.glm <- as.formula( paste(outcome, "~", group, "*(", paste("bs(", cov.quanti, ", df=3)",
            collapse = " + "), " + ", paste(cov.quali, collapse = " + "),  ")", collapse = " ") ) }

if(!is.null(cov.quanti) & is.null(cov.quali)) {
.f.glm <- as.formula( paste(outcome, "~", group, "*(", paste("bs(", cov.quanti, ", df=3)",
            collapse = " + "), ")", collapse = " ") ) }

if(is.null(cov.quanti) & !is.null(cov.quali)) {
.f.glm <- as.formula( paste(outcome, "~", group, "*(", paste(cov.quali, collapse = " + "),
            ")", collapse = " ") ) }

full <- glm( .f.glm, family = binomial(link = logit),  data = data)

.l <- length(full$coefficients)

elasticnet <-  train(.f.glm, data = data.obs.caret, method = 'glmnet', tuneLength = tuneLength,
    metric = "ROC", trControl = control, family = "binomial", penalty.factor = c(0, rep(1, .l-1)) )

lasso <- train(.f.glm, data = data.obs.caret, method = 'glmnet',
               tuneGrid = expand.grid(.alpha = 1, .lambda = unique(elasticnet$results$lambda)),
               metric = "ROC", trControl = control, family = "binomial",
               penalty.factor = c(0, rep(1, .l-1)) )

svmRadial <-  train(.f.nnet, data = data.obs.caret, method = 'svmRadialSigma', trace = FALSE,
               tuneLength = tuneLength, metric = "ROC", trControl = control)

N <- dim(data)[1]

if(!is.null(cov.quanti) & !is.null(cov.quali)) {
.f.caret <- as.formula( paste("~ -1 +", group, "*(", paste("bs(", cov.quanti, ", df=3)", collapse = " + "), " + ", paste(cov.quali, collapse = " + "),  ")", collapse = " ") ) }

if(!is.null(cov.quanti) & is.null(cov.quali)) {
.f.caret <- as.formula( paste("~ -1 +", group, "*(", paste("bs(", cov.quanti, ", df=3)", collapse = " + "), ")", collapse = " ") ) }

if(is.null(cov.quanti) & !is.null(cov.quali)) {
.f.caret <- as.formula( paste("~ -1 +", group, "*(", paste(cov.quali, collapse = " + "),  ")", collapse = " ") ) }

if(n.cluster==1){
	
	SL.elasticnet.caret <- function (Y, X, newX, family,  ...)
	{
	
	X <- model.matrix(.f.caret, X)
  
	newX <- model.matrix(.f.caret,  newX)
  
	fitElastic <- glmnet(x = X, y = Y, family = "binomial", alpha = elasticnet$bestTune$alpha, lambda = elasticnet$bestTune$lambda, penalty.factor = c(0, rep(1, .l-1)) )
  
	pred <- predict(fitElastic, newx = newX, type = "response")
	fit <- list(object = fitElastic)
	class(fit) <- "SL.elasticnet.caret"
	out <- list(pred = pred, fit = fit)
	return(out)
	}

	SL.lasso.caret <- function (Y, X, newX, family,  ...)
	{
	
	X <-  model.matrix(.f.caret,  X)
  
	newX <-  model.matrix(.f.caret,  newX)
  
	fitLasso <- glmnet(x = X, y = Y, family = "binomial", alpha = 1,  lambda = lasso$bestTune$lambda,  penalty.factor = c(0, rep(1, .l-1)))
  
	pred <- predict(fitLasso, newx = newX, type = "response")
	fit <- list(object = fitLasso)
	class(fit) <- "SL.lasso.caret"
	out <- list(pred = pred, fit = fit)
	return(out)
	}

	SL.nnet.caret <- function (Y, X, newX, family, ...)
	{
	fit.nnet <- nnet(x = X,  y = Y, size = nnet$bestTune$size, decay = nnet$bestTune$decay, trace = FALSE)
  
	pred <- predict(fit.nnet, newdata = newX, type = "raw")
	fit <- list(object = fit.nnet)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.nnet.caret")
	return(out)
	}

	SL.ksvm.caret <- function (Y, X, newX, family, type = "C-svc", kernel = "rbfdot", kpar = "automatic", nu = 0.2, epsilon = 0.1, cross = 0, prob.model = family$family == "binomial", class.weights = NULL, cache = 40, tol = 0.001, shrinking = T, ...)
	{
	if (!is.matrix(X)) {X <- model.matrix( ~ ., data = X); X <- X[,-1]}
  
	Y <- as.factor(Y)
  
	predict_type <- "probabilities"
  
	model <- ksvm(X, Y, type = "C-svc", kernel = "rbfdot",  kpar = list(sigma = svmRadial$bestTune$sigma), C = svmRadial$bestTune$C, nu = nu, epsilon = epsilon, prob.model = prob.model, class.weights = class.weights)
  
	if (!is.matrix(newX)) {newX <-  model.matrix( ~ ., data = newX); newX <- newX[,-1, drop = FALSE]}
  
	pred <- kernlab::predict(model, newX, predict_type)[, 2]
	fit <- list(object = model, family = family)
	out <- list(pred = pred, fit = fit)
	class(out$fit) = "SL.ksvm.caret"
	return(out)
	}
	
p1 <- p0 <- logOR <- delta <- rep(-99, iterations)

if(effect=="ATE"){
		
for(i in 1:iterations){
			  
  id <- sample(1:N, size = N, replace = TRUE)
  learn <- data[id,]
  valid <- data[-sort(unique(id)),]
  
  learn.caret <- learn
  learn.caret[learn[,outcome] == 1, outcome] <- "yes"
  learn.caret[learn[,outcome] == 0, outcome] <- "no"
  learn.caret[ , outcome] <- as.factor(learn.caret[ , outcome])
  
  valid.caret <- valid
  valid.caret[valid[,outcome] == 1, outcome] <- "yes"
  valid.caret[valid[,outcome] == 0, outcome] <- "no"
  valid.caret[ , outcome] <- as.factor(valid.caret[ , outcome])
  
  valid.caret1 <- valid.caret0 <- valid.caret;
  valid.caret1[ , group] <- 1; valid.caret0[ , group] <- 0
  
  y.sl <-  as.numeric(learn.caret[ , outcome]) - 1
  x.sl <-  learn.caret[ , c(group, cov.quanti, cov.quali)]
  
  N.valid <- length(valid.caret[ , outcome])
  
  y.sl.valid <-  as.numeric(valid.caret[ , outcome]) - 1
  x.sl.valid <-  valid.caret[ , c(group, cov.quanti, cov.quali)]
  x.sl.valid <- data.frame(rbind(x.sl.valid, x.sl.valid))
  x.sl.valid[1:N.valid, group] <- 1
  x.sl.valid[(N.valid + 1):(2 * N.valid), group] <- 0
  
  sl <- SuperLearner(Y = y.sl, X = x.sl, newX = x.sl.valid, family = binomial(), SL.library = list("SL.nnet.caret", "SL.elasticnet.caret", "SL.lasso.caret", "SL.ksvm.caret"),  "method.AUC", cvControl = list(V = cv) )
  
  .pred <- predict(sl, onlySL = TRUE)$pred
  
  p1[i] <- mean(.pred[1:N.valid])
  p0[i] <- mean(.pred[(N.valid + 1):(2 * N.valid)])
  logOR[i] <- log((p1[i] / (1 - p1[i])) / (p0[i] / (1 - p0[i])))
  delta[i] <- p1[i] - p0[i]
}
}

if(effect=="ATT"){

for(i in 1:iterations){
			  
  id <- sample(1:N, size = N, replace = TRUE)
  learn <- data[id,]
  valid <- data[-sort(unique(id)),]
  
  learn.caret <- learn
  learn.caret[learn[,outcome] == 1, outcome] <- "yes"
  learn.caret[learn[,outcome] == 0, outcome] <- "no"
  learn.caret[ , outcome] <- as.factor(learn.caret[ , outcome])
  
  valid.caret <- valid[valid[,group]==1,]
  valid.caret[valid[valid[,group]==1,outcome] == 1, outcome] <- "yes"
  valid.caret[valid[valid[,group]==1,outcome] == 0, outcome] <- "no"
  valid.caret[ , outcome] <- as.factor(valid.caret[ , outcome])
  
  valid.caret1 <- valid.caret0 <- valid.caret;
  valid.caret1[ , group] <- 1; valid.caret0[ , group] <- 0
  
  y.sl <-  as.numeric(learn.caret[ , outcome]) - 1
  x.sl <-  learn.caret[ , c(group, cov.quanti, cov.quali)]
  
  N.valid <- length(valid.caret[ , outcome])
  
  y.sl.valid <-  as.numeric(valid.caret[ , outcome]) - 1
  x.sl.valid <-  valid.caret[ , c(group, cov.quanti, cov.quali)]
  x.sl.valid <- data.frame(rbind(x.sl.valid, x.sl.valid))
  x.sl.valid[1:N.valid, group] <- 1
  x.sl.valid[(N.valid + 1):(2 * N.valid), group] <- 0
  
  sl <- SuperLearner(Y = y.sl, X = x.sl, newX = x.sl.valid, family = binomial(), SL.library = list("SL.nnet.caret", "SL.elasticnet.caret", "SL.lasso.caret", "SL.ksvm.caret"),  "method.AUC", cvControl = list(V = cv) )
  
  .pred <- predict(sl, onlySL = TRUE)$pred
  
  p1[i] <- mean(.pred[1:N.valid])
  p0[i] <- mean(.pred[(N.valid + 1):(2 * N.valid)])
  logOR[i] <- log((p1[i] / (1 - p1[i])) / (p0[i] / (1 - p0[i])))
  delta[i] <- p1[i] - p0[i]
}

}


if(effect=="ATU"){

for(i in 1:iterations){
			  
  id <- sample(1:N, size = N, replace = TRUE)
  learn <- data[id,]
  valid <- data[-sort(unique(id)),]
  
  learn.caret <- learn
  learn.caret[learn[,outcome] == 1, outcome] <- "yes"
  learn.caret[learn[,outcome] == 0, outcome] <- "no"
  learn.caret[ , outcome] <- as.factor(learn.caret[ , outcome])
  
  valid.caret <- valid[valid[,group]==0,]
  valid.caret[valid[valid[,group]==0,outcome] == 1, outcome] <- "yes"
  valid.caret[valid[valid[,group]==0,outcome] == 0, outcome] <- "no"
  valid.caret[ , outcome] <- as.factor(valid.caret[ , outcome])
  
  valid.caret1 <- valid.caret0 <- valid.caret;
  valid.caret1[ , group] <- 1; valid.caret0[ , group] <- 0
  
  y.sl <-  as.numeric(learn.caret[ , outcome]) - 1
  x.sl <-  learn.caret[ , c(group, cov.quanti, cov.quali)]
  
  N.valid <- length(valid.caret[ , outcome])
  
  y.sl.valid <-  as.numeric(valid.caret[ , outcome]) - 1
  x.sl.valid <-  valid.caret[ , c(group, cov.quanti, cov.quali)]
  x.sl.valid <- data.frame(rbind(x.sl.valid, x.sl.valid))
  x.sl.valid[1:N.valid, group] <- 1
  x.sl.valid[(N.valid + 1):(2 * N.valid), group] <- 0
  
  sl <- SuperLearner(Y = y.sl, X = x.sl, newX = x.sl.valid, family = binomial(), SL.library = list("SL.nnet.caret", "SL.elasticnet.caret", "SL.lasso.caret", "SL.ksvm.caret"),  "method.AUC", cvControl = list(V = cv) )
  
  .pred <- predict(sl, onlySL = TRUE)$pred
  
  p1[i] <- mean(.pred[1:N.valid])
  p0[i] <- mean(.pred[(N.valid + 1):(2 * N.valid)])
  logOR[i] <- log((p1[i] / (1 - p1[i])) / (p0[i] / (1 - p0[i])))
  delta[i] <- p1[i] - p0[i]
}

}

}


if(n.cluster>1){
	
if(effect=="ATE"){
		
bsim <- function() {
  id <- sample(1:N, size = N, replace = TRUE)
  learn <- data[id,]
  valid <- data[-sort(unique(id)),]
  
  learn.caret <- learn
  learn.caret[learn[,outcome] == 1, outcome] <- "yes"
  learn.caret[learn[,outcome] == 0, outcome] <- "no"
  learn.caret[ , outcome] <- as.factor(learn.caret[ , outcome])
  
  valid.caret <- valid
  valid.caret[valid[,outcome] == 1, outcome] <- "yes"
  valid.caret[valid[,outcome] == 0, outcome] <- "no"
  valid.caret[ , outcome] <- as.factor(valid.caret[ , outcome])
  
  valid.caret1 <- valid.caret0 <- valid.caret;
  valid.caret1[ , group] <- 1; valid.caret0[ , group] <- 0
  
  y.sl <-  as.numeric(learn.caret[ , outcome]) - 1
  x.sl <-  learn.caret[ , c(group, cov.quanti, cov.quali)]
  
  N.valid <- length(valid.caret[ , outcome])
  
  y.sl.valid <-  as.numeric(valid.caret[ , outcome]) - 1
  x.sl.valid <-  valid.caret[ , c(group, cov.quanti, cov.quali)]
  x.sl.valid <- data.frame(rbind(x.sl.valid, x.sl.valid))
  x.sl.valid[1:N.valid, group] <- 1
  x.sl.valid[(N.valid + 1):(2 * N.valid), group] <- 0
  
	SL.elasticnet.caret <- function (Y, X, newX, family,  ...)
	{
	X <- model.matrix(.f.caret, X)
  
	newX <- model.matrix(.f.caret,  newX)
  
	fitElastic <- glmnet(x = X, y = Y, family = "binomial", alpha = elasticnet$bestTune$alpha, lambda = elasticnet$bestTune$lambda, penalty.factor = c(0, rep(1, .l-1)) )
  
	pred <- predict(fitElastic, newx = newX, type = "response")
	fit <- list(object = fitElastic)
	class(fit) <- "SL.elasticnet.caret"
	out <- list(pred = pred, fit = fit)
	return(out)
	}

	SL.lasso.caret <- function (Y, X, newX, family,  ...)
	{
	X <-  model.matrix(.f.caret,  X)
  
	newX <-  model.matrix(.f.caret,  newX)
  
	fitLasso <- glmnet(x = X, y = Y, family = "binomial", alpha = 1,  lambda = lasso$bestTune$lambda,  penalty.factor = c(0, rep(1, .l-1)))
  
	pred <- predict(fitLasso, newx = newX, type = "response")
	fit <- list(object = fitLasso)
	class(fit) <- "SL.lasso.caret"
	out <- list(pred = pred, fit = fit)
	return(out)
	}

	SL.nnet.caret <- function (Y, X, newX, family, ...)
	{
	fit.nnet <- nnet::nnet(x = X,  y = Y, size = nnet$bestTune$size, decay = nnet$bestTune$decay, trace = FALSE)
  
	pred <- predict(fit.nnet, newdata = newX, type = "raw")
	fit <- list(object = fit.nnet)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.nnet.caret")
	return(out)
	}

	SL.ksvm.caret <- function (Y, X, newX, family, type = "C-svc", kernel = "rbfdot", kpar = "automatic", nu = 0.2, epsilon = 0.1, cross = 0, prob.model = family$family == "binomial", class.weights = NULL, cache = 40, tol = 0.001, shrinking = T, ...)
	{
	if (!is.matrix(X)) {X <- model.matrix( ~ ., data = X); X <- X[,-1]}
  
	Y <- as.factor(Y)
  
	predict_type <- "probabilities"
  
	model <- kernlab::ksvm(X, Y, type = "C-svc", kernel = "rbfdot",  kpar = list(sigma = svmRadial$bestTune$sigma), C = svmRadial$bestTune$C, nu = nu, epsilon = epsilon, prob.model = prob.model, class.weights = class.weights)
  
	if (!is.matrix(newX)) {newX <-  model.matrix( ~ ., data = newX); newX <- newX[,-1, drop = FALSE]}
  
	pred <- kernlab::predict(model, newX, predict_type)[, 2]
	fit <- list(object = model, family = family)
	out <- list(pred = pred, fit = fit)
	class(out$fit) = "SL.ksvm.caret"
	return(out)
	}

  sl <- SuperLearner(Y = y.sl, X = x.sl, newX = x.sl.valid, family = binomial(), SL.library = list("SL.nnet.caret", "SL.elasticnet.caret", "SL.lasso.caret", "SL.ksvm.caret"),  "method.AUC", cvControl = list(V = cv) )
  
  .pred <- predict(sl, onlySL = TRUE)$pred
  
  p1 <- mean(.pred[1:N.valid])
  p0 <- mean(.pred[(N.valid + 1):(2 * N.valid)])
  logOR <- log((p1 / (1 - p1)) / (p0 / (1 - p0)))
  delta <- p1 - p0
  
  return(c(p1=p1, p0=p0, logOR=logOR, delta=delta))
}
}

if(effect=="ATT"){
			
bsim <- function(){
  id <- sample(1:N, size = N, replace = TRUE)
  learn <- data[id,]
  valid <- data[-sort(unique(id)),]
  
  learn.caret <- learn
  learn.caret[learn[,outcome] == 1, outcome] <- "yes"
  learn.caret[learn[,outcome] == 0, outcome] <- "no"
  learn.caret[ , outcome] <- as.factor(learn.caret[ , outcome])
  
  valid.caret <- valid[valid[,group]==1,]
  valid.caret[valid[valid[,group]==1,outcome] == 1, outcome] <- "yes"
  valid.caret[valid[valid[,group]==1,outcome] == 0, outcome] <- "no"
  valid.caret[ , outcome] <- as.factor(valid.caret[ , outcome])
  
  valid.caret1 <- valid.caret0 <- valid.caret;
  valid.caret1[ , group] <- 1; valid.caret0[ , group] <- 0
  
  y.sl <-  as.numeric(learn.caret[ , outcome]) - 1
  x.sl <-  learn.caret[ , c(group, cov.quanti, cov.quali)]
  
  N.valid <- length(valid.caret[ , outcome])
  
  y.sl.valid <-  as.numeric(valid.caret[ , outcome]) - 1
  x.sl.valid <-  valid.caret[ , c(group, cov.quanti, cov.quali)]
  x.sl.valid <- data.frame(rbind(x.sl.valid, x.sl.valid))
  x.sl.valid[1:N.valid, group] <- 1
  x.sl.valid[(N.valid + 1):(2 * N.valid), group] <- 0
  
	SL.elasticnet.caret <- function (Y, X, newX, family,  ...)
	{
	X <- model.matrix(.f.caret, X)
  
	newX <- model.matrix(.f.caret,  newX)
  
	fitElastic <- glmnet(x = X, y = Y, family = "binomial", alpha = elasticnet$bestTune$alpha, lambda = elasticnet$bestTune$lambda, penalty.factor = c(0, rep(1, .l-1)) )
  
	pred <- predict(fitElastic, newx = newX, type = "response")
	fit <- list(object = fitElastic)
	class(fit) <- "SL.elasticnet.caret"
	out <- list(pred = pred, fit = fit)
	return(out)
	}

	SL.lasso.caret <- function (Y, X, newX, family,  ...)
	{

	X <-  model.matrix(.f.caret,  X)
  
	newX <-  model.matrix(.f.caret,  newX)
  
	fitLasso <- glmnet(x = X, y = Y, family = "binomial", alpha = 1,  lambda = lasso$bestTune$lambda,  penalty.factor = c(0, rep(1, .l-1)))
  
	pred <- predict(fitLasso, newx = newX, type = "response")
	fit <- list(object = fitLasso)
	class(fit) <- "SL.lasso.caret"
	out <- list(pred = pred, fit = fit)
	return(out)
	}

	SL.nnet.caret <- function (Y, X, newX, family, ...)
	{
	fit.nnet <- nnet::nnet(x = X,  y = Y, size = nnet$bestTune$size, decay = nnet$bestTune$decay, trace = FALSE)
  
	pred <- predict(fit.nnet, newdata = newX, type = "raw")
	fit <- list(object = fit.nnet)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.nnet.caret")
	return(out)
	}

	SL.ksvm.caret <- function (Y, X, newX, family, type = "C-svc", kernel = "rbfdot", kpar = "automatic", nu = 0.2, epsilon = 0.1, cross = 0, prob.model = family$family == "binomial", class.weights = NULL, cache = 40, tol = 0.001, shrinking = T, ...)
	{
	if (!is.matrix(X)) {X <- model.matrix( ~ ., data = X); X <- X[,-1]}
  
	Y <- as.factor(Y)
  
	predict_type <- "probabilities"
  
	model <- kernlab::ksvm(X, Y, type = "C-svc", kernel = "rbfdot",  kpar = list(sigma = svmRadial$bestTune$sigma), C = svmRadial$bestTune$C, nu = nu, epsilon = epsilon, prob.model = prob.model, class.weights = class.weights)
  
	if (!is.matrix(newX)) {newX <-  model.matrix( ~ ., data = newX); newX <- newX[,-1, drop = FALSE]}
  
	pred <- kernlab::predict(model, newX, predict_type)[, 2]
	fit <- list(object = model, family = family)
	out <- list(pred = pred, fit = fit)
	class(out$fit) = "SL.ksvm.caret"
	return(out)
	}

  sl <- SuperLearner(Y = y.sl, X = x.sl, newX = x.sl.valid, family = binomial(), SL.library = list("SL.nnet.caret", "SL.elasticnet.caret", "SL.lasso.caret", "SL.ksvm.caret"),  "method.AUC", cvControl = list(V = cv) )
  
  .pred <- predict(sl, onlySL = TRUE)$pred
  
  p1 <- mean(.pred[1:N.valid])
  p0 <- mean(.pred[(N.valid + 1):(2 * N.valid)])
  logOR <- log((p1 / (1 - p1)) / (p0 / (1 - p0)))
  delta <- p1 - p0
  
  return(c(p1=p1, p0=p0, logOR=logOR, delta=delta))
}
}

if(effect=="ATU"){			

bsim <- function(){
  id <- sample(1:N, size = N, replace = TRUE)
  learn <- data[id,]
  valid <- data[-sort(unique(id)),]
  
  learn.caret <- learn
  learn.caret[learn[,outcome] == 1, outcome] <- "yes"
  learn.caret[learn[,outcome] == 0, outcome] <- "no"
  learn.caret[ , outcome] <- as.factor(learn.caret[ , outcome])
  
  valid.caret <- valid[valid[,group]==0,]
  valid.caret[valid[valid[,group]==0,outcome] == 1, outcome] <- "yes"
  valid.caret[valid[valid[,group]==0,outcome] == 0, outcome] <- "no"
  valid.caret[ , outcome] <- as.factor(valid.caret[ , outcome])
  
  valid.caret1 <- valid.caret0 <- valid.caret;
  valid.caret1[ , group] <- 1; valid.caret0[ , group] <- 0
  
  y.sl <-  as.numeric(learn.caret[ , outcome]) - 1
  x.sl <-  learn.caret[ , c(group, cov.quanti, cov.quali)]
  
  N.valid <- length(valid.caret[ , outcome])
  
  y.sl.valid <-  as.numeric(valid.caret[ , outcome]) - 1
  x.sl.valid <-  valid.caret[ , c(group, cov.quanti, cov.quali)]
  x.sl.valid <- data.frame(rbind(x.sl.valid, x.sl.valid))
  x.sl.valid[1:N.valid, group] <- 1
  x.sl.valid[(N.valid + 1):(2 * N.valid), group] <- 0
  
	SL.elasticnet.caret <- function (Y, X, newX, family,  ...)
	{

	X <- model.matrix(.f.caret, X)
  
	newX <- model.matrix(.f.caret,  newX)
  
	fitElastic <- glmnet(x = X, y = Y, family = "binomial", alpha = elasticnet$bestTune$alpha, lambda = elasticnet$bestTune$lambda, penalty.factor = c(0, rep(1, .l-1)) )
  
	pred <- predict(fitElastic, newx = newX, type = "response")
	fit <- list(object = fitElastic)
	class(fit) <- "SL.elasticnet.caret"
	out <- list(pred = pred, fit = fit)
	return(out)
	}

	SL.lasso.caret <- function (Y, X, newX, family,  ...)
	{
	X <-  model.matrix(.f.caret,  X)
  
	newX <-  model.matrix(.f.caret,  newX)
  
	fitLasso <- glmnet(x = X, y = Y, family = "binomial", alpha = 1,  lambda = lasso$bestTune$lambda,  penalty.factor = c(0, rep(1, .l-1)))
  
	pred <- predict(fitLasso, newx = newX, type = "response")
	fit <- list(object = fitLasso)
	class(fit) <- "SL.lasso.caret"
	out <- list(pred = pred, fit = fit)
	return(out)
	}

	SL.nnet.caret <- function (Y, X, newX, family, ...)
	{
	fit.nnet <- nnet::nnet(x = X,  y = Y, size = nnet$bestTune$size, decay = nnet$bestTune$decay, trace = FALSE)
  
	pred <- predict(fit.nnet, newdata = newX, type = "raw")
	fit <- list(object = fit.nnet)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.nnet.caret")
	return(out)
	}

	SL.ksvm.caret <- function (Y, X, newX, family, type = "C-svc", kernel = "rbfdot", kpar = "automatic", nu = 0.2, epsilon = 0.1, cross = 0, prob.model = family$family == "binomial", class.weights = NULL, cache = 40, tol = 0.001, shrinking = T, ...)
	{
	if (!is.matrix(X)) {X <- model.matrix( ~ ., data = X); X <- X[,-1]}
  
	Y <- as.factor(Y)
  
	predict_type <- "probabilities"
  
	model <- kernlab::ksvm(X, Y, type = "C-svc", kernel = "rbfdot",  kpar = list(sigma = svmRadial$bestTune$sigma), C = svmRadial$bestTune$C, nu = nu, epsilon = epsilon, prob.model = prob.model, class.weights = class.weights)
  
	if (!is.matrix(newX)) {newX <-  model.matrix( ~ ., data = newX); newX <- newX[,-1, drop = FALSE]}
  
	pred <- kernlab::predict(model, newX, predict_type)[, 2]
	fit <- list(object = model, family = family)
	out <- list(pred = pred, fit = fit)
	class(out$fit) = "SL.ksvm.caret"
	return(out)
	}

  sl <- SuperLearner(Y = y.sl, X = x.sl, newX = x.sl.valid, family = binomial(), SL.library = list("SL.nnet.caret", "SL.elasticnet.caret", "SL.lasso.caret", "SL.ksvm.caret"),  "method.AUC", cvControl = list(V = cv) )
  
  .pred <- predict(sl, onlySL = TRUE)$pred
  
  p1 <- mean(.pred[1:N.valid])
  p0 <- mean(.pred[(N.valid + 1):(2 * N.valid)])
  logOR <- log((p1 / (1 - p1)) / (p0 / (1 - p0)))
  delta <- p1 - p0
  
  return(c(p1=p1, p0=p0, logOR=logOR, delta=delta))
}
}
		
res <- NULL
res <- foreach(i = 1:iterations, .combine=rbind, .inorder=TRUE) %dopar% {set.seed(i); bsim()}

registerDoSEQ()
stopCluster(cl)

p1 <- as.numeric(res[,"p1"])
p0 <- as.numeric(res[,"p0"])
logOR <- as.numeric(res[,"logOR"])
delta <- as.numeric(res[,"delta"])

}

se.delta <- sd(delta, na.rm=TRUE)
se.logOR <- sd(logOR, na.rm=TRUE)
se.p0 <- sd(p0, na.rm=TRUE)
se.p1 <- sd(p1, na.rm=TRUE)

mean.delta <- mean(delta, na.rm=TRUE)
mean.logOR <- mean(logOR, na.rm=TRUE)
mean.p0 <- mean(p0, na.rm=TRUE)
mean.p1 <- mean(p1, na.rm=TRUE)

pv <- function(m, x){
  ztest <- m/sd(x)
  return(ifelse(ztest<0,2*pnorm(ztest),2*(1-pnorm(ztest))))
}

p.value.OR <- pv(m=mean.logOR, x=logOR)
  

if(p.value.OR==0){ p.value.OR <- "<0.001" }

ci.low.logOR <- mean.logOR - qnorm(0.975, 0, 1)*se.logOR
ci.upp.logOR <- mean.logOR + qnorm(0.975, 0, 1)*se.logOR

ci.low.delta <- mean.delta - qnorm(0.975, 0, 1)*se.delta
ci.upp.delta <- mean.delta + qnorm(0.975, 0, 1)*se.delta

ci.low.p0 <- mean.p0 - qnorm(0.975, 0, 1)*se.p0
ci.upp.p0 <- mean.p0 + qnorm(0.975, 0, 1)*se.p0

ci.low.p1 <- mean.p1 - qnorm(0.975, 0, 1)*se.p1
ci.upp.p1 <- mean.p1 + qnorm(0.975, 0, 1)*se.p1

return( list(
missing=sum(.na),
effect=effect,
p0=data.frame(estimate=mean.p0, ci.lower=ci.low.p0, ci.upper=ci.upp.p0, row.names = NULL),
p1=data.frame(estimate=mean.p1, ci.lower=ci.low.p1, ci.upper=ci.upp.p1, row.names = NULL),
delta=data.frame(estimate=mean.delta, std.error = se.delta, ci.lower=ci.low.delta, ci.upper=ci.upp.delta, row.names = NULL),
logOR=data.frame(estimate=mean.logOR, std.error = se.logOR, ci.lower=ci.low.logOR, ci.upper=ci.upp.logOR, row.names = NULL),
p.value=p.value.OR) )

}

