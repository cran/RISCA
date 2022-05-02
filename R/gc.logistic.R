gc.logistic<-function(glm.obj, data, group, effect="ATE", var.method="simulations", iterations=1000,
                      n.cluster=1, cluster.type="PSOCK")
{

if(glm.obj$method != "glm.fit"){  
	stop("The argument \"glm.obj\" need to be a glm object") }
	
if(!is.data.frame(data) & !is.matrix(data)){
	stop("The argument \"data\" need to be a data.frame or a matrix") }
	
if(!is.character(group) & !is.numeric(group)){
	stop("The argument \"group\" need to be scalar or a character string") }
  
  mod <- paste(sort(unique(data[,group])), collapse = ",")
  if(!(mod %in% c(paste(c(0,1) , collapse = ",")))){
    stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and
1 (for treated/exposed patients) are required in the argument \"group\" ") }
  
if(is.na(match(effect,c("ATT","ATU","ATE"))) | length(effect)!=1){
	stop("Incorrect modality specified in the argument \"effect\": only one option among \"ATE\", \"ATT\" or \"ATU\" is possible") }
	
if(is.na(match(var.method,c("simulations","bootstrap"))) | length(var.method)!=1){
	stop("Incorrect modality specified in the argument \"var.method\": only one option among \"simulations\" and \"bootstrap\" is possible") }

if(length(grep("$", names(glm.obj$model), fixed = TRUE)) > 0 | length(grep("[", names(glm.obj$model), fixed = TRUE)) > 0){
	stop("Incorrect formula specified in the argument \"glm.obj\": don't use the syntax data$var or data[,var]") }
	
if(length(names(glm.obj$model)) < 3){
	stop("Incorrect formula specified in the argument \"glm.obj\": 2 covariables, including exposure, are requiered") }

	
if(is.na(match(group,names(glm.obj$model))) & is.na(match(names(data)[group],names(glm.obj$model)))){
	stop("Incorrect formula specified in the argument \"glm.obj\": exposure is requiered") } 

if(effect=="ATE"){ .ttt <- which(data[,group] == 0 | data[,group]==1) }

if(effect=="ATT"){ .ttt <- which(data[,group] == 1) }

if(effect=="ATU"){  .ttt <- which(data[,group] == 0) }

d1 <- d0 <- data[.ttt,]
d1[,group] <- 1
d0[,group] <- 0

p1 <- mean(predict(glm.obj, newdata=d1, type="response", na.action=na.omit))
p0 <- mean(predict(glm.obj, newdata=d0, type="response", na.action=na.omit))
logOR <- log((p1/(1-p1))/(p0/(1-p0)))

logOR.s <- p1.s <- p0.s <- rep(-99, iterations)

if(var.method=="simulations")
{
simul <- mvrnorm(n = iterations, mu=as.numeric(glm.obj$coef), Sigma=summary(glm.obj)$cov.unscaled)

	if(n.cluster==1){
		for(i in 1:iterations)	{
		glm.obj$coefficients <- simul[i,]
		p1.s[i] <- mean(predict(glm.obj, newdata=d1, type="response", na.action=na.omit))
		p0.s[i] <- mean(predict(glm.obj, newdata=d0, type="response", na.action=na.omit))
		logOR.s[i] <- log((p1.s[i]/(1-p1.s[i]))/(p0.s[i]/(1-p0.s[i])))
		}
	}
	
	if(n.cluster>1){
		ssimul <- function(i, glm.obj){
		
			glm.obj$coefficients <- simul[i,]
			p1.s <- mean(predict(glm.obj, newdata=d1, type="response", na.action=na.omit))
			p0.s <- mean(predict(glm.obj, newdata=d0, type="response", na.action=na.omit))
			
			return(c(logOR.s = log((p1.s/(1-p1.s))/(p0.s/(1-p0.s))),  p1.s = p1.s, p0.s = p0.s))
		}
		
		cl <- makeCluster(n.cluster, type=cluster.type) 
		registerDoParallel(cl)
    clusterEvalQ(cl, {library(splines)})
		res <- NULL
		res <- foreach(i = 1:iterations, .combine=rbind, .inorder=TRUE) %dopar% {ssimul(i, glm.obj)}
		registerDoSEQ()
		stopCluster(cl) 

		logOR.s <- as.numeric(res[,"logOR.s"])
		p1.s <- as.numeric(res[,"p1.s"])
		p0.s <- as.numeric(res[,"p0.s"])
	}
}
  
  
if(var.method=="bootstrap"){

	if(n.cluster==1){
		if(effect=="ATE"){
			for(i in 1:iterations){
			  dboot <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]
			  glm.s <- glm(glm.obj$form, data=dboot, family=binomial(link=logit))
			  .d1 <- .d0 <- dboot
			  .d1[,group] <- 1; .d0[,group] <- 0
			  p1.s[i] <- mean(predict(glm.s, newdata=.d1, type="response", na.action=na.omit))
			  p0.s[i] <- mean(predict(glm.s, newdata=.d0, type="response", na.action=na.omit))
			  logOR.s[i] <- log((p1.s[i]/(1-p1.s[i]))/(p0.s[i]/(1-p0.s[i])))
			}
		}else{
			if(effect=="ATT"){
				for(i in 1:iterations){
				  dboot <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]
				  glm.s <- glm(glm.obj$form, data=dboot, family=binomial(link=logit))
				  .d1 <- .d0 <- dboot[dboot[,group] == 1,]
				  .d1[,group] <- 1; .d0[,group] <- 0
				  p1.s[i] <- mean(predict(glm.s, newdata=.d1, type="response", na.action=na.omit))
				  p0.s[i] <- mean(predict(glm.s, newdata=.d0, type="response", na.action=na.omit))
				  logOR.s[i] <- log((p1.s[i]/(1-p1.s[i]))/(p0.s[i]/(1-p0.s[i])))
				}
			}else{
				for(i in 1:iterations){
				  dboot <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]
				  glm.s <- glm(glm.obj$form, data=dboot, family=binomial(link=logit))
				  .d1 <- .d0 <- dboot[dboot[,group] == 0,]
				  .d1[,group] <- 1; .d0[,group] <- 0
				  p1.s[i] <- mean(predict(glm.s, newdata=.d1, type="response", na.action=na.omit))
				  p0.s[i] <- mean(predict(glm.s, newdata=.d0, type="response", na.action=na.omit))
				  logOR.s[i] <- log((p1.s[i]/(1-p1.s[i]))/(p0.s[i]/(1-p0.s[i])))
				}
			}
		}
	}


	if(n.cluster>1){
	
		if(effect=="ATE"){
			bsimul <- function(glm.obj){
			dboot <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]
			glm.s <- glm(glm.obj$form, data=dboot, family=binomial(link=logit))
			.d1 <- .d0 <- dboot
			.d1[,group] <- 1; .d0[,group] <- 0
			p1.s <- mean(predict(glm.s, newdata=.d1, type="response", na.action=na.omit))
			p0.s <- mean(predict(glm.s, newdata=.d0, type="response", na.action=na.omit))
			return(c(logOR.s = log((p1.s/(1-p1.s))/(p0.s/(1-p0.s))),  p1.s = p1.s, p0.s = p0.s))
			}
		}else{
			if(effect=="ATT"){
				bsimul <- function(glm.obj){
				dboot <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]
				glm.s <- glm(glm.obj$form, data=dboot, family=binomial(link=logit))
				.d1 <- .d0 <- dboot[dboot[,group] == 1,]
				.d1[,group] <- 1; .d0[,group] <- 0
				p1.s <- mean(predict(glm.s, newdata=.d1, type="response", na.action=na.omit))
				p0.s <- mean(predict(glm.s, newdata=.d0, type="response", na.action=na.omit))
				return(c(logOR.s = log((p1.s/(1-p1.s))/(p0.s/(1-p0.s))),  p1.s = p1.s, p0.s = p0.s))
				}
			}else{
				bsimul <- function(glm.obj){
				dboot <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]
				glm.s <- glm(glm.obj$form, data=dboot, family=binomial(link=logit))
				.d1 <- .d0 <- dboot[dboot[,group] == 0,]
				.d1[,group] <- 1; .d0[,group] <- 0
				p1.s <- mean(predict(glm.s, newdata=.d1, type="response", na.action=na.omit))
				p0.s <- mean(predict(glm.s, newdata=.d0, type="response", na.action=na.omit))
				return(c(logOR.s = log((p1.s/(1-p1.s))/(p0.s/(1-p0.s))),  p1.s = p1.s, p0.s = p0.s))
				}
			}
		}
	  
	  cl <- makeCluster(n.cluster, type=cluster.type) 
		registerDoParallel(cl)
    clusterEvalQ(cl, {library(splines)})
		res <- NULL
		res <- foreach(i = 1:iterations, .combine=rbind, .inorder=TRUE) %dopar% {bsimul(glm.obj)}
		registerDoSEQ()
		stopCluster(cl) 

		logOR.s <- as.numeric(res[,"logOR.s"])
		p1.s <- as.numeric(res[,"p1.s"])
		p0.s <- as.numeric(res[,"p0.s"])
	}
}

se.logOR <- sd(logOR.s, na.rm=TRUE)

pv <- function(m, s){
  ztest <- m/s
  return(ifelse(ztest<0,2*pnorm(ztest),2*(1-pnorm(ztest))))
}

p.value.OR <- pv(m=logOR, s=se.logOR)

if(p.value.OR==0){ p.value.OR <- "<0.001" }

ci.low.logOR <- logOR - qnorm(0.975, 0, 1)*se.logOR
ci.upp.logOR <- logOR + qnorm(0.975, 0, 1)*se.logOR

ci.low.p0 <-  p0 -  qnorm(0.975, 0, 1)*sd(p0.s, na.rm=TRUE)
ci.upp.p0 <-  p0 +  qnorm(0.975, 0, 1)*sd(p0.s, na.rm=TRUE)

ci.low.p1 <- p1 -  qnorm(0.975, 0, 1)*sd(p1.s, na.rm=TRUE)
ci.upp.p1 <- p1 +  qnorm(0.975, 0, 1)*sd(p1.s, na.rm=TRUE)

delta.s <- p1.s - p0.s
se.delta <- sd(delta.s, na.rm=TRUE)

ci.low.delta <- (p1-p0) - qnorm(0.975, 0, 1)*se.delta
ci.upp.delta <- (p1-p0) + qnorm(0.975, 0, 1)*se.delta

return(list(
effect=effect,
p0=data.frame(estimate=p0, ci.lower=ci.low.p0, ci.upper=ci.upp.p0),
p1=data.frame(estimate=p1, ci.lower=ci.low.p1, ci.upper=ci.upp.p1),
delta=data.frame(estimate=p1-p0, std.error = se.delta, ci.lower=ci.low.delta,
 ci.upper=ci.upp.delta),
logOR=data.frame(estimate=logOR, std.error = se.logOR, ci.lower=ci.low.logOR,
 ci.upper=ci.upp.logOR),
p.value=p.value.OR) )

}

