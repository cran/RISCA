
gc.logistic<-function(glm.obj, data, group, effect="ATE", var.method="simulations", iterations=1000)
{

if(glm.obj$method != "glm.fit"){  
	stop("The argument \"glm.obj\" need to be a glm object") }
	
if(!is.data.frame(data) & !is.matrix(data)){
	stop("The argument \"data\" need to be a data.frame or a matrix") }
	
if(!is.character(group) & !is.numeric(group)){
	stop("The argument \"group\" need to be scalar or a character string") }
	
mod <- unique(data[,group]) 
if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
	stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the argument \"group\" ") }

if(is.na(match(effect,c("ATT","ATU","ATE"))) | length(effect)!=1){
	stop("Incorrect modality specified in the argument \"effect\": only one option among \"ATE\", \"ATT\" or \"ATU\" is possible") }
	
if(is.na(match(var.method,c("simulations","bootstrap"))) | length(var.method)!=1){
	stop("Incorrect modality specified in the argument \"var.method\": only one option among \"simulations\" and \"bootstrap\" is possible") }

if(length(grep("$", names(glm.obj$model), fixed = TRUE)) > 0 | length(grep("[", names(glm.obj$model), fixed = TRUE)) > 0){
	stop("Incorrect formula specified in the argument \"glm.obj\": don't use the syntax data$var or data[,var]") } #a reformuler
	
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

for(i in 1:iterations)
{
  glm.obj$coefficients <- simul[i,]
  p1.s[i] <- mean(predict(glm.obj, newdata=d1, type="response", na.action=na.omit))
  p0.s[i] <- mean(predict(glm.obj, newdata=d0, type="response", na.action=na.omit))
  logOR.s[i] <- log((p1.s[i]/(1-p1.s[i]))/(p0.s[i]/(1-p0.s[i])))
}
}
  
if(var.method=="bootstrap")
{

for(i in 1:iterations)
{
  j <- sample(1:nrow(data), size=nrow(data), replace = TRUE)
  .d1 <- data[j,]
  glm.s <- glm(glm.obj$form, data=.d1, family=binomial(link=logit))
  if(effect=="ATE"){ .ttt <- which(data[,group] == 0 | data[,group]==1) }
  if(effect=="ATT"){ .ttt <- which(data[,group] == 1) }
  if(effect=="ATU"){ .ttt <- which(data[,group] == 0) }
  .d1 <- .d0 <- .d1[.ttt,]
  .d1[,group] <- 1; .d0[,group] <- 0
  p1.s[i] <- mean(predict(glm.s, newdata=.d1, type="response", na.action=na.omit))
  p0.s[i] <- mean(predict(glm.s, newdata=.d0, type="response", na.action=na.omit))
  logOR.s[i] <- log((p1.s[i]/(1-p1.s[i]))/(p0.s[i]/(1-p0.s[i])))
}
}

se.logOR <- sd(logOR.s, na.rm=TRUE)


pv <- function(x, iterations){
	if(mean(x)<0){ 2*(sum(x>0)/iterations) } else {2*(sum(x<0)/iterations) } }
	
p.value.OR <- pv(logOR.s, iterations=iterations)

ci.low.logOR <- quantile(logOR.s, probs=c(0.025), na.rm=T)
ci.upp.logOR <- quantile(logOR.s, probs=c(0.975), na.rm=T)

ci.low.p0 <- quantile(p0.s, probs=c(0.025), na.rm=T)
ci.upp.p0 <- quantile(p0.s, probs=c(0.975), na.rm=T)

ci.low.p1 <- quantile(p1.s, probs=c(0.025), na.rm=T)
ci.upp.p1 <- quantile(p1.s, probs=c(0.975), na.rm=T)

delta.s <- p1.s - p0.s
se.delta <- sd(delta.s, na.rm=TRUE)

ci.low.delta <- quantile(delta.s, probs=c(0.025), na.rm=T)
ci.upp.delta <- quantile(delta.s, probs=c(0.975), na.rm=T)


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

