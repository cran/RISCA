
gc.survival<-function(coxph.obj, data, group, times, failures, max.time, effect="ATE", iterations=1000){

### Controls - Start

if(coxph.obj$call[1] != "coxph()"){  
	stop("The argument \"coxph.obj\" needs to be a coxph object") }

if(length(unlist(strsplit(as.character(coxph.obj$formula[2]), ","))) != 2){
	stop("The argument \"coxph.obj\" needs to be a proportional hazard model")}	

form <- coxph.obj$formula
	
if(!is.data.frame(data) & !is.matrix(data)){
	stop("The argument \"data\" needs to be a data.frame or a matrix") }
	
if(!is.character(group) & !is.numeric(group)){
	stop("The argument \"group\" needs to be scalar or a character string") }
	
if(length(grep("$", form, fixed = TRUE)) > 0 | length(grep("[", form, fixed = TRUE)) > 0){
	stop("Incorrect formula specified in the argument \"coxph.obj\": don't use the syntax data$var or data[,var]") } #a reformuler
	
if(!is.character(times)){
	if(is.numeric(times)){ 
	times <- colnames(data)[times]
	}else{ stop("The argument \"times\" needs to be scalar or a character string")  } }
	
mod <- unique(data[,group]) 
if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
	stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the argument \"group\" ") }
	
mod <- unique(data[,failures]) 
if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
	stop("Two modalities encoded 0 (for censored patients) and 1 (for uncensored patients) are required in the argument \"failures\" ") }

if(is.na(match(effect,c("ATT","ATU","ATE"))) | length(effect)!=1){
	stop("Incorrect modality specified in the argument \"effect\": only one option among \"ATE\", \"ATT\" or \"ATU\" are possible") }
	
if(length(names(coxph.obj$coef)) < 2){
	stop("Incorrect formula specified in the argument \"coxph.obj\": at least 2 covariables, including exposure, are requiered") }

if(is.na(match(group,names(coxph.obj$coef))) & is.na(match(names(data)[group],names(coxph.obj$coef)))){
	stop("Incorrect formula specified in the argument \"coxph.obj\": exposure is requiered") } 
	
if(max.time > max(data[,times], na.rm=T) | max.time < min(data[,times], na.rm=T)){
	stop("The argument \"max.time\" needs to be a number includes between the min and the max of the vector \"times\" ") }	

	
### Controls - End

finite.differences <- function(x, y) {
  n <- length(x)
  fdx0 <- vector(length = n)
  for (i in 2:n) {
    fdx0[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i])   }
    fdx0[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  return(fdx0)  }

  
rmst <- function(times, surv, max.time) {
        .t <- c(0, times[times <= max.time], min(max.time, max(times)))
        .s <- c(1, surv[times <= max.time], surv[length(surv)])
        return(sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)]))
    }


 
if(effect=="ATE"){ .ttt <- which(data[,group] == 0 | data[,group]==1)
}else if(effect=="ATT"){ .ttt <- which(data[,group] == 1)
}else .ttt <- which(data[,group] == 0 ) 

d1 <- d0 <- data[.ttt,]
d1[,group] <- 1 
d0[,group] <- 0

D1 <- data.frame(
	S = colMeans(riskRegression::predictCox(object=coxph.obj, newdata=d1, times=sort(unique(data[data[,failures]==1,times])), se = FALSE, iid = FALSE, keep.times = FALSE, type = "survival")$survival),
	T = sort(unique(data[data[,failures]==1,times]))	)
	
D0 <- data.frame(
	S = colMeans(riskRegression::predictCox(object=coxph.obj, newdata=d0, times=sort(unique(data[data[,failures]==1,times])), se = FALSE, iid = FALSE, keep.times = FALSE, type = "survival")$survival),
	T = sort(unique(data[data[,failures]==1,times]))	)

D1$h <- finite.differences(D1$T, -1*log(D1$S))
D0$h <- finite.differences(D0$T, -1*log(D0$S))

logHR <- log(mean(D1$h/D0$h))

RMST1 <- rmst(D1$T, D1$S, max.time)
RMST0 <- rmst(D0$T, D0$S, max.time)


logHR.s <- RMST1.s <- RMST0.s <- rep(-99, iterations)
for(i in 1:iterations){
	j <- sample(1:nrow(data), size=nrow(data), replace = TRUE)
	.d1 <- data[j,]
	cox.s <- coxph(form, data=.d1, x=TRUE)
	if(effect=="ATE"){ .ttt <- which(data[,group] == 0 | data[,group]==1)
	}else if(effect=="ATT"){ .ttt <- which(data[,group] == 1)
	}else .ttt <- which(data[,group] == 0 )
	.d1 <- .d0 <- .d1[.ttt,]
    .d1[,group] <- 1; .d0[,group] <- 0
	
	.D1 <- data.frame(
		S = colMeans(riskRegression::predictCox(object=cox.s, newdata=.d1, times=sort(unique(.d1[.d1[,failures]==1,times])), se = FALSE, iid = FALSE, keep.times = FALSE, type = "survival")$survival),
		T = sort(unique(.d1[.d1[,failures]==1,times]))	)
	
	.D0 <- data.frame(
		S = colMeans(riskRegression::predictCox(object=cox.s, newdata=.d0, times=sort(unique(.d0[.d0[,failures]==1,times])), se = FALSE, iid = FALSE, keep.times = FALSE, type = "survival")$survival),
		T = sort(unique(.d0[.d0[,failures]==1,times]))	)

	.D1$h <- finite.differences(.D1$T, -1*log(.D1$S))
	.D0$h <- finite.differences(.D0$T, -1*log(.D0$S))

	logHR.s[i] <- log(mean(.D1$h/.D0$h))

	RMST1.s[i] <- rmst(.D1$T, .D1$S, max.time)
	RMST0.s[i] <- rmst(.D0$T, .D0$S, max.time)
}

se.logHR <- sd(logHR.s, na.rm=TRUE)

pv <- function(x, iterations){
	if(mean(x)<0){ 2*(sum(x>0)/iterations) } else {2*(sum(x<0)/iterations) } }
	
p.value.HR <- pv(logHR.s, iterations=iterations)

ci.low.logHR <- quantile(logHR.s, probs=c(0.025), na.rm=T)
ci.upp.logHR <- quantile(logHR.s, probs=c(0.975), na.rm=T)

ci.low.RMST1 <- quantile(RMST1.s, probs=c(0.025), na.rm=T)
ci.upp.RMST1 <- quantile(RMST1.s, probs=c(0.975), na.rm=T)

ci.low.RMST0 <- quantile(RMST0.s, probs=c(0.025), na.rm=T)
ci.upp.RMST0 <- quantile(RMST0.s, probs=c(0.975), na.rm=T)

delta.s <- RMST1.s - RMST0.s
se.delta <- sd(delta.s, na.rm=TRUE)

p.value.delta <- pv(delta.s, iterations=iterations)

ci.low.delta <- quantile(delta.s, probs=c(0.025), na.rm=T)
ci.upp.delta <- quantile(delta.s, probs=c(0.975), na.rm=T)

return(list(
effect=effect,
max.time=max.time,
RMST0=data.frame(estimate=RMST0, ci.lower=ci.low.RMST0, ci.upper=ci.upp.RMST0),
RMST1=data.frame(estimate=RMST1, ci.lower=ci.low.RMST1, ci.upper=ci.upp.RMST1),
delta=data.frame(estimate=RMST1-RMST0, std.error = se.delta, ci.lower=ci.low.delta,
 ci.upper=ci.upp.delta, p.value=p.value.delta),
logHR=data.frame(estimate=logHR, std.error = se.logHR, ci.lower=ci.low.logHR,
 ci.upper=ci.upp.logHR, p.value=p.value.HR) ) )
}



