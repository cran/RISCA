gc.survival<-function(coxph.obj, data, group, times, failures, max.time,
effect="ATE", iterations=1000, n.cluster=1){

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

if(length(grep("$", form, fixed = TRUE)) > 0 | length(grep("[", form,
fixed = TRUE)) > 0){
	stop("Incorrect formula specified in the argument \"coxph.obj\": don't
use the syntax data$var or data[,var]") } #a reformuler

if(!is.character(times)){
	if(is.numeric(times)){
	times <- colnames(data)[times]
	}else{ stop("The argument \"times\" needs to be scalar or a character
string")  } }

mod <- unique(data[,group])
if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2]
!= 0))){
	stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and
1 (for treated/exposed patients) are required in the argument \"group\"
") }

mod <- unique(data[,failures])
if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2]
!= 0))){
	stop("Two modalities encoded 0 (for censored patients) and 1 (for
uncensored patients) are required in the argument \"failures\" ") }

if(is.na(match(effect,c("ATT","ATU","ATE"))) | length(effect)!=1){
	stop("Incorrect modality specified in the argument \"effect\": only one
option among \"ATE\", \"ATT\" or \"ATU\" are possible") }

if(length(names(coxph.obj$coef)) < 2){
	stop("Incorrect formula specified in the argument \"coxph.obj\": at least
2 covariables, including exposure, are requiered") }

if(is.na(match(group,names(coxph.obj$coef))) &
is.na(match(names(data)[group],names(coxph.obj$coef)))){
	stop("Incorrect formula specified in the argument \"coxph.obj\": exposure
is requiered") }

if(max.time > max(data[,times], na.rm=T) | max.time < min(data[,times],
na.rm=T)){
	stop("The argument \"max.time\" needs to be a number includes between the
min and the max of the vector \"times\" ") }

if(!is.null(attr(coxph.obj$terms,"specials")$strata) |
!is.null(attr(coxph.obj$terms,"specials")$tt) |
!is.null(attr(coxph.obj$terms,"specials")$cluster)){
	stop("Incorrect argument in the argument \"coxph.obj\": time-transform
functions, stratification and clustering are not implemented") }

ties <- coxph.obj$method

### Controls - End

explp <-function(coxph.obj, newdata){
	return( exp(predict(coxph.obj, newdata, type="lp") +
sum(setNames(coxph.obj$means, names(coef(coxph.obj)))*coef(coxph.obj))) )
}

if(effect=="ATE"){ ttt <- which(data[,group] %in% c(0,1))
}else if(effect=="ATT"){ ttt <- which(data[,group] == 1)
}else ttt <- which(data[,group] == 0 )

d1 <- d0 <- d <- data[ttt,]
d1[,group] <- 1
d0[,group] <- 0

b <- basehaz(coxph.obj, centered=FALSE)
H0.multi <- b$hazard[b$time %in% sort(unique(d[d[,failures]==1,times]))]
T.multi  <- b$time[b$time %in% sort(unique(d[d[,failures]==1,times]))]

D1 <- data.frame(
	S = colMeans(exp(explp(coxph.obj,d1) %o% -H0.multi )),
	T = T.multi)

D0 <- data.frame(
	S = colMeans(exp(explp(coxph.obj,d0) %o% -H0.multi )),
	T = T.multi)

D1$h <- differentiation(x=D1$T, fx=-1*log(D1$S))
D0$h <- differentiation(x=D0$T, fx=-1*log(D0$S))

logHR <- log(mean(D1$h/D0$h))

table.surv <- data.frame(
	times=D1$T,
	survival=c(D1$S,D0$S),
	n.risk=c(
		sapply(D1$T, function(temps) {
			sum(d[d[,group]==1,times] >= temps)} ),
		sapply(D0$T, function(temps) {
			sum(d[d[,group]==0,times] >= temps)})),
	variable=c(
		rep(1,dim(D1)[1]),
		rep(0,dim(D0)[1])) )

RMST1 <- rmst(D1$T, D1$S, max.time)
RMST0 <- rmst(D0$T, D0$S, max.time)

logHR.s <- RMST1.s <- RMST0.s <- rep(-99, iterations)

if(n.cluster==1){

	if(effect=="ATE"){
		for(i in 1:iterations){
			data.b <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]

			cox.s <- coxph(form, data=data.b, ties=ties, x=TRUE, y=TRUE)
			sfit <- survfit(cox.s, data=data.b, se.fit=FALSE)

			.d1 <- .d0 <- data.b
			.d1[,group] <- 1
			.d0[,group] <- 0

			.H0.multi <- sfit$cumhaz[sfit$n.event>0] * exp(-sum(cox.s$means *
				coef(cox.s)))
			.T.multi  <- sfit$time[sfit$n.event>0]

			.D1 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d1) %o% -.H0.multi )),	T = .T.multi)
			.D0 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d0) %o% -.H0.multi )),	T = .T.multi)

			logHR.s[i] <- log( mean(
				differentiation(x=.D1$T, fx=-1*log(.D1$S)) /
				differentiation(x=.D0$T, fx=-1*log(.D0$S))) )
			RMST1.s[i] <- rmst(.D1$T, .D1$S, max.time)
			RMST0.s[i] <- rmst(.D0$T, .D0$S, max.time)
		}
	}else{
		if(effect=="ATT"){
			for(i in 1:iterations){
			data.b <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]

			cox.s <- coxph(form, data=data.b, ties=ties, x=TRUE, y=TRUE)
			sfit <- survfit(cox.s, data=data.b, se.fit=FALSE)

			.d1 <- .d0 <- data.b[data.b[,group] == 1,]
			.d1[,group] <- 1
			.d0[,group] <- 0

			.H0.multi <- sfit$cumhaz[sfit$n.event>0] * exp(-sum(cox.s$means *
				coef(cox.s)))
			.T.multi  <- sfit$time[sfit$n.event>0]

			.D1 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d1) %o% -.H0.multi )),	T = .T.multi)
			.D0 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d0) %o% -.H0.multi )),	T = .T.multi)

			logHR.s[i] <- log( mean(
				differentiation(x=.D1$T, fx=-1*log(.D1$S)) /
				differentiation(x=.D0$T, fx=-1*log(.D0$S))) )
			RMST1.s[i] <- rmst(.D1$T, .D1$S, max.time)
			RMST0.s[i] <- rmst(.D0$T, .D0$S, max.time)
			}
		}else{
			for(i in 1:iterations){
			data.b <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]

			cox.s <- coxph(form, data=data.b, ties=ties, x=TRUE, y=TRUE)
			sfit <- survfit(cox.s, data=data.b, se.fit=FALSE)

			.d1 <- .d0 <- data.b[data.b[,group] == 0,]
			.d1[,group] <- 1
			.d0[,group] <- 0

			.H0.multi <- sfit$cumhaz[sfit$n.event>0] * exp(-sum(cox.s$means *
				coef(cox.s)))
			.T.multi  <- sfit$time[sfit$n.event>0]

			.D1 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d1) %o% -.H0.multi )),	T = .T.multi)
			.D0 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d0) %o% -.H0.multi )),	T = .T.multi)

			logHR.s[i] <- log( mean(
				differentiation(x=.D1$T, fx=-1*log(.D1$S)) /
				differentiation(x=.D0$T, fx=-1*log(.D0$S))) )
			RMST1.s[i] <- rmst(.D1$T, .D1$S, max.time)
			RMST0.s[i] <- rmst(.D0$T, .D0$S, max.time)
			}
		}
	}
}

if(n.cluster>1){

	if(effect=="ATE"){
		simul <- function(){
			data.b <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]

			cox.s <- coxph(form, data=data.b, ties=ties, x=TRUE, y=TRUE)
			sfit <- survfit(cox.s, data=data.b, se.fit=FALSE)

			.d1 <- .d0 <- data.b
			.d1[,group] <- 1
			.d0[,group] <- 0

			.H0.multi <- sfit$cumhaz[sfit$n.event>0] * exp(-sum(cox.s$means *
coef(cox.s)))
			.T.multi  <- sfit$time[sfit$n.event>0]

			.D1 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d1) %o% -.H0.multi )),	T = .T.multi)
			.D0 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d0) %o% -.H0.multi )),	T = .T.multi)

			return(c(logHR.s = log(mean(differentiation(x=.D1$T,
fx=-1*log(.D1$S))/differentiation(x=.D0$T, fx=-1*log(.D0$S)))), RMST1.s
= rmst(.D1$T, .D1$S, max.time), RMST0.s = rmst(.D0$T, .D0$S,
max.time)))
		}
	}else{
		if(effect=="ATT"){
			simul <- function(){
			data.b <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]

			cox.s <- coxph(form, data=data.b, ties=ties, x=TRUE, y=TRUE)
			sfit <- survfit(cox.s, data=data.b, se.fit=FALSE)

			.d1 <- .d0 <- data.b[data.b[,group] == 1,]
			.d1[,group] <- 1
			.d0[,group] <- 0

			.H0.multi <- sfit$cumhaz[sfit$n.event>0] * exp(-sum(cox.s$means *
coef(cox.s)))
			.T.multi  <- sfit$time[sfit$n.event>0]

			.D1 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d1) %o% -.H0.multi )),	T = .T.multi)
			.D0 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d0) %o% -.H0.multi )),	T = .T.multi)

			return(c(logHR.s = log(mean(differentiation(x=.D1$T,
fx=-1*log(.D1$S))/differentiation(x=.D0$T, fx=-1*log(.D0$S)))), RMST1.s
= rmst(.D1$T, .D1$S, max.time), RMST0.s = rmst(.D0$T, .D0$S,
max.time)))
			}
		}else{
			simul <- function(){
			data.b <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]

			cox.s <- coxph(form, data=data.b, ties=ties, x=TRUE, y=TRUE)
			sfit <- survfit(cox.s, data=data.b, se.fit=FALSE)

			.d1 <- .d0 <- data.b[data.b[,group] == 0,]
			.d1[,group] <- 1
			.d0[,group] <- 0

			.H0.multi <- sfit$cumhaz[sfit$n.event>0] * exp(-sum(cox.s$means *
coef(cox.s)))
			.T.multi  <- sfit$time[sfit$n.event>0]

			.D1 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d1) %o% -.H0.multi )),	T = .T.multi)
			.D0 <- data.frame(
				S = colMeans(exp(explp(cox.s,.d0) %o% -.H0.multi )),	T = .T.multi)

			return(c(logHR.s = log(mean(differentiation(x=.D1$T,
fx=-1*log(.D1$S))/differentiation(x=.D0$T, fx=-1*log(.D0$S)))), RMST1.s
= rmst(.D1$T, .D1$S, max.time), RMST0.s = rmst(.D0$T, .D0$S,
max.time)))
			}
		}
	}
# ouverture des coeurs
if(.Platform[[1]] == "windows"){
	cl <- makeCluster(n.cluster, type="PSOCK")
}else{
	cl <- makeCluster(n.cluster, type="FORK")
}
registerDoParallel(cl)
clusterEvalQ(cl, {library(splines); library(survival)})
res <- NULL
res <- foreach(i = 1:iterations, .combine=rbind, .inorder=TRUE) %dopar%
{simul()}

# fermeture des coeurs
registerDoSEQ()
stopCluster(cl)

logHR.s <- as.numeric(res[,"logHR.s"])
RMST1.s <- as.numeric(res[,"RMST1.s"])
RMST0.s <- as.numeric(res[,"RMST0.s"])

}

se.logHR <- sd(logHR.s, na.rm=TRUE)

pv <- function(x, iterations){
	if(mean(x)<0){ 2*(sum(x>0)/iterations) } else {2*(sum(x<0)/iterations) } }

p.value.HR <- pv(logHR.s, iterations=iterations)
if(p.value.HR==0){ p.value.HR <- "<0.001" }

ci.low.logHR <- quantile(logHR.s, probs=c(0.025), na.rm=T)
ci.upp.logHR <- quantile(logHR.s, probs=c(0.975), na.rm=T)

ci.low.RMST1 <- quantile(RMST1.s, probs=c(0.025), na.rm=T)
ci.upp.RMST1 <- quantile(RMST1.s, probs=c(0.975), na.rm=T)

ci.low.RMST0 <- quantile(RMST0.s, probs=c(0.025), na.rm=T)
ci.upp.RMST0 <- quantile(RMST0.s, probs=c(0.975), na.rm=T)

delta.s <- RMST1.s - RMST0.s
se.delta <- sd(delta.s, na.rm=TRUE)

p.value.delta <- pv(delta.s, iterations=iterations)
if(p.value.delta==0){ p.value.delta <- "<0.001" }

ci.low.delta <- quantile(delta.s, probs=c(0.025), na.rm=T)
ci.upp.delta <- quantile(delta.s, probs=c(0.975), na.rm=T)

.obj <- list(
	table.surv=table.surv,
	effect=effect,
	max.time=data.frame(max.time=max.time),
	RMST0=data.frame(estimate=RMST0, ci.lower=ci.low.RMST0,
		ci.upper=ci.upp.RMST0),
	RMST1=data.frame(estimate=RMST1, ci.lower=ci.low.RMST1,
		ci.upper=ci.upp.RMST1),
	delta=data.frame(estimate=RMST1-RMST0, std.error = se.delta,
		ci.lower=ci.low.delta, ci.upper=ci.upp.delta, p.value=p.value.delta),
	logHR=data.frame(estimate=logHR, std.error = se.logHR, ci.lower=ci.low.logHR,
		ci.upper=ci.upp.logHR, p.value=p.value.HR) )

class(.obj) <- "survival"
return(.obj)
}

