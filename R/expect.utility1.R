
expect.utility1 <- function(times, failures, variable, pro.time, u.A0, u.A1, u.B0, u.B1, n.boot=NULL, rmst.change)
{
	
    if ( pro.time> max(times, na.rm=TRUE) ) {
        stop("The prognostic time is higher than the maximum observed survival time.")
    }
    if ( is.null(u.A0)|is.null(u.A1)|is.null(u.B0)|is.null(u.B1) ) {
        stop("At least one utility is missing.")
    }
    if ( ((length(times)+length(failures)+length(variable))/3 ) != min(length(times),length(failures),length(variable))) {
        stop("The lengths of the arguments times, failures and variable have to be equalled.")
    }
	if  ( length(u.A0)!=1 | length(u.A1)!=1 | length(u.B0)!=1 | length(u.B1)!=1)  {
        stop("The lengths of utilities have to be 1.")
    }
	if (is.na(rmst.change))  {
        stop("The Restricted Mean Survival Time (RMST) change is missing.")
    }
	
cut.off <- sort(unique(variable))

.temp.data <- data.frame(times, failures, variable)
.temp.data$failures[.temp.data$times > pro.time] <- 0
.temp.data$times[.temp.data$times > pro.time] <- pro.time + 0.01
.na <- is.na(.temp.data$times + .temp.data$failures + .temp.data$variable)
.n.na <- sum(.na)
.temp.data <- .temp.data[.na==FALSE, ]

meanA <- function(x) {
    .data <- .temp.data[.temp.data$variable > x , ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)]==0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- ( sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)]) ) * (1 + rmst.change)
         return(.res)  }  }  }

meanB <- function(x) {
    .data <- .temp.data[.temp.data$variable <= x , ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)] == 0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])
         return(.res)  }  } }

EA <- pmin(sapply(cut.off, FUN = "meanA"), pro.time)
EA[cut.off >= min(c(cut.off[is.na(EA)], (max(cut.off) + 1)))] <- NA

EB <- pmin(sapply(cut.off, FUN = "meanB"), pro.time)
EB[cut.off < max( c( min(cut.off)-1 , cut.off[is.na(EB)] ))] <- NA

qA <- u.A0 * EA + u.A1 * (pro.time - EA)
qB <- u.B0 * EB + u.B1 * (pro.time - EB)

p <- function(x) {mean(.temp.data$variable > x)}
pA <- sapply(cut.off, FUN = "p")
pB <- 1-pA

esp.res <- data.frame(cut.off = cut.off, eA=EA, eB=EB)
tab.res <- data.frame(cut.off = cut.off, utility = pA * qA + pB * qB)
prob.res <- data.frame(cut.off = cut.off, pA = pA, pB = pB)
qaly.res <- data.frame(cut.off = cut.off, qA = qA, qB = qB)

esp.res <- esp.res[!is.na(tab.res$utility), ]
prob.res <- prob.res[!is.na(tab.res$utility), ]
qaly.res <- qaly.res[!is.na(tab.res$utility), ]
tab.res <- tab.res[!is.na(tab.res$utility), ]

est <- tab.res$cut.off[max(tab.res$utility, na.rm = TRUE)==tab.res$utility & !is.na(tab.res$utility)][1]
m <- round(max(tab.res$utility, na.rm = TRUE), 10)

meanALL <- function(times, failures, pro.time) {
.km <- summary(survfit(Surv(times, failures) ~ 1))
.t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(times)))
.s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
return(sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])) }

meanALL.A <-pmin( meanALL(.temp.data$times, .temp.data$failures, pro.time) * (1+rmst.change), pro.time)
eu.allA <- u.A0 * meanALL.A + u.A1 * (pro.time - meanALL.A)

meanALL.B <- meanALL(.temp.data$times, .temp.data$failures, pro.time)
eu.allB <- u.B0 * meanALL.B + u.B1 * (pro.time - meanALL.B)


m.correct <- max(m, eu.allA, eu.allB)
if (m.correct!=m){
    if (eu.allB < eu.allA){ est <- min(.temp.data$variable)  }
      else{est <- max(.temp.data$variable)}
} else {
    if (m==eu.allA){est <- min(.temp.data$variable)  }
      else{if (m==eu.allB){est <- max(.temp.data$variable)  }}}

esp.res <- rbind(c(cut.off=-Inf, eA=meanALL.A, eB=NA), esp.res, c(cut.off=+Inf, eA=NA, eB=meanALL.B))
prob.res <- rbind(c(cut.off=-Inf, pA=1, pB=0), prob.res, c(cut.off=+Inf, pA=0, pB=1))
qaly.res <- rbind(c(cut.off=-Inf, qA=eu.allA , qB=NA), qaly.res, c(cut.off=+Inf, qA=NA, qB=eu.allB))
tab.res <- rbind(c(cut.off=-Inf, utility=eu.allA), tab.res, c(cut.off=+Inf, utility=eu.allB))

meanAp <- function(x) {
    .data <- .temp.data[.temp.data$variable > x , ]
    if (dim(.data)[1] == 0) {return(NA)}
    else {
      if ((max(.data$times) < pro.time) && (.data$failures[.data$times == max(.data$times)]==0)) {return(NA)}
       else {
         .km <- summary(survfit(Surv(times, failures) ~ 1, data = .data))
         .t <- c(0, .km$time[.km$time <= pro.time], min(pro.time, max(.data$times)))
         .s <- c(1, .km$surv[.km$time <= pro.time], .km$surv[length(.km$surv)])
         .res <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])
         return(.res)  }  }  }

EAp <- meanAp(est)

qAp <- u.B0 * EAp + u.B1 * (pro.time - EAp)

if(est==min(.temp.data$variable)){
    esp<-esp.res[1,]$eA
    qaly<-qaly.res[1,]$qA
}else{
if(est==max(.temp.data$variable)){
    esp<-esp.res[length(esp.res$eA),]$eA
    qaly<-qaly.res[length(qaly.res$qA),]$qA

}else{
    esp<-esp.res$eA[qaly.res$cut.off==est]
    qaly<-qaly.res$qA[qaly.res$cut.off==est]
}}


est.res<-est

if (!is.null(n.boot)){

.temp.data$ident<-1:(dim(.temp.data)[1])
sgde<-.temp.data
res.boot<-rep(-99, n.boot)

for (b in 1:n.boot)
{

.temp.data<-sgde[sample(sgde$ident, size=dim(sgde)[1], replace = TRUE),]

.temp.x.boot<-sort(unique(.temp.data$variable))

EA.b <- pmin(sapply(.temp.x.boot, FUN = "meanA"), pro.time)
EA.b[.temp.x.boot >= min(c(.temp.x.boot[is.na(EA.b)], (max(.temp.x.boot) + 1)))] <- NA

EB.b <- pmin(sapply(.temp.x.boot, FUN = "meanB"), pro.time)
EB.b[.temp.x.boot < max( c( min(.temp.x.boot)-1 ,.temp.x.boot[is.na(EB.b)] ))] <- NA

qA.b <- u.A0 * EA.b + u.A1 * (pro.time - EA.b)
qB.b <- u.B0 * EB.b + u.B1 * (pro.time - EB.b)

pA.b <- sapply(.temp.x.boot, FUN = "p")
pB.b <- 1-pA.b

tab.res.b <- data.frame(cut.off = .temp.x.boot, utility = pA.b * qA.b + pB.b * qB.b)
tab.res.b <- tab.res.b[!is.na(tab.res.b$utility), ]

est.b <- tab.res.b$cut.off[max(tab.res.b$utility, na.rm = TRUE)==tab.res.b$utility & !is.na(tab.res.b$utility)][1]
m.b <- round(max(tab.res.b$utility, na.rm = TRUE), 10)

meanALL.A.b <-pmin( meanALL(.temp.data$times, .temp.data$failures, pro.time) * (1+rmst.change),pro.time)
eu.allA.b <- u.A0 * meanALL.A.b + u.A1 * (pro.time - meanALL.A.b)

meanALL.B.b <- meanALL(.temp.data$times, .temp.data$failures, pro.time)
eu.allB.b <- u.B0 * meanALL.B.b + u.B1 * (pro.time - meanALL.B.b)

m.correct.b <- max(m.b, eu.allA.b, eu.allB.b)
if (m.correct.b!=m.b){
    if (eu.allB.b < eu.allA.b){ est.b <- min(.temp.data$variable)  }
      else{est.b <- max(.temp.data$variable)} 
}else{
	if (m.b==eu.allA.b){est.b <- min(.temp.data$variable)  }
      else{if (m.b==eu.allB.b){est.b <- max(.temp.data$variable)  }}}
	  
res.boot[b] <- est.b

}

est.res<-c(estimation = est, binf = quantile(res.boot, probs=0.025), bsup = quantile(res.boot, probs=0.975))

}

return(list(
 estimation = est.res,
 max.eu = m.correct,
 table = cbind(tab.res, prob.res[,c("pA", "pB")], qaly.res[,c("qA", "qB")], esp.res[,c("eA", "eB")]),
 delta.rmst = esp - EAp,
 delta.qaly = qaly - qAp,
 missing = .n.na ) )
}