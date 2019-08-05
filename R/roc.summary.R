utils::globalVariables(c(".s1", ".s2", ".s3", ".s4", "d.int", "d.surv4", "d.surv5", "p.joint.model4", "p.joint.model5", "se.sp.marg4", "se.sp.marg5"))

roc.summary<-function(study.num, classe, n, year, surv, nrisk, proba, marker.min, marker.max, init.nlme1=c(0,1), precision=10^-6, pro.time, time.cutoff)
{
if (missing(study.num)) stop("Argument 'study.num' is missing with no default")
if (missing(classe)) stop("Argument 'classe' is missing with no default")
if (missing(n)) stop("Argument 'n' is missing with no default")
if (missing(year)) stop("Argument 'year' is missing with no default")
if (missing(surv)) stop("Argument 'surv' is missing with no default")
if (missing(nrisk)) stop("Argument 'nrisk' is missing with no default")
if (missing(proba)) stop("Argument 'proba' is missing with no default")
if (missing(marker.min)) stop("Argument 'marker.min' is missing with no default")
if (missing(marker.max)) stop("Argument 'marker.max' is missing with no default")
if (missing(pro.time)) stop("Argument 'pro.time' is missing with no default")

if (!is.vector(time.cutoff) | !is.numeric(time.cutoff))
        stop("Argument 'time.cutoff' must be a numeric vector")

if (length(unique(time.cutoff))<3 | length(unique(time.cutoff))>4) stop("Argument 'time.cutoff' must contain 3 or 4 distinct values")

if (length(unique(time.cutoff))!=length(time.cutoff)) stop("Argument 'time.cutoff' must contain distinct values")

if (!is.vector(study.num) | !is.numeric(study.num))
        stop("Argument 'study.num' must be a numeric vector")

if (!is.vector(classe) | !is.numeric(classe))
        stop("Argument 'classe' must be a numeric vector")

if (!is.vector(n) | !is.numeric(n))
        stop("Argument 'n' must be a numeric vector")

if (!is.vector(year) | !is.numeric(year))
        stop("Argument 'year' must be a numeric vector")

if (!is.vector(surv) | !is.numeric(surv))
        stop("Argument 'surv' must be a numeric vector")

if (!is.vector(nrisk) | !is.numeric(nrisk))
        stop("Argument 'nrisk' must be a numeric vector")

if (!is.vector(proba) | !is.numeric(proba))
        stop("Argument 'study.num' must be a numeric vector")

if (!is.numeric(marker.min)) stop("Argument 'marker.min' must be a numeric vector")

if (!is.numeric(marker.max)) stop("Argument 'marker.max' must be a numeric vector")

if (!is.numeric(precision)) stop("Argument 'precision' must be a numeric vector")
if (length(precision)!=1) stop("Argument 'precision' must be a single numeric value")

if (!is.numeric(pro.time)) stop("Argument 'pro.time' must be a numeric vector")
if (length(pro.time)!=1) stop("Argument 'pro.time' must be a single numeric value")

if(length(time.cutoff)==4)
{
assign(".s1", time.cutoff[1], envir = parent.frame()) 
assign(".s2", time.cutoff[2], envir = parent.frame()) 
assign(".s3", time.cutoff[3], envir = parent.frame()) 
assign(".s4", time.cutoff[4], envir = parent.frame()) 

assign("d.int", function(mu, sigma, min1, max1)
{
noeuds<-gauss.quad(30, kind="legendre", alpha=0, beta=0)$nodes
poids<-gauss.quad(30, kind="legendre", alpha=0, beta=0)$weights

.expr1<-(max1-min1)
.expr2<-(max1+min1)
return(
    0.5*(.expr1)*(
    poids[1]*dnorm(0.5*.expr1*noeuds[1]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[2]*dnorm(0.5*.expr1*noeuds[2]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[3]*dnorm(0.5*.expr1*noeuds[3]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[4]*dnorm(0.5*.expr1*noeuds[4]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[5]*dnorm(0.5*.expr1*noeuds[5]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[6]*dnorm(0.5*.expr1*noeuds[6]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[7]*dnorm(0.5*.expr1*noeuds[7]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[8]*dnorm(0.5*.expr1*noeuds[8]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[9]*dnorm(0.5*.expr1*noeuds[9]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[10]*dnorm(0.5*.expr1*noeuds[10]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[11]*dnorm(0.5*.expr1*noeuds[11]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[12]*dnorm(0.5*.expr1*noeuds[12]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[13]*dnorm(0.5*.expr1*noeuds[13]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[14]*dnorm(0.5*.expr1*noeuds[14]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[15]*dnorm(0.5*.expr1*noeuds[15]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[16]*dnorm(0.5*.expr1*noeuds[16]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[17]*dnorm(0.5*.expr1*noeuds[17]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[18]*dnorm(0.5*.expr1*noeuds[18]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[19]*dnorm(0.5*.expr1*noeuds[19]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[20]*dnorm(0.5*.expr1*noeuds[20]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[21]*dnorm(0.5*.expr1*noeuds[21]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[22]*dnorm(0.5*.expr1*noeuds[22]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[23]*dnorm(0.5*.expr1*noeuds[23]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[24]*dnorm(0.5*.expr1*noeuds[24]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[25]*dnorm(0.5*.expr1*noeuds[25]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[26]*dnorm(0.5*.expr1*noeuds[26]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[27]*dnorm(0.5*.expr1*noeuds[27]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[28]*dnorm(0.5*.expr1*noeuds[28]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[29]*dnorm(0.5*.expr1*noeuds[29]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
	poids[30]*dnorm(0.5*.expr1*noeuds[30]+0.5*.expr2, mean=mu, sd=exp(sigma)) ) )
}, envir = parent.frame())


assign("d.surv5", function(z, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)
{
.l1 <- .s1
.l2 <- .s2 - .s1
.l3 <- .s3 - .s2
.l4 <- .s4 - .s3

.expr0<-beta0.1+beta0.2*z
.expr1<-beta1.1+beta1.2*z
.expr2<-beta2.1+beta2.2*z
.expr3<-beta3.1+beta3.2*z
.expr4<-beta4.1+beta4.2*z

.expr5<-(-t*exp(.expr0))*(t<=.s1) +
(-.l1*exp(.expr0) -(t-.s1)*exp(.expr0+.expr1))*(t>.s1)*(t<=.s2) +
(-.l1*exp(.expr0) - .l2*exp(.expr0+.expr1) -(t-.s2)*exp(.expr0+.expr1+.expr2))*(t>.s2)*(t<=.s3) +
(-.l1*exp(.expr0) - .l2*exp(.expr0+.expr1) - .l3*exp(.expr0+.expr1+.expr2) - (t-.s3)*exp(.expr0+.expr1+.expr2+.expr3))*(t>.s3)*(t<=.s4)+
(-.l1*exp(.expr0) - .l2*exp(.expr0+.expr1) - .l3*exp(.expr0+.expr1+.expr2) - .l4*exp(.expr0+.expr1+.expr2+.expr3) -(t-.s4)*exp(.expr0+.expr1+.expr2+.expr3+.expr4) )*(t>.s4)
return(exp(.expr5) * dnorm(z, mean=mu, sd=exp(sigma)))
}, envir = parent.frame())


assign("p.joint.model5", function(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma, t, min.z, max.z)
{
noeuds<-gauss.quad(30, kind="legendre", alpha=0, beta=0)$nodes
poids<-gauss.quad(30, kind="legendre", alpha=0, beta=0)$weights

.exprA<-0.5*(max.z-min.z)
.exprB<-0.5*(max.z+min.z)
return( .exprA*(
 poids[1]*d.surv5(.exprA*noeuds[1]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[2]*d.surv5(.exprA*noeuds[2]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[3]*d.surv5(.exprA*noeuds[3]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[4]*d.surv5(.exprA*noeuds[4]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[5]*d.surv5(.exprA*noeuds[5]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[6]*d.surv5(.exprA*noeuds[6]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[7]*d.surv5(.exprA*noeuds[7]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[8]*d.surv5(.exprA*noeuds[8]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[9]*d.surv5(.exprA*noeuds[9]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[10]*d.surv5(.exprA*noeuds[10]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[11]*d.surv5(.exprA*noeuds[11]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[12]*d.surv5(.exprA*noeuds[12]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[13]*d.surv5(.exprA*noeuds[13]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[14]*d.surv5(.exprA*noeuds[14]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[15]*d.surv5(.exprA*noeuds[15]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[16]*d.surv5(.exprA*noeuds[16]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[17]*d.surv5(.exprA*noeuds[17]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[18]*d.surv5(.exprA*noeuds[18]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[19]*d.surv5(.exprA*noeuds[19]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[20]*d.surv5(.exprA*noeuds[20]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[21]*d.surv5(.exprA*noeuds[21]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[22]*d.surv5(.exprA*noeuds[22]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[23]*d.surv5(.exprA*noeuds[23]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[24]*d.surv5(.exprA*noeuds[24]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[25]*d.surv5(.exprA*noeuds[25]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[26]*d.surv5(.exprA*noeuds[26]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[27]*d.surv5(.exprA*noeuds[27]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[28]*d.surv5(.exprA*noeuds[28]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[29]*d.surv5(.exprA*noeuds[29]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma)+ 
 poids[30]*d.surv5(.exprA*noeuds[30]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma) ) )
}, envir = parent.frame()) 

assign("se.sp.marg5", function(cut, tps)
{
.marker.min<-qnorm(0.000001, mean = nlme1$coefficients$fixed[1], sd = exp(nlme1$coefficients$fixed[2]), lower.tail = TRUE, log.p = FALSE)

.marker.max<-qnorm(0.999999, mean = nlme1$coefficients$fixed[1], sd = exp(nlme1$coefficients$fixed[2]), lower.tail = TRUE, log.p = FALSE)

est.mu    <- nlme1$coefficients$fixed[1]
est.sigma <- nlme1$coefficients$fixed[2]

est.beta0.1 <- nlme2$coefficients$fixed[1]
est.beta0.2 <- nlme2$coefficients$fixed[2]
est.beta1.1 <- nlme2$coefficients$fixed[3]
est.beta1.2 <- nlme2$coefficients$fixed[4]
est.beta2.1 <- nlme2$coefficients$fixed[5]
est.beta2.2 <- nlme2$coefficients$fixed[6]
est.beta3.1 <- nlme2$coefficients$fixed[7]
est.beta3.2 <- nlme2$coefficients$fixed[8]
est.beta4.1 <- nlme2$coefficients$fixed[9]
est.beta4.2 <- nlme2$coefficients$fixed[10]

exprA.<-d.int(mu=est.mu, sigma=est.sigma, cut, .marker.max)

exprB.<-p.joint.model5(est.beta0.1, est.beta0.2, est.beta1.1, est.beta1.2, est.beta2.1, est.beta2.2, est.beta3.1, est.beta3.2, est.beta4.1, est.beta4.2, est.mu, est.sigma, tps, cut, .marker.max)

exprC.<-p.joint.model5(est.beta0.1, est.beta0.2, est.beta1.1, est.beta1.2, est.beta2.1, est.beta2.2, est.beta3.1, est.beta3.2, est.beta4.1, est.beta4.2, est.mu, est.sigma, tps, .marker.min, .marker.max)

return( c(( exprA. - exprB.  ) / ( 1 - exprC. ) , 1- exprB. / exprC.) )

}, envir = parent.frame()) 

.data <- data.frame(classe, n, year, surv, nrisk,  proba, marker.min, marker.max, study.num)

.data$p.joint <- .data$surv * .data$proba

.eff <- .data[(.data$year==1), c("classe", "n", "marker.min", "marker.max", "study.num")]

.temp.func <- function(x) { sum(.eff$n[.eff$study.num==x]) }

.centre <- data.frame(study.num = unique(.eff$study.num), n.centre=sapply(unique(.eff$study.num), FUN=".temp.func"))

.eff <- merge(.eff, .centre, by.x="study.num", by.y="study.num")

.eff$poids <- .eff$n/sum(.eff$n.centre[.eff$classe==1])

.eff$proba <- .eff$n/.eff$n.centre

nlme1 <- suppressWarnings(nlme(proba ~ d.int(mu, sigma, marker.min, marker.max),
  fixed = mu + sigma ~ 1, random =  mu ~ 1|study.num,
  start = init.nlme1, data = .eff))

nlme1 <- suppressWarnings(update(nlme1, weights=varFixed(~n), data=.eff))

assign("nlme1", nlme1, envir = parent.frame())

.data$mu<-nlme1$coefficients$fixed[1] +
    nlme1$coefficients$random$study.num[.data$study.num]
.data$sigma<-nlme1$coefficients$fixed[2]
.data$poids <- .data$n / sum(.data$n)

.data$beta1.1.fix<-0
.data$beta1.2.fix<-0
.data$beta2.1.fix<-0
.data$beta2.2.fix<-0
.data$beta3.1.fix<-0
.data$beta3.2.fix<-0
.data$beta4.1.fix<-0
.data$beta4.2.fix<-0

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model5(beta0.1, beta0.2, beta1.1.fix, beta1.2.fix, beta2.1.fix, beta2.2.fix, beta3.1.fix, beta3.2.fix, beta4.1.fix, beta4.2.fix, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 ~ 1, random = beta0.1 ~ 1|study.num, start = c(0, 0), data = .data))

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model5(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1.fix, beta2.2.fix, beta3.1.fix, beta3.2.fix, beta4.1.fix, beta4.2.fix, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 + beta1.1 + beta1.2 ~ 1, random = beta0.1 ~ 1|study.num, start = c(nlme2$coefficients$fixed, 0, 0), data = .data))

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model5(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1.fix, beta3.2.fix, beta4.1.fix, beta4.2.fix, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 + beta1.1 + beta1.2 + beta2.1 + beta2.2~ 1, random = beta0.1 ~ 1|study.num, start = c(nlme2$coefficients$fixed, 0, 0), data = .data))

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model5(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1.fix, beta4.2.fix, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 + beta1.1 + beta1.2 + beta2.1 + beta2.2 + beta3.1 + beta3.2 ~ 1, random = beta0.1 ~ 1|study.num, start = c(nlme2$coefficients$fixed, 0, 0), data = .data))

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model5(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, beta4.1, beta4.2, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 + beta1.1 + beta1.2 + beta2.1 + beta2.2 + beta3.1 + beta3.2 + beta4.1 + beta4.2 ~ 1, random = beta0.1 ~ 1|study.num, start = c(nlme2$coefficients$fixed, 0, 0), weights=varPower(form= ~nrisk), data = .data))

assign("nlme2", nlme2, envir = parent.frame())

.roc<-data.frame( cut=seq(
qnorm(0.01, mean = nlme1$coefficients$fixed[1], sd = exp(nlme1$coefficients$fixed[2]), lower.tail = TRUE, log.p = FALSE), 
qnorm(0.99, mean = nlme1$coefficients$fixed[1], sd = exp(nlme1$coefficients$fixed[2]),
lower.tail = TRUE, log.p = FALSE), 
 length.out=precision) )

.se.sp<-apply(.roc, 1, se.sp.marg5, tps=pro.time)

.auc <- auc(sens=.se.sp[1,], spec=.se.sp[2,])

.eff <- .eff[, c("classe", "n", "proba", "marker.min", "marker.max", "study.num")]
.eff$fitted <- nlme1$fitted[,2]
.eff$resid <- .eff$proba - .eff$fitted 

.data <- .data[,c("classe", "year", "surv", "n", "proba", "marker.min", "marker.max", "study.num", "p.joint")]
.data$fitted <- nlme2$fitted[,2]
.data$resid <- .data$p.joint - .data$fitted

return(list(
nlme1 = nlme1,
nlme2 = nlme2,
table = data.frame(cut.off=.roc$cut, se=.se.sp[1,], sp=.se.sp[2,]),
auc = .auc,
data.marker = .eff,
data.surv = .data ))

rm(nlme2, nlme1, se.sp.marg5, p.joint.model5, d.surv5, d.int)
}

if(length(time.cutoff)==3)
{
assign(".s1", time.cutoff[1], envir = parent.frame())
assign(".s2", time.cutoff[2], envir = parent.frame())
assign(".s3", time.cutoff[3], envir = parent.frame())

assign("d.int", function(mu, sigma, min1, max1)
{
noeuds<-gauss.quad(30, kind="legendre", alpha=0, beta=0)$nodes
poids<-gauss.quad(30, kind="legendre", alpha=0, beta=0)$weights

.expr1<-(max1-min1)
.expr2<-(max1+min1)
return(
    0.5*(.expr1)*(
    poids[1]*dnorm(0.5*.expr1*noeuds[1]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[2]*dnorm(0.5*.expr1*noeuds[2]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[3]*dnorm(0.5*.expr1*noeuds[3]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[4]*dnorm(0.5*.expr1*noeuds[4]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[5]*dnorm(0.5*.expr1*noeuds[5]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[6]*dnorm(0.5*.expr1*noeuds[6]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[7]*dnorm(0.5*.expr1*noeuds[7]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[8]*dnorm(0.5*.expr1*noeuds[8]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[9]*dnorm(0.5*.expr1*noeuds[9]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[10]*dnorm(0.5*.expr1*noeuds[10]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[11]*dnorm(0.5*.expr1*noeuds[11]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[12]*dnorm(0.5*.expr1*noeuds[12]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[13]*dnorm(0.5*.expr1*noeuds[13]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[14]*dnorm(0.5*.expr1*noeuds[14]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[15]*dnorm(0.5*.expr1*noeuds[15]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[16]*dnorm(0.5*.expr1*noeuds[16]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[17]*dnorm(0.5*.expr1*noeuds[17]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[18]*dnorm(0.5*.expr1*noeuds[18]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[19]*dnorm(0.5*.expr1*noeuds[19]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[20]*dnorm(0.5*.expr1*noeuds[20]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[21]*dnorm(0.5*.expr1*noeuds[21]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[22]*dnorm(0.5*.expr1*noeuds[22]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[23]*dnorm(0.5*.expr1*noeuds[23]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[24]*dnorm(0.5*.expr1*noeuds[24]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[25]*dnorm(0.5*.expr1*noeuds[25]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[26]*dnorm(0.5*.expr1*noeuds[26]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[27]*dnorm(0.5*.expr1*noeuds[27]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[28]*dnorm(0.5*.expr1*noeuds[28]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
    poids[29]*dnorm(0.5*.expr1*noeuds[29]+0.5*.expr2, mean=mu, sd=exp(sigma))+ 
	poids[30]*dnorm(0.5*.expr1*noeuds[30]+0.5*.expr2, mean=mu, sd=exp(sigma)) ) )
} , envir = parent.frame()) 

assign("d.surv4", function(z, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)
{

.l1 <- .s1
.l2 <- .s2 - .s1
.l3 <- .s3 - .s2

.expr0<-beta0.1+beta0.2*z
.expr1<-beta1.1+beta1.2*z
.expr2<-beta2.1+beta2.2*z
.expr3<-beta3.1+beta3.2*z

.expr5<-(-t*exp(.expr0))*(t<=.s1) +
(-.l1*exp(.expr0) -(t-.s1)*exp(.expr0+.expr1))*(t>.s1)*(t<=.s2) +
(-.l1*exp(.expr0) -.l2*exp(.expr0+.expr1) -(t-.s2)*exp(.expr0+.expr1+.expr2))*(t>.s2)*(t<=.s3) +
(-.l1*exp(.expr0) -.l2*exp(.expr0+.expr1) -.l3*exp(.expr0+.expr1+.expr2) - (t-.s3)*exp(.expr0+.expr1+.expr2+.expr3))*(t>.s3)

return(exp(.expr5) * dnorm(z, mean=mu, sd=exp(sigma)))
} , envir = parent.frame())

assign("p.joint.model4", function(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma, t, min.z, max.z)
{

noeuds<-gauss.quad(30, kind="legendre", alpha=0, beta=0)$nodes
poids<-gauss.quad(30, kind="legendre", alpha=0, beta=0)$weights

.exprA<-0.5*(max.z-min.z)
.exprB<-0.5*(max.z+min.z)
return( .exprA*(
 poids[1]*d.surv4(.exprA*noeuds[1]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[2]*d.surv4(.exprA*noeuds[2]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[3]*d.surv4(.exprA*noeuds[3]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[4]*d.surv4(.exprA*noeuds[4]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[5]*d.surv4(.exprA*noeuds[5]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[6]*d.surv4(.exprA*noeuds[6]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[7]*d.surv4(.exprA*noeuds[7]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[8]*d.surv4(.exprA*noeuds[8]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[9]*d.surv4(.exprA*noeuds[9]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[10]*d.surv4(.exprA*noeuds[10]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[11]*d.surv4(.exprA*noeuds[11]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[12]*d.surv4(.exprA*noeuds[12]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[13]*d.surv4(.exprA*noeuds[13]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[14]*d.surv4(.exprA*noeuds[14]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[15]*d.surv4(.exprA*noeuds[15]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[16]*d.surv4(.exprA*noeuds[16]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[17]*d.surv4(.exprA*noeuds[17]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[18]*d.surv4(.exprA*noeuds[18]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[19]*d.surv4(.exprA*noeuds[19]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[20]*d.surv4(.exprA*noeuds[20]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[21]*d.surv4(.exprA*noeuds[21]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[22]*d.surv4(.exprA*noeuds[22]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[23]*d.surv4(.exprA*noeuds[23]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[24]*d.surv4(.exprA*noeuds[24]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[25]*d.surv4(.exprA*noeuds[25]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[26]*d.surv4(.exprA*noeuds[26]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[27]*d.surv4(.exprA*noeuds[27]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[28]*d.surv4(.exprA*noeuds[28]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[29]*d.surv4(.exprA*noeuds[29]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma)+ 
 poids[30]*d.surv4(.exprA*noeuds[30]+.exprB, t, beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma) ) )
} , envir = parent.frame()) 

assign("se.sp.marg4", function(cut, tps)
{
.marker.min<-qnorm(0.000001, mean = nlme1$coefficients$fixed[1], sd = exp(nlme1$coefficients$fixed[2]), lower.tail = TRUE, log.p = FALSE)

.marker.max<-qnorm(0.999999, mean = nlme1$coefficients$fixed[1], sd = exp(nlme1$coefficients$fixed[2]), lower.tail = TRUE, log.p = FALSE)

est.mu    <- nlme1$coefficients$fixed[1]
est.sigma <- nlme1$coefficients$fixed[2]
est.beta0.1 <- nlme2$coefficients$fixed[1]
est.beta0.2 <- nlme2$coefficients$fixed[2]
est.beta1.1 <- nlme2$coefficients$fixed[3]
est.beta1.2 <- nlme2$coefficients$fixed[4]
est.beta2.1 <- nlme2$coefficients$fixed[5]
est.beta2.2 <- nlme2$coefficients$fixed[6]
est.beta3.1 <- nlme2$coefficients$fixed[7]
est.beta3.2 <- nlme2$coefficients$fixed[8]

exprA.<-d.int(mu=est.mu, sigma=est.sigma, cut, .marker.max)

exprB.<-p.joint.model4(est.beta0.1, est.beta0.2, est.beta1.1, est.beta1.2, est.beta2.1, est.beta2.2, est.beta3.1, est.beta3.2, est.mu, est.sigma, tps, cut, .marker.max)

exprC.<-p.joint.model4(est.beta0.1, est.beta0.2, est.beta1.1, est.beta1.2, est.beta2.1, est.beta2.2, est.beta3.1, est.beta3.2, est.mu, est.sigma, tps, .marker.min, .marker.max)

return( c(( exprA. - exprB.  ) / ( 1 - exprC. ) , 1- exprB. / exprC.) )
} , envir = parent.frame())

.data <- data.frame(classe, n, year, surv, nrisk, proba, marker.min, marker.max, study.num)

.data$p.joint <- .data$surv * .data$proba

.eff <- .data[(.data$year==1), c("classe", "n", "marker.min", "marker.max", "study.num")]

.temp.func <- function(x) { sum(.eff$n[.eff$study.num==x]) }

.centre <- data.frame(study.num = unique(.eff$study.num), n.centre=sapply(unique(.eff$study.num), FUN=".temp.func"))

.eff <- merge(.eff, .centre, by.x="study.num", by.y="study.num")

.eff$poids <- .eff$n/sum(.eff$n.centre[.eff$classe==1])

.eff$proba <- .eff$n/.eff$n.centre

nlme1 <- suppressWarnings(nlme(proba ~ d.int(mu, sigma, marker.min, marker.max),
  fixed = mu + sigma ~ 1, random =  mu ~ 1|study.num,
  start = init.nlme1, data = .eff))

nlme1 <- suppressWarnings(update(nlme1, weights=varFixed(~n), data=.eff))

assign("nlme1", nlme1, envir = parent.frame()) 

.data$mu<-nlme1$coefficients$fixed[1] +
    nlme1$coefficients$random$study.num[.data$study.num]
.data$sigma<-nlme1$coefficients$fixed[2]
.data$poids <- .data$n / sum(.data$n)

.data$beta1.1.fix<-0
.data$beta1.2.fix<-0
.data$beta2.1.fix<-0
.data$beta2.2.fix<-0
.data$beta3.1.fix<-0
.data$beta3.2.fix<-0

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model4(beta0.1, beta0.2, beta1.1.fix, beta1.2.fix, beta2.1.fix, beta2.2.fix, beta3.1.fix, beta3.2.fix, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 ~ 1, random = beta0.1 ~ 1|study.num, start = c(0, 0), data = .data))

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model4(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1.fix, beta2.2.fix, beta3.1.fix, beta3.2.fix, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 + beta1.1 + beta1.2 ~ 1, random = beta0.1 ~ 1|study.num, start = c(nlme2$coefficients$fixed, 0, 0), data = .data))

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model4(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1.fix, beta3.2.fix, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 + beta1.1 + beta1.2 + beta2.1 + beta2.2~ 1, random = beta0.1 ~ 1|study.num, start = c(nlme2$coefficients$fixed, 0, 0), data = .data))

nlme2<-suppressWarnings(nlme(p.joint ~ p.joint.model4(beta0.1, beta0.2, beta1.1, beta1.2, beta2.1, beta2.2, beta3.1, beta3.2, mu, sigma, year, marker.min, marker.max), fixed = beta0.1 + beta0.2 + beta1.1 + beta1.2 + beta2.1 + beta2.2 + beta3.1 + beta3.2~ 1, random = beta0.1 ~ 1|study.num, start = c(nlme2$coefficients$fixed, 0, 0), weights=varPower(form= ~nrisk), data = .data))

assign("nlme2", nlme2, envir = parent.frame())

.roc<-data.frame( cut=seq(
qnorm(0.01, mean = nlme1$coefficients$fixed[1], sd = exp(nlme1$coefficients$fixed[2]), lower.tail = TRUE, log.p = FALSE), 
qnorm(0.99, mean = nlme1$coefficients$fixed[1], sd = exp(nlme1$coefficients$fixed[2]),
lower.tail = TRUE, log.p = FALSE), 
 length.out=precision) )

.se.sp<-apply(.roc, 1, se.sp.marg4, tps=pro.time)

.auc <- auc(sens=.se.sp[1,], spec=.se.sp[2,])

.eff <- .eff[, c("classe", "n", "proba", "marker.min", "marker.max", "study.num")]
.eff$fitted <- nlme1$fitted[,2]
.eff$resid <- .eff$proba - .eff$fitted 

.data <- .data[,c("classe", "year", "surv", "n", "proba", "marker.min", "marker.max", "study.num", "p.joint")]
.data$fitted <- nlme2$fitted[,2]
.data$resid <- .data$p.joint - .data$fitted

return(list(
nlme1 = nlme1,
nlme2 = nlme2,
table = data.frame(cut.off=.roc$cut, se=.se.sp[1,], sp=.se.sp[2,]),
auc = .auc,
data.marker = .eff,
data.surv = .data ))

rm(nlme2, nlme1, se.sp.marg4, p.joint.model4, d.surv4, d.int)
}

}


 