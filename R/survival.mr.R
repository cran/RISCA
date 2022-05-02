
survival.mr <- function(times, failures, cov.relative, data, cox.reference, cov.reference,  ini, iterations)
{

coef.simul <- data.frame(rmvnorm(iterations, mean=cox.reference$coefficient, sigma=cox.reference$var))
rel.data <- data
n <- dim(rel.data)[1]
rel.data$ident <- 1:n
res.theta <- matrix(-99, nrow=iterations, ncol=length(cov.relative))
dimnames(res.theta)[[2]] <- cov.relative

for (b in 1:iterations)
{
num<-sample(rel.data$ident, size=n, replace = TRUE)
data.boot<-rel.data[num,]
data.boot<-data.boot[order(data.boot[,times]),]

data.boot$nb<-1
for (i in unique(data.boot[data.boot[,failures], times])) {
data.boot$nb[data.boot[,times]==i & data.boot[,failures]==1] <- 1:sum(data.boot[,times]==i & data.boot[,failures]==1) }
data.boot$Tps.Evt.Cor<-data.boot[,times]+data.boot$nb/1000
data.boot<-data.boot[order(data.boot$Tps.Evt.Cor, 1-data.boot[,failures]),]
data.boot$indic<-1:n

covarAtt<-data.boot[,cov.reference]
covarRel<-data.boot[,cov.relative]
a <- as.matrix(covarAtt) %*% as.vector(as.numeric(coef.simul[b,]))

logVP.boot <-     
    
model <- optim(par=ini,
             fn = function(x){
                return(sum(sapply(data.boot$indic[data.boot[,failures]==1], FUN = function(y) {
                a[y] + as.matrix(covarRel[y,]) %*% as.vector(x[1:(dim(covarRel)[2])]) -
                log(sum(exp(a[y:n]) * exp(as.matrix(covarRel[y:n,]) %*% as.vector(x[1:(dim(covarRel)[2])])))) } ) ) ) },
    method = "Nelder-Mead", hessian = TRUE,  control=list(fnscale=-1, maxit=100000) )

res.theta[b,]<-model$par
}

return(list(
matrix.coef = res.theta,
estim.coef = apply(res.theta, FUN="mean", MARGIN=2),
lower95.coef = apply(res.theta, FUN= function(x) {quantile(x, probs=0.025)}, MARGIN=2),
upper95.coef = apply(res.theta, FUN= function(x) {quantile(x, probs=0.975)}, MARGIN=2)
))

}