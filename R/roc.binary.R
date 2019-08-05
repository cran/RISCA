
roc.binary <- function(status, variable, confounders, data, precision=seq(0.05, 0.95, by=0.025), estimator="ipw")
{

cut.off <- quantile(data[,variable], probs=precision, na.rm=TRUE)

if((max(precision)==1) | (min(precision)==0)){ stop("The cut-off values have to be different from the minimum or the maximum of the variable") }

if(estimator != "ipw" & estimator != "pv"){stop("Error: the argument used in estimator is incorrect")}

data$temp <- data[,status] + data[,variable]

form0 <- update.formula(confounders, temp ~ .)

if((length(data[,status]) - summary(glm(form0, data=data))$df.null - 1) > 0) {stop("Error: missing values are not allowed")}

if(estimator == "ipw"){

data$YyY54 <- data[, status]
form <- update.formula(confounders, YyY54 ~ .)
W <- glm(form, family = binomial(link = logit), data = data)$fitted.values

se <- function(x) {
	sum(1 * (data[, variable] > cut.off[x]) * data[, 
		status] * pmin(1/(W),1/(mean(W)^2*0.1)))/sum(data[, status] * 
		pmin(1/(W),1/(mean(W)^2*0.1)))
}
sp <- function(x) {
	sum(1 * (data[, variable] <= cut.off[x]) * (data[, 
		status] - 1) * pmin(1/(1-W),1/(mean(1-W)^2*0.1)))/sum((data[, status] - 
		1) * pmin(1/(1-W),1/(mean(1-W)^2*0.1)))
}

temp.se <- sapply(1:length(cut.off), FUN = "se")
temp.sp <- sapply(1:length(cut.off), FUN = "sp")

.tab <-data.frame( cut.off = cut.off, se = temp.se, sp1 = 1-temp.sp)

.tab$se[.tab$se > 1] <- NA
.tab$sp1[.tab$sp1 > 1] <- NA

.tab.res.temp <- .tab[!is.na(.tab$sp1 + .tab$se), ]
.tab.res.temp <- .tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se), ]
.tab.res.temp <- rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1)) 
colnames(.tab.res.temp) <- colnames(.tab)

.tab$sp <- 1 - .tab$sp1
.tab <- rbind(c(min(data[,variable]),1,1,0) ,.tab, c(max(data[,variable]),0,0,1))

if(dim(.tab.res.temp)[1]>2){
.auc <- auc(.tab.res.temp$se, 1-.tab.res.temp$sp1)
}else{.auc<-NA}

.obj <- list(table=.tab[,c("cut.off", "se", "sp")], auc = .auc)

class(.obj) <- "roc"

return(.obj)
}


if(estimator == "pv")
{
data$YyY54 <- data[, variable]
data0 <- data[data[, status] == 0, ]
form <- update.formula(confounders, YyY54 ~ .)
lm.0 <- lm(form, data = data0)

chaine <- "(data[,variable] - (lm.0$coef[1]"
if( length(names((lm.0)$coef)) > 1){
    for (i in 2:length(names((lm.0)$coef))){
    chaine <- paste(chaine, " + lm.0$coef[",i,"]*data[,names((lm.0)$coef[",i,"])]",sep="") 
    }
}

chaine <- paste(chaine,"))/sd(lm.0$residuals)",sep="")
.y0 <- eval(parse(text=chaine))

.pv <- sapply(.y0[data[,status]==1], FUN = function(x) {sum(.y0[data[,status]==0]<=x)/length(.y0[data[,status]==0])})

.ROCu <- sapply(precision, FUN=function(x) {sum((1-.pv)<=x)/length(.pv)})

.tab <-data.frame(se = .ROCu, sp = 1-precision)

.tab$se[.tab$se > 1] <- NA
.tab$sp[.tab$sp > 1] <- NA

.tab <- .tab[order(.tab$sp, .tab$se), ]
.tab <- rbind(c(1,0) ,.tab, c(0,1))

.auc <- auc(.tab$se, .tab$sp)

.obj <- list(table=.tab[,c("se", "sp")], auc = .auc)

class(.obj) <- "roc"

return(.obj)
}

}
