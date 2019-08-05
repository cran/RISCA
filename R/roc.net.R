
roc.net <- function(times, failures, variable, p.age, p.sex, p.year, rate.table, pro.time, cut.off, knn=FALSE, prop=NULL) {

.warnings <- "none"

.temp.data <- data.frame(times=times, failures=failures, variable=variable, age=p.age, sex=p.sex, year=p.year)

.n.na<-sum(is.na(.temp.data$times + .temp.data$failures + .temp.data$variable +
   .temp.data$age + 1*(.temp.data$sex=="male") + .temp.data$year))

.temp.data <- .temp.data[!is.na(.temp.data$times + .temp.data$failures + .temp.data$variable + .temp.data$age + 1*(.temp.data$sex=="male") + .temp.data$year),]
   
.rs.model <- try(summary(relsurv::rs.surv(Surv(times, failures)~ 1 + ratetable(age=age, sex=sex, year=year),  ratetable=rate.table, data=.temp.data, method="pohar-perme")), silent = TRUE)

  if(inherits(.rs.model, "try-error")) {return(NA)}
  
  else {
  
  if(sum(.rs.model$surv>1)>0) { .warnings <- "Esimations of the net survival are higher than 1." }
  
  
if (knn==FALSE) {

.surv.total <- data.frame(time = c(1-0.99999999999, .rs.model$time),
     surv = c(0.99999999999, .rs.model$surv))
.surv.total$surv[.surv.total$surv==0] <- 1-0.99999999999
.surv.total$indic <- (.surv.total$time<=pro.time)
  
.temp.se<-function(x) {

  .rs.model.cond <- try(summary(relsurv::rs.surv(Surv(times, failures)~ 1 + ratetable(age=age, sex=sex, year=year),  ratetable=rate.table, data=.temp.data[.temp.data$variable>x,], method="pohar-perme")), silent = TRUE)
  
  if(inherits(.rs.model.cond, "try-error")) {return(NA)}

  else {
  
  .surv.cond <- data.frame(time = c(1-0.99999999999, .rs.model.cond$time),
     surv = c(0.99999999999, .rs.model.cond$surv))
  .surv.cond$surv[.surv.cond$surv==0] <- 1-0.99999999999
  .surv.cond$indic <- (.surv.cond$time<=pro.time)

   return( (1-.surv.cond$surv[.surv.cond$indic][sum(.surv.cond$indic)]) * 
           (1-(sum(.temp.data$variable<=x)/length(times))) / 
           (1-.surv.total$surv[.surv.total$indic][sum(.surv.total$indic)]) ) } }
		   
.temp.sp<-function(x) {
  
  .rs.model.cond <- try(summary(relsurv::rs.surv(Surv(times, failures)~ 1 + ratetable(age=age, sex=sex, year=year),  ratetable=rate.table, data=.temp.data[.temp.data$variable<=x,], method="pohar-perme")), silent = TRUE)

  if(inherits(.rs.model.cond, "try-error")) {return(NA)}

  else {
  
  .surv.cond <- data.frame(time = c(1-0.99999999999, .rs.model.cond$time),
     surv = c(0.99999999999, .rs.model.cond$surv))
  .surv.cond$surv[.surv.cond$surv==0] <- 1-0.99999999999
  .surv.cond$indic <- (.surv.cond$time<=pro.time)

    return( (.surv.cond$surv[.surv.cond$indic][sum(.surv.cond$indic)]) *
	        (sum(.temp.data$variable<=x)/length(times)) /
            .surv.total$surv[.surv.total$indic][sum(.surv.total$indic)] ) } }

.tab.res<-data.frame(
   cut.off=cut.off,
   se=suppressWarnings(sapply(cut.off, FUN=".temp.se")),
   sp1=1-suppressWarnings(sapply(cut.off, FUN=".temp.sp")))

.tab.res$se[.tab.res$se>1]<-NA
.tab.res$sp1[.tab.res$sp1>1]<-NA

.tab.res<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se),]

.tab.res.temp<-rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1))

.tab.res$sp<-1-.tab.res$sp1

return(list(
table = .tab.res[,c("cut.off", "se", "sp")],
auc = sum((.tab.res.temp$sp1[2:length(.tab.res.temp$sp1)]-.tab.res.temp$sp1[1:(length(.tab.res.temp$sp1)-1)])*(.tab.res.temp$se[2:length(.tab.res.temp$se)])),
missing = .n.na  ) ) }


if (knn==TRUE) {

.temp.data <- data.frame(times=times, failures=failures, variable=variable,
   age = p.age, sex = p.sex, year = p.year)

.n.na<-sum(is.na(.temp.data$times + .temp.data$failures + .temp.data$variable +
   .temp.data$age + 1*(.temp.data$sex=="male") + .temp.data$year))
   
.temp.data <- .temp.data[!is.na(.temp.data$times + .temp.data$failures + .temp.data$variable + .temp.data$age + 1*(.temp.data$sex=="male") + .temp.data$year),]
   
.temp.data <- .temp.data[order(.temp.data$variable),]

.pas <- round(prop * dim(.temp.data)[1])

.n <- dim(.temp.data)[1]

.survie.x<-function(x)  {

   .tmp<-try(summary(relsurv::rs.surv(Surv(times, failures)~ 1 + ratetable(age=age, sex=sex, year=year),  ratetable=rate.table, data=.temp.data[max(1, (x-.pas)):min((x+.pas), .n),], method="pohar-perme")), silent = TRUE)

  if(inherits(.tmp, "try-error")) {return(NA)}

  else {
  
  .tmp <- data.frame(time = c(1-0.99999999999, .tmp$time),
     surv = c(0.99999999999, .tmp$surv))
  .tmp$surv[.tmp$surv==0] <- 1-0.99999999999
  .tmp$indic <- (.tmp$time<=pro.time)

   return(.tmp$surv[.tmp$indic][sum(.tmp$indic)]) } }

.survie.temp<-sapply(1:.n, FUN = ".survie.x")

.survie.prop<-function(x) {
 mean(.survie.temp * (.temp.data$variable>x), na.rm=TRUE) }

.survie.marginale<-mean(.survie.temp, na.rm=TRUE)

.temp.se<-function(x) {
 (1-sum(.temp.data$variable<=x)/.n-.survie.prop(x))/(1-.survie.marginale) }
 
.temp.sp<-function(x)
{ 1-(.survie.prop(x)/.survie.marginale) }

.tab.res<-data.frame(
   cut.off=cut.off,
   se=suppressWarnings(sapply(cut.off, FUN=".temp.se")),
   sp1=1-suppressWarnings(sapply(cut.off, FUN=".temp.sp")))

.tab.res$se[.tab.res$se>1]<-NA
.tab.res$sp1[.tab.res$sp1>1]<-NA

.tab.res<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res[!is.na(.tab.res$sp1 + .tab.res$se),]

.tab.res.temp<-.tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se),]

.tab.res.temp<-rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1))

.tab.res$sp<-1-.tab.res$sp1



return(list(
table = .tab.res[,c("cut.off", "se", "sp")],
auc = sum((.tab.res.temp$sp1[2:length(.tab.res.temp$sp1)] - .tab.res.temp$sp1[1:(length(.tab.res.temp$sp1)-1)]) * 0.5 * (.tab.res.temp$se[2:length(.tab.res.temp$sp1)] + .tab.res.temp$se[1:(length(.tab.res.temp$sp1)-1)])),
missing = .n.na,
warning = .warnings ) )

}
}
}
