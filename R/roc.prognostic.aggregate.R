
roc.prognostic.aggregate<-function(time.lr, surv.lr, time.hr, surv.hr)
{
.time<-sort(unique(c(time.lr, time.hr)))
.time<-.time[.time<=min(max(time.lr), max(time.hr))]
.n<-length(.time)
.s.LR<-rep(1, .n)
.s.HR<-rep(1, .n)

for (i in .time) {
.s.LR[.time==i]<-min(1, surv.lr[time.lr<=i][sum(time.lr<=i)])
.s.HR[.time==i]<-min(1, surv.hr[time.hr<=i][sum(time.hr<=i)])
}

.table<-data.frame(time=.time, x=1-.s.LR, y=1-.s.HR)

.auc <- sum((.table$x[2:.n] - .table$x[1:(.n-1)]) * 0.5 * (.table$y[2:.n] + .table$y[1:(.n-1)]))

.s.HR.last <- .s.HR[.n]
.s.LR.last <- .s.LR[.n]

.lower <- .auc + .s.LR.last * (1-.s.HR.last)

.pessimist <- .lower +  (.s.HR.last) * (.s.HR.last) * 0.5

.upper <- .auc + .s.LR.last

.non.informative <- .lower +  .s.LR.last * .s.HR.last * 0.5

.optimist <- .lower / (1- .s.LR.last * .s.HR.last) 

return(list(max.time=min(max(time.lr), max(time.hr)),
 table=.table,
 auc=data.frame(
  lower=.lower,
  pessimist=.pessimist,
  noninformative=.non.informative,
  optimist=.optimist,
  upper=.upper ) ) )
}

