
gc.sl.time <- function(object, max.time=NULL, effect="ATE", estim.tune=TRUE,
        estim.weights=TRUE, methods=NULL, conditional=FALSE,  iterations=100, length.out=100)
{
    # iterations <- 20
    # object <- sl1
    # pro.time <- object$pro.time
    # effect <- "ATE"
    # n.cluster <- 20
    # cluster.type <- "PSOCK"
    # conditional <- TRUE
    # length.out <- 100
  
  if(!is.null(max.time)){ rt <- max.time }
  if(is.null(max.time) & is.null(object$pro.time)) { rt <- median(object$data[,"times"])  }  
  if(is.null(max.time) & !is.null(object$pro.time)) { rt <- object$pro.time } 
  
  if(is.null(object$predictors$group)) { stop("The argument \"group\" in the sl.time function is NULL")  }
  
  if(class(object) != "sl.time") { stop("The argument \"object\" needs to be obtain by the sl.time function") }
  
  if( !is.null(methods) ) { if(min(methods %in% object$methods) == 0) stop("Only methods included in the sl.time object can be asked in the argument \"method\".") }
  
  .methods <- c(methods, "sl")
  
  mod <- unique(object$data[,object$predictors$group])
  if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
    stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and
    1 (for treated/exposed patients) are required in the argument \"group\"") }
  
  mod <- unique(object$data[,"failures"])
  if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2]!= 0))){
    stop("Two modalities encoded 0 (for censored patients) and 1 (for
    uncensored patients) are required in the argument \"failures\" ") }
  
  if(is.na(match(effect,c("ATT","ATU","ATE"))) | length(effect)!=1){
    stop("Incorrect modality specified in the argument \"effect\": only one
    option among \"ATE\", \"ATT\" or \"ATU\" are possible") }
  
  if(rt > max(object$data[,"times"], na.rm=T) | rt < min(object$data[,"times"], na.rm=T)){
    stop("The argument \"max.time\" needs to be a number includes between the
    min and the max of the vector \"times\" ") }
  
  data1 <- data0 <- object$data
  data1[,object$predictors$group] <- 1; data0[,object$predictors$group] <- 0
  
  #ti <- sort(unique(object$data[object$data[,"failures"]==1, "times"]))
  
  times <- seq(0, max(object$times), length.out=length.out)
  
  s1 <- s0 <- logHR <- RMST1 <- RMST0 <- delta <- list() 
  
  for(i in .methods){
    
    surv.sl1 <- predict(object, newdata=data1, newtimes = times)$predictions[[i]]
    surv.sl0 <- predict(object, newdata=data0, newtimes = times)$predictions[[i]]
    
    s1[[i]] <- colMeans(surv.sl1)
    s0[[i]] <- colMeans(surv.sl0)
    
    h1 <- differentiation(x=times, fx=-1*log(s1[[i]]))
    h0 <- differentiation(x=times, fx=-1*log(s0[[i]]))
    
    .h1 <- spliner(h1 ~ times, monotonic = FALSE)
    .h0 <- spliner(h0 ~ times, monotonic = FALSE)
    
    .h <- function(x) { .h1(x) + .h0(x) }
    
    .f <- function(x) {
      .h(x) * exp(-1*cubintegrate(f=.h, lower=0, upper=x,  method="pcubature")$integral) }
    
    .ah1 <- function(x) { (.h1(x)/.h(x)) * .f(x) }
    .ah0 <- function(x) { (.h0(x)/.h(x)) * .f(x) }
    
    logHR[[i]] <- log( cubintegrate(f =.ah1, lower=0, upper=max(object$times),
                method="pcubature")$integral / cubintegrate(f =.ah0, lower=0,
               upper=max(object$times),  method="pcubature")$integral )
    
    RMST1[[i]] <- rmst(times, s1[[i]], max.time=rt, type="l")
    RMST0[[i]] <- rmst(times, s0[[i]], max.time=rt, type="l")
    delta[[i]] <- RMST1[[i]] - RMST0[[i]]
  
  }
    
  #logHR <- log(mean(D1$h/D0$h))
  
  # if(conditional)
  #   {
  #   h.sl1 <- h.sl0 <- surv.sl1; logHRcond.i <- rep(-99, dim(surv.sl0)[1])
  #   
  #   for (j in 1:(dim(surv.sl0)[1]))  {
  #     h.sl1[j,] <- differentiation(x=ti, fx=-1*log(surv.sl1[j,]))
  #     h.sl0[j,] <- differentiation(x=ti, fx=-1*log(surv.sl0[j,]))
  #     logHRcond.i[j] <- mean(h.sl1[j,]/h.sl0[j,])  }
  #   
  #   logHRcond <- log(mean(logHRcond.i, na.rm=TRUE))
  #   }
  
  #ti <- object$times
  
  #surv.sl1 <- predict(object, newdata=data1)$predictions$sl
  #surv.sl0 <- predict(object, newdata=data0)$predictions$sl
  
  #D1 <- data.frame(S = colMeans(surv.sl1), T = ti)
  #D0 <- data.frame(S = colMeans(surv.sl0), T = ti)

  # S1.b <- S0.b <- matrix(-99, nrow=iterations, ncol=length(times)) 
  # logHR.b <- RMST1.b <- RMST0.b <- delta.b <- rep(-99, iterations) 
    
  #if(conditional){logHRcond.b  <- rep(-99, iterations)}
    
  S1.b <- S0.b <- logHR.b <- RMST1.b <- RMST0.b <- delta.b <- list() 
  
  for(j in .methods){
    S1.b[[j]] <- S0.b[[j]] <- matrix(-99, nrow=iterations, ncol=length(times)) 
    logHR.b[[j]] <- RMST1.b[[j]] <- RMST0.b[[j]] <- delta.b[[j]] <- rep(-99, iterations) 
  }
  
  for(i in 1:iterations){

      id <- sample(1:(dim(object$data)[1]), size = dim(object$data)[1], replace = TRUE)
      learn <- object$data[id,]
      
      sl.b<-sl.time(
        methods=object$methods, metric=object$metric,
        data=learn,
        times="times",  failures="failures",
        pro.time = object$pro.time,
        group=object$predictors$group,
        cov.quanti=object$predictors$cov.quanti,  cov.quali=object$predictors$cov.quali,
        cv=object$cv,
        param.tune = if(estim.tune==FALSE){object$param.tune$optimal}else{NULL},
        param.weights.fix=if(estim.weights==FALSE){object$weights$coefficients}else{NULL},
        keep.predictions=TRUE,
        verbose=FALSE)
      
      #learn1 <- learn0 <- learn
      #learn1[,object$predictors$group] <- 1; learn0[,object$predictors$group] <- 0
      #surv.sl1.b <- predict(sl.b, newdata=learn1, newtimes=times)$predictions$sl
      #surv.sl0.b <- predict(sl.b, newdata=learn0, newtimes=times)$predictions$sl
      
      # ti.b <- sort(unique(learn[learn[,"failures"]==1, "times"]))
      
      valid1 <- valid0 <- valid <- object$data[-sort(unique(id)),]
      valid1[,object$predictors$group] <- 1; valid0[,object$predictors$group] <- 0
      
      for(j in .methods){
      
        surv.sl1.b <- predict(sl.b, newdata=valid1, newtimes=times)$predictions[[j]]
        surv.sl0.b <- predict(sl.b, newdata=valid0, newtimes=times)$predictions[[j]]
        
        s1.b <- colMeans(surv.sl1.b)
        s0.b <- colMeans(surv.sl0.b)
        
        h1.b <- differentiation(x=times, fx=-1*log(s1.b))
        h0.b <- differentiation(x=times, fx=-1*log(s0.b))
        
        .h1.b <- spliner(h1.b ~ times, monotonic = FALSE)
        .h0.b <- spliner(h0.b ~ times, monotonic = FALSE)
        .h.b <- function(x) { .h1.b(x) + .h0.b(x) }
        .f.b <- function(x) {
          .h.b(x) * exp( -1 * cubintegrate(f = .h.b, lower = 0, upper = x,
                                         method = "pcubature")$integral) }
        .ah1.b <- function(x) { (.h1.b(x)/.h.b(x)) * .f.b(x) }
        .ah0.b <- function(x) { (.h0.b(x)/.h.b(x)) * .f.b(x) }
        
        logHR.b[[j]][i] <- log( cubintegrate(f =.ah1.b, lower=0, upper=max(object$times),
                                   method="pcubature")$integral / cubintegrate(f =.ah0.b,
                                  lower=0, upper=max(object$times),  method="pcubature")$integral )
      
        #logHR.b[i] <- log(mean(D1.b$h/D0.b$h))
        
        # if(conditional)
        #   {
        #   h.sl1.b <- h.sl0.b <- surv.sl1.b; logHRcond.i.b <- rep(-99, dim(surv.sl0.b)[1])
        #   
        #   for (j in 1:(dim(surv.sl0.b)[1])) {
        #     h.sl1.b[j,] <- differentiation(x=ti.b, fx=-1*log(surv.sl1.b[j,]))
        #     h.sl0.b[j,] <- differentiation(x=ti.b, fx=-1*log(surv.sl0.b[j,]))
        #     logHRcond.i.b[j] <- mean(h.sl1.b[j,]/h.sl0.b[j,])  }
        #   
        #   logHRcond.b[i] <- log(mean(logHRcond.i.b))
        #   }
        
        #ti.b <- object$times
        #surv.sl1.b <- predict(sl.b, newdata=learn1, newtimes = ti.b)$predictions$sl
        #surv.sl0.b <- predict(sl.b, newdata=learn0, newtimes = ti.b)$predictions$sl
        
        #D1.b <- data.frame(S = colMeans(surv.sl1.b), T = ti.b)
        #D0.b <- data.frame(S = colMeans(surv.sl0.b), T = ti.b)
        
        S1.b[[j]][i,] <- s1.b
        S0.b[[j]][i,] <- s0.b
        
        RMST1.b[[j]][i] <- rmst(times, s1.b, max.time=rt, type="l")
        RMST0.b[[j]][i] <- rmst(times, s0.b, max.time=rt, type="l")
        delta.b[[j]][i] <- RMST1.b[[j]][i] - RMST0.b[[j]][i]
      
      }
          
  }
  
  .obj.temp <- list()
  
  for(j in .methods){
    
    se.RMST1 <- sd(RMST1.b[[j]], na.rm=TRUE) 
    se.RMST0 <- sd(RMST0.b[[j]], na.rm=TRUE)
    se.logHR <- sd(logHR.b[[j]], na.rm=TRUE)
    se.delta <- sd(delta.b[[j]], na.rm=TRUE)
    
    pv <- function(m, s){ ztest <- m/s;  return(ifelse(ztest<0,2*pnorm(ztest),2*(1-pnorm(ztest)))) }
    
    p.value.logHR <- pv(m=logHR[[j]], s=se.logHR)
    #if(p.value.logHR==0){ p.value.logHR <- "<0.001" }
    
    ci.low.logHR <- logHR[[j]] - qnorm(0.975, 0, 1)*se.logHR
    ci.upp.logHR <- logHR[[j]] + qnorm(0.975, 0, 1)*se.logHR
    
    ci.low.RMST1 <- RMST1[[j]] - qnorm(0.975, 0, 1)*se.RMST1
    ci.upp.RMST1 <- RMST1[[j]] + qnorm(0.975, 0, 1)*se.RMST1
    
    ci.low.RMST0 <- RMST0[[j]] - qnorm(0.975, 0, 1)*se.RMST0
    ci.upp.RMST0 <- RMST0[[j]] + qnorm(0.975, 0, 1)*se.RMST0
    
    p.value.delta <- pv(m=delta[[j]], s=se.delta)
    #if(p.value.delta==0){ p.value.delta <- "<0.001" }
    
    ci.low.delta <- delta[[j]] - qnorm(0.975, 0, 1)*se.delta
    ci.upp.delta <- delta[[j]] + qnorm(0.975, 0, 1)*se.delta
    
    ci.low.S1 <- apply(S1.b[[j]], FUN = function(x) {quantile(x, probs=0.025, na.rm=TRUE)}, MARGIN = 2)
    ci.upp.S1 <- apply(S1.b[[j]], FUN = function(x) {quantile(x, probs=0.975, na.rm=TRUE)}, MARGIN = 2)
    
    ci.low.S0 <- apply(S0.b[[j]], FUN = function(x) {quantile(x, probs=0.025, na.rm=TRUE)}, MARGIN = 2)
    ci.upp.S0 <- apply(S0.b[[j]], FUN = function(x) {quantile(x, probs=0.975, na.rm=TRUE)}, MARGIN = 2)
    
    table.surv <- data.frame(
        times=times,
        survival=c(s1[[j]], s0[[j]]),
        ci.lower=c(ci.low.S1, ci.low.S0), 
        ci.upper=c(ci.upp.S1, ci.upp.S0), 
        n.risk=c(
          sapply(times, function(x) { sum(object$data$times[object$data[,object$predictors$group]==1] >= x) } ),
          sapply(times, function(x) { sum(object$data$times[object$data[,object$predictors$group]==0] >= x) }) ),
        variable=c( rep(1, length(times)), rep(0,length(times)) ) 
    )
    
    if(conditional==FALSE){
      
      .obj.temp[[j]] <- list( 
        method = j,
        table.surv=table.surv,  
        RMST0=data.frame(estimate=RMST0[[j]], std.error = se.RMST0, ci.lower=ci.low.RMST0, 
                         ci.upper=ci.upp.RMST0),
        RMST1=data.frame(estimate=RMST1[[j]], std.error = se.RMST1, ci.lower=ci.low.RMST1, 
                         ci.upper=ci.upp.RMST1),
        delta=data.frame(estimate=delta[[j]], std.error = se.delta, ci.lower=ci.low.delta, 
                         ci.upper=ci.upp.delta, p.value=p.value.delta),
        logHR=data.frame(estimate=logHR[[j]], std.error = se.logHR, ci.lower=ci.low.logHR, 
                       ci.upper=ci.upp.logHR, p.value=p.value.logHR) 
      ) 
    }
  
  }
  
  .obj<- list()

  .obj$effect <- effect
  .obj$max.time <- rt
  .obj$RMST0 <- .obj.temp[["sl"]]$RMST0
  .obj$RMST1 <- .obj.temp[["sl"]]$RMST1
  .obj$delta <- .obj.temp[["sl"]]$delta
  .obj$logHR <- .obj.temp[["sl"]]$logHR
  .obj$table.surv <- .obj.temp[["sl"]]$table.surv

  if(length(.methods)>1) {.obj$learner <- .obj.temp[names(.obj.temp) != "sl"]}
  #if(length(.methods)>1) {.obj$learner <- .obj.temp}
  
  # if(conditional==TRUE){
  #   
  #   se.logHRcond <- sd(logHRcond.b, na.rm=TRUE)
  #   
  #   pv <- function(m, s){ ztest <- m/s;  return(ifelse(ztest<0,2*pnorm(ztest),2*(1-pnorm(ztest)))) }
  #   
  #   p.value.logHRcond <- pv(m=logHRcond, s=se.logHRcond)
  #   if(p.value.logHRcond==0)  { p.value.logHRcond <- "<0.001" }
  #   
  #   ci.low.logHRcond <- logHRcond - qnorm(0.975, 0, 1)*se.logHRcond
  #   ci.upp.logHRcond <- logHRcond + qnorm(0.975, 0, 1)*se.logHRcond
  # 
  #   .obj <- list(  table.surv=table.surv,  effect=effect, max.time=rt,
  #   RMST0=data.frame(estimate=RMST0, std.error = se.RMST0,
  #                    ci.lower=ci.low.RMST0, ci.upper=ci.upp.RMST0),
  #   RMST1=data.frame(estimate=RMST1, std.error = se.RMST1,
  #                    ci.lower=ci.low.RMST1, ci.upper=ci.upp.RMST1),
  #   delta=data.frame(estimate=delta, std.error = se.delta,
  #                    ci.lower=ci.low.delta, ci.upper=ci.upp.delta, p.value=p.value.delta),
  #   logHR=data.frame(estimate=logHR, std.error = se.logHR,
  #                    ci.lower=ci.low.logHR, ci.upper=ci.upp.logHR, p.value=p.value.logHR),
  #   logHR.conditional = data.frame(estimate=logHRcond, std.error = se.logHRcond,
  #                    ci.lower=ci.low.logHRcond, ci.upper=ci.upp.logHRcond, p.value=p.value.logHRcond),
  #   logHR.conditional.values = logHRcond.i ) }
  
  class(.obj) <- "survrisca"
  return(.obj)
}


# library(mosaic)
# library(mosaicCalc)
# library(cubature)
# 
# data(dataDIVAT2)
# 
# sl1<-sl.time( methods=c("aft.gamma", "ph.gompertz"),  metric="ibs",
#               data=dataDIVAT2,  times="times", failures="failures", group="ecd",
#               cov.quanti=c("age"),  cov.quali=c("hla", "retransplant"), cv=5)
# 
# #Marginal effect of the treatment (ATE): use 1000 iterations instead of 3
# gc.ate <- gc.sl.time(sl1, max.time=12, effect="ATE", iterations=6, length.out=100)
# 
# #Plot the survival curves
# plot(gc.ate, ylab="Confounder-adjusted survival",
#      xlab="Time post-transplantation (years)", col=c(1,2))
