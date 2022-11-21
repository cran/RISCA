
hr.sl.time <- function(object, pro.time=NULL, iterations=1000,
                  n.cluster = 1, cluster.type="PSOCK", conda.env=NULL)
{
 
 #  iterations <- 20
 #  object <- sl3
 #  pro.time <- object$pro.time
 #  n.cluster <- 20
 #  cluster.type <- "PSOCK"

  if(is.null(pro.time) & is.null(object$pro.time)) { pro.time <- median(object$data[,"times"])  }  
  if(is.null(pro.time) & !is.null(object$pro.time)) { pro.time <- object$pro.time } 
  
  if(is.null(object$predictors$group)) { stop("The argument \"group\" in the sl.time function is NULL")  }
  
  if(is(object) != "sl.time") { stop("The argument \"object\" needs to be obtain by the sl.time function") }
  
  mod <- unique(object$data[,object$predictors$group])
  if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
    stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and
    1 (for treated/exposed patients) are required in the argument \"group\"") }
  
  mod <- unique(object$data[,"failures"])
  if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2]!= 0))){
    stop("Two modalities encoded 0 (for censored patients) and 1 (for
    uncensored patients) are required in the argument \"failures\" ") }
  
  if(pro.time > max(object$data[,"times"], na.rm=T) | pro.time < min(object$data[,"times"], na.rm=T)){
    stop("The argument \"pro.time\" needs to be a number includes between the
    min and the max of the vector \"times\" ") }

        data1 <- data0 <- object$data
        data1[,object$predictors$group] <- 1; data0[,object$predictors$group] <- 0

        # pas sur de utilite de BCV par rapport a Bootstrap simple ???
        pred1 <- predict(object, newdata=data1, newtimes=object$times)
        pred0 <- predict(object, newdata=data0, newtimes=object$times)
        
        # faut-il vraiment aller prendre les temps de echantillon apprentissage ???
        surv.sl1 <- pred1$predictions[[length(object$methods)+1]][,pred1$times %in%
         sort(unique(object$data[object$data[,"failures"]==1, "times"]))]
         
        surv.sl0 <- pred0$predictions[[length(object$methods)+1]][,pred0$times %in%
          sort(unique(object$data[object$data[,"failures"]==1, "times"]))]
        
        ti  <- pred1$times[pred1$times %in%
         sort(unique(object$data[object$data[,"failures"]==1, "times"]))]
        
        h.sl1 <- h.sl0 <- surv.sl1
        
        logHRcond.i <- rep(-99, dim(surv.sl0)[1])
        
        for (j in 1:(dim(surv.sl0)[1])) {
        h.sl1[j,] <- differentiation(x=ti, fx=-1*log(surv.sl1[j,]))
        h.sl0[j,] <- differentiation(x=ti, fx=-1*log(surv.sl0[j,]))
        logHRcond.i[j] <- log(mean(h.sl1[j,]/h.sl0[j,]))  }
  
  
  if(n.cluster==1){

    logHRcond.b  <- rep(-99, iterations)

  for(i in 1:iterations){

        id <- sample(1:(dim(object$data)[1]), size = dim(object$data)[1], replace = TRUE)
        learn <- object$data[id,]
        valid <- object$data[-sort(unique(id)),]
        
        sl.b<-sl.time( # YOHANN2: Il faudra une option avec ou sans hyperparametres fixes
          methods=object$methods, metric=object$metric, data=learn, times="times",
          failures="failures", pro.time = pro.time, group=object$predictors$group,
          cov.quanti=object$predictors$cov.quanti,  cov.quali=object$predictors$cov.quali,
          cv=object$cv,  param.tune = object$param.tune$optimal,
          param.weights.fix=object$weights$coefficients, keep.predictions=FALSE, # Il faudra passer a FALSE
          verbose=FALSE)

        valid1 <- valid0 <- valid
        valid1[,object$predictors$group] <- 1; valid0[,object$predictors$group] <- 0

        # pas sur de utilite de BCV par rapport a Bootstrap simple ???
        pred1.b <- predict(sl.b, newdata=valid1, newtimes=object$times)  
        pred0.b <- predict(sl.b, newdata=valid0, newtimes=object$times)
        
        # faut-il vraiment aller prendre les temps de l echantillon d apprentissage ???
        surv.sl1.b <- pred1.b$predictions[[length(object$methods)+1]][,pred1.b$times %in%
         sort(unique(learn[learn[,"failures"]==1, "times"]))]
         
        surv.sl0.b <- pred0.b$predictions[[length(object$methods)+1]][,pred0.b$times %in%
          sort(unique(learn[learn[,"failures"]==1, "times"]))]
        
        ti.b  <- pred1.b$times[pred1.b$times %in%
         sort(unique(learn[learn[,"failures"]==1, "times"]))]
        
        h.sl1.b <- h.sl0.b <- surv.sl1.b
        
        logHRcond.i.b <- rep(-99, dim(surv.sl0.b)[1])
        
        for (j in 1:(dim(surv.sl0.b)[1])) {
        h.sl1.b[j,] <- differentiation(x=ti.b, fx=-1*log(surv.sl1.b[j,]))
        h.sl0.b[j,] <- differentiation(x=ti.b, fx=-1*log(surv.sl0.b[j,]))
        logHRcond.i.b[j] <- log(mean(h.sl1.b[j,]/h.sl0.b[j,]))  }
        
        logHRcond.b[i] <- mean(logHRcond.i.b)
        }   }

  if(n.cluster>1){

        simul <- function() {
          
        id <- sample(1:(dim(object$data)[1]), size = dim(object$data)[1], replace = TRUE)
        learn <- object$data[id,]; valid <- object$data[-sort(unique(id)),]
     
        sl.b<-sl.time( # YOHANN2: Il faudra car les hyper-parametres soient fixes a chaque iteration
          methods=object$methods, metric=object$metric, data=learn, times="times",
          failures="failures", pro.time = pro.time, group=object$predictors$group,
          cov.quanti=object$predictors$cov.quanti,  cov.quali=object$predictors$cov.quali,
          cv=object$cv,  param.tune = object$param.tune$optimal,
          param.weights.fix=object$weights$coefficients, keep.predictions=FALSE, # Il faudra passer a FALSE
          verbose=FALSE)
       
          valid1 <- valid0 <- valid
          valid1[,object$predictors$group] <- 1; valid0[,object$predictors$group] <- 0
        
        # pas sur de utilite de BCV par rapport a Bootstrap simple ???
        pred1.b <- predict(sl.b, newdata=valid1, newtimes=object$times)
        pred0.b <- predict(sl.b, newdata=valid0, newtimes=object$times)
        
        # faut-il vraiment aller prendre les temps de echantillon apprentissage ???
        surv.sl1.b <- pred1.b$predictions[[length(object$methods)+1]][,pred1.b$times %in%
           sort(unique(learn[learn[,"failures"]==1, "times"]))]
        
        surv.sl0.b <- pred0.b$predictions[[length(object$methods)+1]][,pred0.b$times %in%
           sort(unique(learn[learn[,"failures"]==1, "times"]))]
        
        ti.b  <- pred1.b$times[pred1.b$times %in%
            sort(unique(learn[learn[,"failures"]==1, "times"]))]
        
        h.sl1.b <- h.sl0.b <- surv.sl1.b
        
    	logHRcond.i.b <- rep(-99, dim(surv.sl0.b)[1])

        for (j in 1:(dim(surv.sl0.b)[1])) {
        h.sl1.b[j,] <- differentiation(x=ti.b, fx=-1*log(surv.sl1.b[j,]))
        h.sl0.b[j,] <- differentiation(x=ti.b, fx=-1*log(surv.sl0.b[j,]))
        logHRcond.i.b[j] <- log(mean(h.sl1.b[j,]/h.sl0.b[j,]))  }
        
        logHRcond.b <- mean(logHRcond.i.b)
        
        return( list(logHRcond.b = logHRcond.b) ) 
        }
        
        if(!is.null(conda.env)) {
          use_condaenv(condaenv = conda.env, conda = "auto", required = FALSE) }
        
        cl <- makeCluster(n.cluster, type=cluster.type)
        registerDoParallel(cl)
        clusterEvalQ(cl, {library(RISCA)})
        res <- NULL
        res <- foreach(i = 1:iterations, .inorder=TRUE) %dopar% {set.seed(i); simul()}
        registerDoSEQ()
        stopCluster(cl)
   
    x<-paste0(rep("res[[", iterations), 1:iterations, rep("]]$logHRcond.b", iterations),
    		 collapse=", ")
    x<-paste0("c(", x, ")", collapse = "")
    logHRcond.b <- eval(parse(text=x))
  }

  mean.logHRcond <- mean(logHRcond.b, na.rm=TRUE)
 
  se.logHRcond <- sd(logHRcond.b, na.rm=TRUE)

  pv <- function(m, s){ ztest <- m/s;  return(ifelse(ztest<0,2*pnorm(ztest),2*(1-pnorm(ztest)))) }
  
  p.value.logHRcond <- pv(m=mean.logHRcond, s=se.logHRcond); if(p.value.logHRcond==0)
  				{ p.value.logHRcond <- "<0.001" }
  
  ci.low.logHRcond <- mean.logHRcond - qnorm(0.975, 0, 1)*se.logHRcond
  ci.upp.logHRcond <- mean.logHRcond + qnorm(0.975, 0, 1)*se.logHRcond
  
  .obj <- list(
    logHR.conditional.values = logHRcond.i,
    logHR.conditional = data.frame(estimate=mean.logHRcond, std.error = se.logHRcond,
      ci.lower=ci.low.logHRcond, ci.upper=ci.upp.logHRcond, p.value=p.value.logHRcond) )
      
  return(.obj)
}
