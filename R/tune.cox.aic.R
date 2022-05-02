# CS splines ok
tune.cox.aic<- function(times, failures, group=NULL, cov.quanti=NULL, cov.quali=NULL, data, mini.model.cov=NULL, 
                        maxi.model.cov=NULL){
  # CS YF prevoir tous les checks possibles pour mes modeles en mini et maxi
  
  .outcome <- paste("Surv(", times, ",", failures, ")")
  if(is.null(mini.model.cov)==TRUE){
    if(!(is.null(group))){
      .f0<-as.formula( paste(.outcome, "~",group ))
    }
    else{
      .f0<-as.formula( paste(.outcome, "~1" ))
    }
  }
  else{
    .f0<-as.formula( paste(.outcome, "~",  paste(mini.model.cov, collapse = " + "),collapse = " ") )
  }
  if(is.null(maxi.model.cov)==TRUE){
    if(!(is.null(group))){
      if(is.null(cov.quanti)==F & is.null(cov.quali)==F){
        .f <- as.formula( paste(.outcome, "~", group, "*(",paste(cov.quanti, collapse = " + "), " + ",
                                paste(cov.quali, collapse = " + "),  ")",collapse = " ") )
      }
      if(is.null(cov.quanti)==F & is.null(cov.quali)==T){
        .f <- as.formula( paste(.outcome, "~", group, "*(",paste(cov.quanti, collapse = " + "),")" ))
      }
      if(is.null(cov.quanti)==T & is.null(cov.quali)==F){
        .f <- as.formula( paste(.outcome, "~", group, "*(",paste(cov.quali, collapse = " + "),  ")",collapse = " ") )
      }
      if(is.null(cov.quanti)==T & is.null(cov.quali)==T){
        .f <- as.formula( paste(.outcome, "~", group) )
      }
    }
    else{
      if(is.null(cov.quanti)==F & is.null(cov.quali)==F){
        .f <- as.formula( paste(.outcome, "~", paste(cov.quanti, collapse = " + "), " + ", paste(cov.quali, collapse = " + "),
                                collapse = " ") )
      }
      if(is.null(cov.quanti)==F & is.null(cov.quali)==T){
        .f <- as.formula( paste(.outcome, "~", paste( cov.quanti, collapse = " + ")))
      }
      if(is.null(cov.quanti)==T & is.null(cov.quali)==F){
        .f <- as.formula( paste(.outcome, "~",  paste(cov.quali, collapse = " + "),collapse = " ") )
      }
    }
    
  }
  else{
    .f<-as.formula( paste(.outcome, "~",  paste(maxi.model.cov, collapse = " + "),collapse = " ") )
  }





  
  .fit0<-coxph(.f0, data=data)
  
  .fit<-coxph(.f, data=data)
  
  .res.step<-stepAIC(.fit0, scope=formula(.fit), direction="forward", k=2, trace=FALSE)
  

  return(list(optimal=list(final.model.cov=names(.res.step$assign)),
              results = list(res.step=.res.step))
  )
  
  
}



