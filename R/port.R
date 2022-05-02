#library(rpart)
#library(rpart.plot)


port<-function(group, cov.quanti, cov.quali, data, alpha=0.05, beta=0.05, gamma=2, pruning=FALSE,
                minbucket=6, minsplit=20, maxdepth=30) {
  
  if(!is.data.frame(data)){ 
    stop("The argument \"data\" need to be a data.frame or a matrix")
  }
  
  if(alpha>0.5 | alpha<0) stop("The argument \"alpha\" must be a proportion (e.g., 0.05 for 5%).")
  if(beta>1 | beta<=0) stop("The argument \"beta\" must be a non-null proportion (e.g., 0.05 for 5%).")
  
  ## Change to numeric if exposure is coded as a factor
  if(!is.numeric(data[,group]) & !is.integer(data[,group])){
    stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the argument \"group\".")
  }
  
  data[,cov.quali] <- apply(data[,cov.quali], 2, as.factor)
  
  ########### Internal functions
  
  
  
  ### Remove superfluous cutoffs : subgrp is the problematic subgroup, type is continuous or 
  ### categorical
  define_cutoff<-function(subgrp, var, treatment, type, data, pruning){
    sentence<-character(0) 
    ##Extract quantities used in output sentences further down
    n<-subgrp$n[nrow(subgrp)]    
    pourcent<-round(n/subgrp$n[1],digits=3)*100
    proba<-round(subgrp$proba[nrow(subgrp)],digits=3)*100  ## probability in percentage
    ## For consistency reasons we always consider the low-probability case
    treat<-ifelse(proba>=50, "unexposed", "exposed") #respect. 0 and 1
    proba<-ifelse(proba>=50,100-proba,proba)
    
    for(v in var){
      if(which(var == v) == 1){
        
        sentence<-paste0(sentence,"\n - individuals ")
        
      }else if(which(var == v) == length(var)){
        
        sentence<-paste0(sentence," and ")
        
      }else{
        sentence<-paste0(sentence,", ")
      }
      ## Continuous variable 
      if(type[var == v] == "continuous"){
        ## On garde les lignes avec var, var2 est de la forme var > valeur ou var < valeur
        imcs<-subgrp$var2[grepl(v,subgrp$var2,fixed=TRUE)]
        
        ##We need to memorise the positions of the elements in the subgrpset to later remove them
        loc_imcs<-which(subgrp$var2 %in% imcs)
        
        ## We remove the variable name
        imcs<-sapply(imcs,function(i) gsub(v,"",i))
        ## We divide into upper and lower bound
        imcs_sup<-imcs[grepl("<",imcs)];loc_sup<-loc_imcs[grepl("<",imcs)]
        imcs_inf<-imcs[grepl(">",imcs)];loc_inf<-loc_imcs[grepl(">",imcs)]
        
        newcut<-subgrp
        ## Violation removed if it is contained between two values
        if(pruning){
          if(length(imcs_sup) > 0 && length(imcs_inf) > 0){
            newcut<-NULL
            return(list(subgrp=newcut,sentence=character(0)))
          }
        }
        ## Remove character of size 2 (">=" or "> ") to only have cutoff value
        imcs_sup<-as.numeric(vapply(imcs_sup,function(i) substr(i,start=3, stop=nchar(i)),FUN.VALUE=character(1)))
        imcs_inf<-as.numeric(vapply(imcs_inf,function(i) substr(i,start=3, stop=nchar(i)),FUN.VALUE=character(1)))
        
        
        if(length(imcs_sup) > 1 || length(imcs_inf) > 1){
          ## Only keep maximum and minimum values
          ifelse(length(imcs_sup)>0,imc_sup<-min(imcs_sup),imc_sup<-imcs_sup) 
          ifelse(length(imcs_inf)>0,imc_inf<-max(imcs_inf),imc_inf<-imcs_inf)
          keep<-c(loc_sup[which(imcs_sup == imc_sup)],loc_inf[which(imcs_inf == imc_inf)])
          
          ## We readjust because rpart prints the halfway point for bounds
          if(is.integer(data[,v])){
            imc_sup <- ceiling(imc_sup) 
            imc_inf <- floor(imc_inf) 
          }
          
          
          ## Reduce our subgrp set
          newcut<-newcut[-setdiff(loc_imcs,keep),]
          
          ## Write output sentences
          if(length(imc_sup) == 0 & length(imc_inf) > 0){
            
            sentence<-paste0(sentence,"with values greater than or equal than ",imc_inf," on variable ",v)
            
          }
          if(length(imc_sup) > 0 & length(imc_inf) == 0){
            
            sentence<-paste0(sentence,"with values less than ",imc_sup," on variable ",v)
            
          }
          if((length(imc_sup) > 0 & length(imc_inf) > 0) && !pruning){
            
            sentence<-paste0(sentence,"with values greater than or equal than ",imc_inf," and less than ",imc_sup," on variable ",v)
            
          }
        }else{        ## If no shortening is needed
          
          ## We readjust by half the difference between categories because rpart sends the halfway point
          if(is.integer(data[,v])){
            imcs_sup <- ceiling(imcs_sup) 
            imcs_inf <- floor(imcs_inf) 
          }
          
          if(length(imcs_sup) == 0 & length(imcs_inf) == 1){
            
            sentence<-paste0(sentence,"with values greater than or equal to ",imcs_inf,
                             " on variable ",v)
            
          }
          if(length(imcs_sup) == 1 & length(imcs_inf) == 0){
            
            sentence<-paste0(sentence,"with values less than ",imcs_sup," on variable ",v)
            
          }
          if((length(imcs_sup) == 1 & length(imcs_inf) == 1) && !pruning){
            
            sentence<-paste0(sentence,"with values greater than or equal than ",imcs_inf,
                             " and less than ",imcs_sup," on variable ",v)
            
          }
        }
        ## Remove the root (first line) if necessary
        if("root" %in% newcut$var2 && length(newcut$var2) >1){
          newcut<-newcut[-1,]
        }
      }
      ## Categorical variable
      if(type[var==v] == "categorical"){
        ## Select the last knot belonging to this category
        index<-grep(v,subgrp$var2,fixed=TRUE)
        l<-index[-length(index)]
        newcut<-subgrp
        ## Remove all other lines 
        if(length(l) > 0){
          newcut<-newcut[-l,]
        }
        ## Remove the root (first line) if necessary
        if("root" %in% subgrp$var2){
          newcut<-newcut[-1,]
        }
        cat_values<-substr(subgrp$var2[index[length(index)]],start=nchar(v)+2, stop=nchar(subgrp$var2[index[length(index)]]))
        
        sentence<-paste0(sentence,"of categories ",cat_values," on variable ",v)
        
      }
    }    
    
    sentence<-paste0(sentence," (N=",n,", ",pourcent,"% of the sample), have a ",
                     sprintf("%.1f",proba),"% chance of being ", treat,".")
    
    return(list(subgrp=newcut,sentence=sentence))
  }
  
  ## Modified labelling function which put letters for categorical variables
  labels.rpart<-function (object, digits = 4, minlength = 1L, pretty, collapse = TRUE, 
                          ...) 
  {
    if (missing(minlength) && !missing(pretty)) {
      minlength <- if (is.null(pretty)) 
        1L
      else if (is.logical(pretty)) {
        if (pretty) 
          4L
        else 0L
      }
      else 0L
    }
    ff <- object$frame
    n <- nrow(ff)
    if (n == 1L) 
      return("root")
    is.leaf <- (ff$var == "<leaf>")
    whichrow <- !is.leaf
    vnames <- ff$var[whichrow]
    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + !is.leaf))
    irow <- index[c(whichrow, FALSE)]
    ncat <- object$splits[irow, 2L]
    lsplit <- rsplit <- character(length(irow))
    if (any(ncat < 2L)) {
      jrow <- irow[ncat < 2L]
      
      formatg <- function(x, digits = getOption("digits"),
                          format = paste0("%.", digits, "g"))
      {
        if (!is.numeric(x)) stop("'x' must be a numeric vector")
        
        temp <- sprintf(format, x)
        if (is.matrix(x)) matrix(temp, nrow = nrow(x)) else temp
      }
      
      cutpoint <- formatg(object$splits[jrow, 4L], digits)
      temp1 <- (ifelse(ncat < 0, "< ", ">="))[ncat < 
                                                2L]
      temp2 <- (ifelse(ncat < 0, ">=", "< "))[ncat < 
                                                2L]
      lsplit[ncat < 2L] <- paste0(temp1, cutpoint)
      rsplit[ncat < 2L] <- paste0(temp2, cutpoint)
    }
    if (any(ncat > 1L)) {
      xlevels <- attr(object, "xlevels")
      jrow <- seq_along(ncat)[ncat > 1L]
      crow <- object$splits[irow[ncat > 1L], 4L]
      cindex <- (match(vnames, names(xlevels)))[ncat > 1L]
      if (minlength == 1L) {
        if (any(ncat > 52L)) 
          warning("more than 52 levels in a predicting factor, truncated for printout", 
                  domain = NA)
      }
      else if (minlength > 1L) 
        xlevels <- lapply(xlevels, abbreviate, minlength, 
                          ...)
      for (i in seq_along(jrow)) {
        j <- jrow[i]
        splits <- object$csplit[crow[i], ]
        cl <- if (minlength == 1L)
          ""
        else ","
        lsplit[j] <- paste((xlevels[[cindex[i]]])[splits ==
                                                    1L], collapse = ",")
        rsplit[j] <- paste((xlevels[[cindex[i]]])[splits ==
                                                    3L], collapse = ",")
      }
    }
    if (!collapse) {
      ltemp <- rtemp <- rep("<leaf>", n)
      ltemp[whichrow] <- lsplit
      rtemp[whichrow] <- rsplit
      return(cbind(ltemp, rtemp))
    }
    lsplit <- paste0(ifelse(ncat < 2L, "", "="),
                     lsplit)
    rsplit <- paste0(ifelse(ncat < 2L, "", "="),
                     rsplit)
    varname <- (as.character(vnames))
    node <- as.numeric(row.names(ff))
    parent <- match(node%/%2L, node[whichrow])
    odd <- (as.logical(node%%2L))
    labels <- character(n)
    labels[odd] <- paste0(varname[parent[odd]], rsplit[parent[odd]])
    labels[!odd] <- paste0(varname[parent[!odd]], lsplit[parent[!odd]])
    labels[1L] <- "root"
    labels
  }
  
  #####################################
  
  
  covariates <- c(cov.quanti, cov.quali)
  
  ## output lists
  problem_covariates<-problem_cutoffs<-problem_trees<-problem_nodes <- list()
  
  n <- gamma
  covariables <- covariates
  m <- length(covariables)
  combi <- list()
  p <- 1
    for(i in 1:n){
      if(m >= i){
        comb <- combn(m,i)
        for(j in 1:ncol(comb)){
          combi[[as.character(p)]] <- comb[,j]
          p <- p+1
        }  } else{
        combi[[as.character(p)]] <- seq(1,m,by=1)   }
    }

  plus1<-FALSE
  for(q in 1:length(combi)){
    ## Remove variables involved in a one-variable violation
    if(q > length(covariates) & !plus1){
      plus1<-TRUE
      if(length(problem_cutoffs)>0){
        var_prob<-character(0)
        for(k in 1:length(problem_covariates)){
          var_prob<-c(var_prob,problem_covariates[[as.character(k)]])
        }
        bad_cov<-which(covariates %in% var_prob)
      }
    }
    ## Skip variables already involved in a violation
    if(exists("bad_cov") && sum(bad_cov %in% combi[[as.character(q)]]) > 0){
      next
    }
    
    ## Subset of covariates
    covariables<-covariates[combi[[as.character(q)]]]
    ## rpart formula
    regression<-as.formula(paste(group,"~",paste(covariables,collapse="+")))
    
    ## maximal CART tree, with chosen minbucket, minsplit and maxdepth parameters
    cart_max<-rpart(regression,data=data,cp=0,minbucket=minbucket,method="anova",
                    minsplit=minsplit,maxdepth=maxdepth)
    
    
    ## frame contains tree's essential information as a data.frame, one node per line,
    ## with number of individuals and probability
    frame<-cart_max$frame
    ## var2 gives cutoffs written as variable >= or > value
    frame$var2<-labels.rpart(cart_max)
    ## We keep a column with only the name of the variable to make it easier to refer to
    if(nrow(frame)>1){
      for(j in 2:nrow(frame)){
        for(cov in covariables){
          if(grepl(cov,frame$var2[j])){
            frame$var[j]<-cov
          }
        }
      }
    }
    problematic_nodes<-numeric(0)
    ## Look for node with sufficient number of individuals and extreme probability
    for(i in 1:nrow(frame)){
      if((frame$yval[i] >= 1-alpha || frame$yval[i] <= alpha)
         && frame$n[i] >= nrow(data) * beta){
        problematic_nodes<-c(problematic_nodes,as.numeric(rownames(frame)[i]))
      }
    }
    
    ## Gives path from node x to root node
    parent <- function(x) {
      if (x[1] != 1)
        c(Recall(if (x %% 2 == 0L) x / 2 else (x - 1) / 2), x) else x
    }
    problematic_path<-list()
    
    ## Only keep highest node in a path 
    for(n in problematic_nodes){
      for(p in parent(n)[-length(parent(n))]){    ## Remove the last element of the path
        if(p %in% problematic_nodes){
          problematic_nodes <- problematic_nodes[!problematic_nodes == n]
        }
      }
    }
    
    ## Save values from frame for problem nodes
    for(i in problematic_nodes){
      problematic_path[[as.character(i)]]<-frame[which(rownames(frame) %in% parent(i)),c("var","var2","n")]
      problematic_path[[as.character(i)]]$proba<-frame$yval[which(rownames(frame) %in% parent(i))]
    }
    
    problem_names<-1
    for(i in names(problematic_path)){
      problem_names<-c(rownames(problematic_path[[i]]),problem_names)
    }
    problem_names<-unique(problem_names)
    if(length(problematic_path) != 0){
      for(k in 1:length(problematic_path)){
        if(length(problem_cutoffs) == 0 || sum(sapply(problem_cutoffs,
                        problematic_path[k][[names(problematic_path[k])]], FUN=identical)) == 0){
          problem_cutoffs[[as.character(length(problem_cutoffs) + 1)]]<-problematic_path[k][[names(problematic_path[k])]]
          problem_covariates[[as.character(length(problem_covariates)+1)]]<-covariables
          problem_trees[[as.character(length(problem_trees) + 1)]]<-cart_max
          problem_nodes[[as.character(length(problem_nodes) + 1)]]<-problem_names
        }
      }
    }
  }
  
  if(length(problem_cutoffs) > 0){
    
    sentence<-"In order for potential confounding factors to be correctly handled, one must ensure that individual characteristics are sufficiently frequent in each group. We have observed that :"
    
    ### Shortening of problematic paths
    for(i in 1:length(problem_cutoffs)){
      vars<-unique(problem_cutoffs[[as.character(i)]]$var)
      type<-character(length(vars))
      for(k in 1:length(vars)){
        type[k]<-ifelse(vars[k] %in% cov.quanti, "continuous", "categorical")
      }
      problem_covariates[[as.character(i)]]<-paste(problem_covariates[[as.character(i)]],collapse=";")
      
      cut<-define_cutoff(problem_cutoffs[[as.character(i)]],vars,data[,group], type, data, pruning=pruning)
      ### Remove element from lists if it's a violation between two bounds
      if(is.null(cut$data)){
        problem_cutoffs[[as.character(i)]]<-NULL
        problem_covariates[[as.character(i)]]<-NULL
        problem_trees[[as.character(i)]]<-NULL
        problem_nodes[[as.character(i)]]<-NULL
      }else{
        problem_cutoffs[[as.character(i)]]<-cut$data
      }
      if(length(cut$sentence)>0){
        sentence<-paste(sentence,cut$sentence)
      }
    }
  }else{  ## If no issues are detected
    sentence<-"No potential positivity violation identified."
  }
  ## This condition is usually verified if all observed violations were between two bounds
  if(sentence=="In order for potential confounding factors to be correctly handled, one must ensure that individual characteristics are sufficiently frequent in each group. It is observed that :"){ sentence<-"No potential positivity violation identified." }
  

  return(cat(paste(sentence, "\n")))
  
}


