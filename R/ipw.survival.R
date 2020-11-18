
ipw.survival <- function(times, failures, variable, weights=NULL)
    {
	if(sum(times<0)>0){print("Error times must be positive")
	}else{
		if(sum(weights<=0)>0){print("Error weights must be superior to 0")
		}else{
			if(sum(failures!=0 & failures!=1)>0){print("Error failures must be a vector of 0 or 1")
			}else{
				if(is.null(weights)){ .w <- rep(1, length(times)) }	else { .w <- weights }
				.data <- data.frame(t=times, f=failures, v=variable, w=.w)
				.data <- .data[!is.na(.data$v),]
				table.surv <- data.frame(times=NULL, n.risk=NULL, n.event=NULL, survival=NULL, variable=NULL)
					for(i in  unique(variable)){
					.d <- .data[.data$v==i,]
					.tj <- c(0,sort(unique(.d$t[.d$f==1])),max(.d$t))
					.dj <- sapply(.tj, function(x){sum(.d$w[.d$t==x & .d$f==1])})
					.nj <- sapply(.tj, function(x){sum(.d$w[.d$t>=x])})
					.st <- cumprod((.nj-.dj)/.nj)
					table.surv <- rbind(table.surv, data.frame(times = .tj, 
                    n.risk = .nj, n.event = .dj, survival = .st, 
                    variable = i))
					}
				.obj <- list(table.surv=table.surv)
				class(.obj) <- "survrisca"
				return(.obj)
			}
		}
	}
    }
