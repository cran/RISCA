
lrs.multistate <- function(model1, model0)
{
statistic <- 2*(model1$LogLik-model0$LogLik)

if (statistic<0) stop("The statistic value have to be positive")

ddl <- dim(model1$table)[1] - dim(model0$table)[1]

if (ddl<=0) stop("The number of degrees of freedom have to be strictly positive")

pvalue <- 1-pchisq(statistic, ddl)

return(c(statistic=statistic, ddl=ddl, pvalue=pvalue))
}

  