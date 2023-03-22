
print.sl.time <- function (x, ..., digits=7) 
{
  cat("\n", "The contribution of the leaners: ", "\n\n", sep = "")
  
  print(data.frame(leaners = x$methods, weights = round(x$weights$values, digits)))
  
  cat("\n", "The estimations were based on a ", x$cv, "-fold CV of the loss function ",
      x$metric$metric, ". Its value was ", round(x$metric$value, digits), ".", sep = "")
}
