rmst <- function(times, surv.rates, max.time, type) {
  
  if(!(type %in% c("s", "l"))) {stop("Argument 'type' must contained 's' or 'l'")}
  
  if(type=="s") {
    
  if(times[1]!=0) {
        .t <- c(0, times[times <= max.time], min(max.time, max(times)))
        .s <- c(1, surv.rates[times <= max.time], surv.rates[length(surv.rates)]) }
    
  if(times[1]==0) {
    .t <- c(times[times <= max.time], min(max.time, max(times)))
    .s <- c(surv.rates[times <= max.time], surv.rates[length(surv.rates)]) }
  
   return(sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])) }
  
  if(type=="l") {
    
    if(max(times)<max.time) {stop("When type='l', max(time) must be higher than or equal to max.time")}
    
    if(times[1]!=0) {
      .t <- c(0, times)
      .s <- c(1, surv.rates) }
    
    if(times[1]==0) {
      .t <- times
      .s <- surv.rates }
    
    .fun <- spliner(.s ~ .t, monotonic = FALSE)
    
    return(cubintegrate(f = .fun, lower = 0, upper = max.time, method = "pcubature")$integral) }
  
}

# library(cubature)
# library(mosaic)
# toto <- function(x) {x^2}
# cubintegrate(f = toto, lower = 0, upper = 8, method = "pcubature")$integral
# x <- seq(-100:100)
# y <- x^2
# toto <- spliner(y ~ x, monotonic = FALSE)
# cubintegrate(f = toto, lower = 0, upper = 8, method = "pcubature")$integral

# x <- 0:100
# fonc <- function(x) {exp(-0.1*x)}
# y <- fonc(x)
# plot(x, y)
# integrate(fonc, 0, 20)
#rmst(times = x, surv.rates = y, max.time = 20, type="l")
#rmst(times = x, surv.rates = y, max.time = 20, type="s")



