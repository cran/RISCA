rmst <- function(times, surv.rates, max.time) {
        .t <- c(0, times[times <= max.time], min(max.time, max(times)))
        .s <- c(1, surv.rates[times <= max.time], surv.rates[length(surv.rates)])
        return(sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s[1:(length(.s) - 1)])) }
