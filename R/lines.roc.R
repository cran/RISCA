
# lines.roc <- function(x)  UseMethod("lines.roc")

lines.roc <- function(x, ...){
lines(1-x$table$sp, x$table$se, ...)
}
