\name{plot.survrisca}
\alias{plot.survrisca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Method for 'survrisca' Objects
}
\description{
A plot of survival curves is produced. In the RISCA package, it concerns the functions \code{ipw.survival} and \code{gc.survival}.
}
\usage{
\method{plot}{survrisca}(x, ..., col=1, lty=1, lwd=1, type="s", max.time=NULL,
min.y=0, max.y=1, grid.lty=NULL)
}
\arguments{
  \item{x}{An object of class \code{survrisca}, returned by the functions \code{ipw.survival} \code{gc.survival}.}
  \item{...}{Additional arguments affecting the plot.}
  \item{col}{A numeric vector with the color of the survival curves. The default is 1 for black.}
  \item{lty}{A numeric vector with the type of the survival curves. The default is 1.}
  \item{lwd}{A numeric vector with the type of the survival curves. The default is 1.}
  \item{type}{A character string giving the type of plot : "p" for points, "l" for lines, "b" for both points and lines, "c" for empty points joined by lines, "o" for overplotted points and lines. The default is "s" for step function.}
  \item{max.time}{The maximum time of the x-asis. The default is NULL, it corresponds to the maximum follow-up time observed in the database from which the \code{survrisca} object was estimated.}
  \item{min.y}{The minimum of the y-axis. The default is 0.}
  \item{max.y}{The maximum of the y-axis. The default is 1.}
  \item{grid.lty}{A character or (integer) numeric with the line type of the grid lines. The default is NULL for no grid.}
}
\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>

Florent Le Borgne <fleborgne@idbc.fr>
}

\examples{

data(dataDIVAT2)

res.km <- ipw.survival(times=dataDIVAT2$times, failures=dataDIVAT2$failures,
  variable=dataDIVAT2$ecd, weights=NULL)

plot(res.km, ylab="Graft and patient survival",
 xlab="Time post-transplantation (years)", col=c(1,2), grid.lty=1)
}

\keyword{Plotting}