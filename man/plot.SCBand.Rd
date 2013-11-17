\name{plot.SCBand}
\alias{plot.SCBand}

\title{Plot a SCBand Object}

\description{\code{plot} method for class \code{"SCBand"}.}

\usage{
\method{plot}{SCBand}(x , y = NULL, xlim = NULL, ylim = NULL, main = NULL, xlab = NULL, 
     ylab = NULL, col = NULL, cex = NULL, pch = NULL, lty = NULL, lwd = NULL, 
     legend = TRUE, where = NULL, text = NULL, legend.cex = NULL, horiz = TRUE, 
     bty = "n", ...)
}

\arguments{
  \item{x}{
a \code{SCBand} object, typically result of a call to \code{scb.mean}, \code{scb.model}, or \code{scb.equal}.}
  \item{y}{
optional y data.}
  \item{xlim}{
the x limits (x1, x2) of the plot.}
  \item{ylim}{
the y limits of the plot. 
}
  \item{main}{
an overall title for the plot.}
  \item{xlab}{
a title for the x axis.}
  \item{ylab}{
a title for the y axis.}
  \item{col}{
the colors for lines and points.}
  \item{cex}{
a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default.}
  \item{pch}{
a vector of plotting characters or symbols: see \code{\link[graphics]{points}}.}
  \item{lty}{
a vector of line types, see \code{\link[graphics]{par}}.}
  \item{lwd}{
a vector of line widths, see \code{\link[graphics]{par}}.
}
  \item{legend}{
logical; if \code{TRUE}, a legend is added to the plot.}
  \item{where}{
a character vector indicating the legend location: \code{"bottomright", "bottom"}, 
\code{"bottomleft", "left", "topleft", "top", "topright", "right"} 
or \code{"center"}.}
  \item{text}{a character vector to appear in the legend.
}
  \item{legend.cex}{
character expansion factor relative to current par("cex") for the legend. }
  \item{horiz}{
logical; if TRUE, the legend is displayed horizontally rather than vertically.
}
  \item{bty}{
the type of box to be drawn around the legend. The allowed values are \code{"n"} (the default) and \code{"o"}.}
  \item{\dots}{additional arguments passed to the function \code{\link[graphics]{matplot}} for displaying the y data.}
}

\details{
The argument \code{y} can be used to plot subsets of the y data. 
If non null, this argument has priority over the component \code{x$y} for plotting. 

The graphical parameters \code{col}, \code{cex}, \code{pch}, \code{lty}, and \code{lwd} apply to the following components to be plotted: data, parametric estimate, nonparametric estimate(s), normal simultaneous confidence bands (SCB), and bootstrap SCB. More precisely, \code{cex} and \code{pch} must be specified as vectors of length equal to the number of \eqn{y} data sets to be plotted (0, 1, or 2); \code{lty} and \code{lwd} must specified as numeric vectors of length equal to the total number of estimates and SCB components;  \code{col} applies to all components and should be specified accordingly. If necessary, graphical parameters are recycled to match the required length.

By default a legend is plotted horizontally at the bottom of the graph. 
}


\examples{
\dontrun{
## Plasma citrate data 
time <- 8:21
data(plasma)
h <- cv.select(time, plasma, degree = 1, interval = c(.5, 1))
scbplasma <- scb.mean(time, plasma, bandwidth = h, scbtype = "both",
                      gridsize = 100)
plot(scbplasma, cex = .2, legend.cex = .85, xlab = "Time of day", 
     ylab = "Concentration", main = "Plasma citrate data")
}
}

\keyword{methods}
