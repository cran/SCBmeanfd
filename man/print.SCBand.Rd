\name{print.SCBand}
\alias{print.SCBand}

\title{
Print a SCBand Object}
\description{
\code{print} method for class \code{"SCBand"}.
}
\usage{\method{print}{SCBand}(x,...)}

\arguments{
  \item{x}{an object of class \code{"SCBand"}.}
  \item{...}{for compatibility with the generic \code{print} method; argument not currently used.}
}

\details{
The function \code{print.SCBand} concisely displays the information of an object of class \code{"SCBand"}. More precisely it shows the  
data range, bandwidth used in local polynomial estimation, and key information on SCB and statistical tests. 
}


\seealso{
\code{\link{plot.SCBand}}, \code{\link{summary.SCBand}}
}

\examples{
\dontrun{
# Plasma citrate data
data(plasma)
time <- 8:21
h <- cv.select(time, plasma, 1)
scbplasma <- scb.mean(time, plasma, bandwidth = h, scbtype = "both", gridsize = 100)
scbplasma
}
}

\keyword{print}
\keyword{methods}
