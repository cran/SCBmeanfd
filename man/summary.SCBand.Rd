\name{summary.SCBand}
\alias{summary.SCBand}

\title{Summarize a SCBand Object}

\description{\code{summary} method for class \code{"SCBand"} }

\usage{\method{summary}{SCBand}(object, ...)}

\arguments{\item{object}{an object of class \code{"SCBand"}}
\item{...}{additional arguments; not currently used.}
}

\details{
The function \code{summary.SCBand} displays all fields of a \code{SCBand} object at the exception of \code{x}, \code{y}, \code{par}, \code{nonpar}, \code{normscb}, and \code{bootscb} which are potentially big. It provides information on the function call, data, local polynomial fit, SCB, and statistical tests.
}

\seealso{
\code{\link{plot.SCBand}},
\code{\link{print.SCBand}} 
}

\examples{
\dontrun{
## Plasma citrate data 
data(plasma)
time <- 8:21
h <- cv.select(time, plasma, 1)
scbplasma <- scb.mean(time, plasma, bandwidth = h, scbtype = "both", gridsize = 100)
summary(scbplasma)
}
}


\keyword{methods}
