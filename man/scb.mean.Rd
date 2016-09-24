\name{scb.mean}
\alias{scb.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build Simultaneous Confidence Bands for Mean Functions
}

\description{
Fit a local linear estimator and build simultaneous confidence bands (SCB) for the mean of functional data.  
}
\usage{
scb.mean(x, y, bandwidth, level = .95, degree = 1, 
	scbtype = c("normal","bootstrap","both","no"), gridsize = length(x), 
	keep.y = TRUE, nrep = 2e4, nboot = 5e3, parallel = c("no", "multicore", "snow"), 
	ncpus = getOption("boot.ncpus",1L), cl = NULL)
}


\arguments{
\item{x}{a numeric vector of x data. Missing values are not accepted.}

\item{y}{a matrix or data frame with functional observations (= curves) stored in rows. The number of columns of \code{y} must match the length of \code{x}. Missing values are not accepted.}

\item{bandwidth}{the kernel bandwidth smoothing parameter.}

\item{level}{the level of the simultaneous confidence bands.}

\item{degree}{the degree of the local polynomial fit.}

\item{scbtype}{the type of simultaneous confidence bands to build: "normal", "bootstrap", "both", or "no".}

\item{gridsize}{the size of the grid used to evaluate the mean function estimates and SCB. Defaults to \emph{length(x)}.}

\item{keep.y}{logical; if \code{TRUE}, keep \code{y} in the result.}

\item{nrep}{number of replicates for the Gaussian SCB method (20,000 by default).}

\item{nboot}{number of replicates for the bootstrap SCB method (5,000 by default).}

\item{parallel}{the computation method for the SCB. By default, computations are sequential (\code{"no"}). The bootstrap method uses function \code{\link[boot]{boot}} and can be run in parallel using the package \code{\link[parallel]{parallel}}. In this case both options \code{"multicore"} and \code{"snow"} are available.}
 
\item{ncpus}{number of cores to use for parallel computing when \code{parallel = "multicore"}.}
\item{cl}{name of the cluster to use for parallel computing when \code{parallel = "snow"}.}
}


\details{
The local polynomial fitting uses a standard normal kernel and is implemented via the \code{\link[KernSmooth]{locpoly}} function.
Bootstrap SCB are implemented with the \code{\link[boot]{boot}} function and typically require more computation time than normal SCB.
}

\value{
An object of \link[base]{class} \code{"SCBand"}. To accommodate the different functions creating objects of this class (\code{\link{scb.mean}}, \code{\link{scb.model}}, and \code{\link{scb.equal}}), some components of the object are set to \code{NULL}. The component list is:

\item{\code{x}}{the x data.}
\item{\code{y}}{the y data if \code{keep.y} is \code{TRUE}, else \code{NULL}.}
\item{\code{call}}{the function call.}
\item{\code{model}}{\code{NULL}.}  
\item{\code{par}}{\code{NULL}.}
\item{\code{nonpar}}{a nonparametric estimate.}
\item{bandwidth}{the argument \code{bandwidth}.}
\item{degree}{the degree of local polynomial used. 
Currently, only local linear estimation is supported.}
\item{level}{the argument \code{level}.}
\item{scbtype}{the argument \code{type}.}
\item{teststat}{\code{NULL}.}
\item{pnorm}{\code{NULL}.} 
\item{pboot}{\code{NULL}.}  
\item{qnorm}{the quantile used to build the normal SCB.}
\item{qboot}{the quantile used to build the bootstrap SCB.}
\item{normscb}{a matrix containing the normal SCB stored in columns.}
\item{bootscb}{a matrix containing the bootstrap SCB stored in columns.}
\item{gridsize}{the argument \code{gridsize} if nonnull, else  \emph{length(x)}.}
\item{nrep}{the argument \code{nrep}.}
\item{nboot}{the argument \code{nboot}.}


Depending on the value of \code{scbtype}, some of  
the fields \code{qnorm}, \code{normscb}, \code{nrep}, \code{qboot}, \code{normboot} and \code{nboot} may be \code{NULL}.
}

\references{
Degras, D. (2011). Simultaneous confidence bands for nonparametric regression with functional data. \emph{Statistica Sinica}, \bold{21}, 1735--1765.
}




\seealso{
\code{\link{scb.equal}}, \code{\link{scb.model}} 
}


\examples{
\dontrun{
## Plasma citrate data
data(plasma)
time <- 8:21
h <- cv.select(time, plasma, 1, c(.5, 1))
scbplasma <- scb.mean(time, plasma, bandwidth = h, scbtype = "both", gridsize = 100)
scbplasma
plot(scbplasma, cex = .2, legend.cex = .85, xlab = "Time", ylab = "Concentration", 
  main = "Plasma citrate data")
}
}
