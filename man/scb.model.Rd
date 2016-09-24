\name{scb.model}
\alias{scb.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Goodness-Of-Fit of a Model for the Mean Function}
\description{
This is the goodness-of-fit test for parametric models of the mean function described in Degras (2011). The candidate model must be a finite-dimensional function space (curvilinear regression). The test is based on the sup-norm distance between a smoothed parametric estimate and a local linear estimate. Graphically, the candidate model is retained whenever one of the estimates lies within the SCB built around the other.  
}
\usage{
scb.model(x, y, model, bandwidth, level = .05, degree = 1, 
	scbtype = c("normal","bootstrap","both","no"), gridsize = length(x), 
	keep.y = TRUE, nrep = 2e4, nboot = 5e3, parallel = c("no", "multicore", "snow"), 
	ncpus = getOption("boot.ncpus",1L), cl = NULL)
}

\arguments{
\item{x}{a numeric vector of x data. \code{x} must be a uniform grid; missing values are not accepted.}

\item{y}{a matrix or data frame with functional observations (= curves) stored in rows. The number of columns of \code{y} must match the length of \code{x}. Missing values are not accepted.}

\item{model}{an integer specifying the degree of a polynomial basis, or a data frame/matrix containing the basis functions stored in columns. In the latter case, the basis functions must be evaluated on a uniform grid of size \code{gridsize} spanning the range of \code{x}.}
  
\item{bandwidth}{the kernel bandwidth smoothing parameter.}

\item{level}{the significance level of the test (default = .05).}

\item{degree}{the degree of the local polynomial fit.}

\item{scbtype}{the type of simultaneous confidence bands to build: "normal", "bootstrap", "both", or "no".}

\item{gridsize}{the size of the grid over which the mean function is to be estimated. Defaults to \code{length(x)}.}

\item{keep.y}{logical; if \code{TRUE}, keep \code{y} in the result.}

\item{nrep}{the number of replicates for the normal SCB method (default = 20,000).}

\item{nboot}{the number of replicates for the bootstrap SCB method (default = 5,000).}

\item{parallel}{the computation method for the bootstrap SCB. By default, computations are sequential (\code{"no"}).\cr The function \code{\link[boot]{boot}} is used and can be run in parallel using the package \code{\link[parallel]{parallel}}.\cr Both options \code{"multicore"} and \code{"snow"} are available for parallel computing.} 

\item{ncpus}{the number of cores to use for parallel computing when \code{parallel = "multicore"}.}

\item{cl}{the name of the cluster to use for parallel computing when \code{parallel = "snow"}.}}

\value{
An object of \link[base]{class} \code{"SCBand"}. To accommodate the different functions creating objects of this class (\code{\link{scb.mean}}, \code{\link{scb.model}}, and \code{\link{scb.equal}}), some components of the object are set to \code{NULL}. The component list is:

\item{\code{x}}{the argument \code{x}.}
\item{\code{y}}{the argument \code{y} if \code{keep.y} is \code{TRUE}, else \code{NULL}.}
\item{\code{call}}{the function call.}
\item{\code{model}}{the argument \code{model}.}  
\item{\code{par}}{a smoothed parametric estimate.}
\item{\code{nonpar}}{a local linear estimate.}
\item{\code{bandwidth}}{the argument \code{bandwidth}.}
\item{degree}{the degree of the local polynomial. 
Currently, only local linear estimation is supported.}
\item{level}{the argument \code{level}.}
\item{scbtype}{the argument \code{type}.}
\item{teststat}{the test statistic.}
\item{pnorm}{the \emph{p} value for the normal-based statistical test.} 
\item{pboot}{the \emph{p} value for the boostrap-based statistical test.}  
\item{qnorm}{the quantile used to build the normal SCB.}
\item{qboot}{the quantile used to build the bootstrap SCB.}
\item{normscb}{a matrix containing the normal SCB stored in columns.}
\item{bootscb}{a matrix containing the bootstrap SCB stored in columns.}
\item{gridsize}{the argument \code{gridsize}, or \emph{length(x)} if no argument was specified.}
\item{nrep}{the argument \code{nrep}.}
\item{nboot}{the argument \code{nboot}.}


Depending on the value of \code{scbtype}, some or all of  
the fields \code{pnorm}, \code{qnorm}, \code{normscb}, \code{nrep},  \code{pboot}, \code{qboot}, \code{normboot} and \code{nboot} may be   \code{NULL}.

}

\references{
Degras, D. (2011). Simultaneous confidence bands for nonparametric regression with functional data.\cr \emph{Statistica Sinica}, \bold{21}, 1735--1765.
}

\seealso{
\code{\link{scb.equal}}, \code{\link{scb.mean}} 
}




\examples{
## Example from Degras (2011)
## Gaussian process with polynomial mean function 
## and Ornstein-Uhlenbeck covariance function
## The SCB and PLRT tests are compared

set.seed(100)
p    <- 100  	# number of observation points
x 	  <- seq(0, 1, len = p)
mu	  <- 10 * x^3 - 15 * x^4 + 6 * x^5	# mean 
R 	  <- (.25)^2 * exp(20 * log(.9) * abs(outer(x,x,"-"))) # covariance 
eigR <- eigen(R, symmetric = TRUE)  	
simR <- eigR$vectors \%*\% diag(sqrt(eigR$values)) 	 

# Candidate model for mu: polynomial of degree <= 3
# This model, although incorrect, closely approximates mu.
# With n = 50 curves, the SCB and PLRT incorrectly retain the model.
# With n = 70 curves, both tests reject it. 
n <- 50  
y <- mu + simR \%*\% matrix(rnorm(n*p), p, n) 	# simulate data  
y <- t(y) 	# arrange the trajectories in rows
h <- cv.select(x, y, 1)	
scb.model(x, y, 3, bandwidth = h)		  # p value: .652
plrt.model(x, y, 3, verbose = TRUE)	# p value: .450 
n <- 70  	
y <- mu + simR \%*\% matrix(rnorm(n*p), p, n) 		
y <- t(y) 										
h <- cv.select(x, y, 1)	
scb.model(x, y, 3, bandwidth = h)		  # p value: .004
plrt.model(x, y, 3, verbose = TRUE)	# p value: .001

# Correct model: polynomials of degree <= 5
scb.model(x, y, 5, bandwidth = h)  	# p value: .696
plrt.model(x, y, 5, verbose = TRUE)	# p value: .628

}

