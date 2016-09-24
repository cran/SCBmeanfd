\name{scb.equal}
\alias{scb.equal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare Two Mean Functions}
\description{
This two-sample test builds simultaneous confidence bands (SCB) for the difference between two population mean functions and retains the equality assumption if the null function is contained in the bands. Equivalently, SCB are built around one of the local linear estimates (the one for say, population 1), and the equality hypothesis is  accepted if the other estimate (the one for  population 2) lies within the bands.
}
\usage{
scb.equal(x, y, bandwidth, level = .05, degree = 1, 
	scbtype = c("normal","bootstrap","both","no"), gridsize = NULL, 
	keep.y = TRUE, nrep = 2e4, nboot = 1e4, parallel = c("no","multicore","snow"), 
	ncpus = getOption("boot.ncpus",1L), cl = NULL)
}

\arguments{
  \item{x}{a list of length 2 or matrix with 2 columns containing the x values of each sample. If the two samples are observed on the same grid, \code{x} can be specified as a single vector. The range of \code{x} must be the same for each sample; missing values are not accepted.}

\item{y}{a list of length 2 containing matrices or data frames with functional observations (= curves) stored in rows. The number of columns of each component of \code{y} must match the length of the corresponding \code{x}. Missing values are not accepted.}
  
\item{bandwidth}{the kernel bandwidths (numeric vector of length 1 or 2).}

\item{level}{the significance level of the test (default = .05).}

\item{degree}{the degree of the local polynomial fit.}

\item{scbtype}{the type of simultaneous confidence bands to build: "normal", "bootstrap", "both", or "no". }

\item{gridsize}{the size of the grid over which the mean function is to be estimated. Defaults to the length of the smallest \code{x} grid.}

\item{keep.y}{logical; if \code{TRUE}, keep \code{y} in the result.}

\item{nrep}{the number of replicates for the normal SCB method (default = 20,000).}

\item{nboot}{the number of replicates for the bootstrap SCB method (default = 5,000).}

\item{parallel}{the computation method for the bootstrap SCB. By default, computations are sequential (\code{"no"}).\cr The function \code{\link[boot]{boot}} is used and can be run in parallel using the package \code{\link[parallel]{parallel}}.\cr Both options \code{"multicore"} and \code{"snow"} are available for parallel computing.} 

\item{ncpus}{the number of cores to use for parallel computing if \code{parallel = "multicore"}.}

\item{cl}{the name of the cluster to use for parallel computing if \code{parallel = "snow"}.}
}


\value{
A list object of \link[base]{class} \code{"SCBand"}. Depending on the  function used to create the object (\code{\link{scb.mean}}, \code{\link{scb.model}}, or \code{\link{scb.equal}}), some of its components are set to \code{NULL}. For \code{\link{scb.mean}}, the object has components: 

\item{\code{x}}{the argument \code{x}.
}
\item{y}{
if \code{keep.y} is \code{TRUE}, the argument \code{y}, else \code{NULL}.}
\item{\code{call}}{the function call.}
\item{\code{model}}{\code{NULL}.}  
\item{\code{par}}{\code{NULL}.}
\item{\code{nonpar}}{a list of two local linear estimates, one for each population.}
\item{\code{bandwidth}}{the argument \code{bandwidth}.}
\item{degree}{the degree of local polynomial used. 
Currently, only local linear estimation is supported.}
\item{level}{the argument \code{level}.}
\item{scbtype}{the argument \code{scbtype}.}
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
Degras, D. (2011). Simultaneous confidence bands for nonparametric regression with functional data. \emph{Statistica Sinica}, \bold{21}, 1735--1765.
}


\seealso{
\code{\link{scb.mean}}, \code{\link{scb.model}} 
}

\examples{
\dontrun{
# Phoneme data: compare the mean log-periodograms 
# for phonemes "aa" as the vowel in "dark" and "ao" 
# as the first vowel in "water"
data(phoneme)
n <- nrow(phoneme)
N <- ncol(phoneme) 
classes <- split(1:n,phoneme[,N])
names(classes) <- c("sh", "iy", "dcl", "aa", "ao")
freq    <- 1:150
compare.aa.ao <- scb.equal(freq, list(phoneme[classes$aa,-N], 
  phoneme[classes$ao,-N]), bandwidth = c(.75, .75), scbtype = "both", nboot = 2e3)
summary(compare.aa.ao)
}
}


