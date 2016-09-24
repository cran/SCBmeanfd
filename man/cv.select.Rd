\name{cv.select}
\alias{cv.select}

\title{
Cross-Validation Bandwidth Selection for Local Polynomial Estimation}

\description{
Select the cross-validation bandwidth described in Rice and Silverman (1991) for the local polynomial estimation of a mean function based on functional data. 
}

\usage{ cv.select(x, y, degree = 1, interval = NULL, gridsize = length(x), ...) }

\arguments{
\item{x}{observation points. Missing values are not accepted.}
\item{y}{matrix or data frame with functional observations (= curves) stored in rows. The number of columns of \code{y} must match the length of \code{x}. Missing values are not accepted.}
\item{degree}{degree of the local polynomial fit.}
\item{interval}{lower and upper bounds of the search interval (numeric vector of length 2).}
\item{gridsize}{size of evaluation grid for the smoothed data.}  
\item{...}{additional arguments to  pass to the optimization function \code{\link[stats]{optimize}}.}
}

\details{
The cross-validation score is obtained by leaving in turn each curve out and computing the prediction error of the local polynomial smoother based on all other curves. For a bandwith value \eqn{h}, this score is 
\deqn{ CV(h) = \frac{1}{np} \sum_{i=1}^n \sum_{j=1}^p \left( Y_{ij} - \hat{\mu}^{-(i)}(x_j;h) \right)^2, }{CV(h) =  \sum (Y[ij] - \mu^{-i}(x[j];h))^2 / (n*p)}
where \eqn{Y_{ij}}{Y[ij]} is the measurement of the \eqn{i}-th curve at location \eqn{x_j}{x[j]} for \eqn{i=1,\ldots,n} and \eqn{j=1,\ldots,p}, and \eqn{\hat{\mu}^{-(i)}(x_j;h)}{\mu^{-i}(x[j];h)} is the local polynomial estimator with bandwidth \eqn{h} based on all curves except the \eqn{i}-th.

If the \code{x} values are not equally spaced, the data are first smoothed and evaluated on a grid of length \code{gridsize} spanning the range of \code{x}. The smoothed data are then interpolated back to \code{x}. 

\code{cv.select} uses the standard R function \code{\link[stats]{optimize}} to optimize \code{\link{cv.score}}. If the argument \code{interval} is not specified, the lower bound of the search interval is by default \eqn{(x_2-x_1)/2}{(x[2] - x[1])/2} if \eqn{degree < 2} and \eqn{x_2-x_1}{x[2] - x[1]} if \eqn{degree >= 2}. The default value of the upper bound is \eqn{(\max(x)-\min(x))/2}{(max(x) - min(x))/2}. These values guarantee in most cases that the local polynomial estimator is well defined. It is often useful to plot the function to be optimized for a range of argument values (grid search) before applying a numerical optimizer. In this way, the search interval can be narrowed down and the optimizer is more likely to find a global solution.  


}


\value{
a bandwidth that minimizes the cross-validation score. 
}

\references{
Rice, J. A. and Silverman, B. W. (1991). Estimating the mean and covariance structure nonparametrically when the data are curves. \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, \bold{53}, 233--243.
}



\seealso{
\code{\link{cv.score}}, \code{\link{plugin.select}}
}

\examples{
\dontrun{
## Plasma citrate data
## Compare cross-validation scores and bandwidths  
## for local linear and local quadratic smoothing

data(plasma)
time <- 8:21   				
## Local linear smoothing						
cv.select(time, plasma, 1)	# local solution h = 3.76, S(h) = 463.08			
cv.select(time, plasma, 1, interval = c(.5, 1))	# global solution = .75, S(h) = 439.54

## Local quadratic smoothing						
cv.select(time, plasma, 2)	# global solution h = 1.15, S(h) = 432.75			
cv.select(time, plasma, 2, interval = c(1, 1.5))	# same
}
}


