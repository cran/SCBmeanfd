\name{cv.score}
\alias{cv.score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leave-One-Curve-out Cross-Validation Score
}
\description{
Compute the cross-validation score of Rice and Silverman (1991) for  the local polynomial estimation of a mean function.}
\usage{
cv.score(bandwidth, x, y, degree = 1, gridsize = length(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{bandwidth}{kernel bandwidth.}
\item{x}{observation points. Missing values are not accepted.}
\item{y}{matrix or data frame with functional observations (= curves) stored in rows. The number of columns of \code{y} must match the length of \code{x}. Missing values are not accepted.}
\item{degree}{degree of the local polynomial fit.}
\item{gridsize}{size of evaluation grid for the smoothed data.}  
}

\details{
The cross-validation score is obtained by leaving in turn each curve out and computing the prediction error of the local polynomial smoother based on all other curves. For a bandwith value \eqn{h}, this score is 
\deqn{ CV(h) = \frac{1}{np} \sum_{i=1}^n \sum_{j=1}^p \left( Y_{ij} - \hat{\mu}^{-(i)}(x_j;h) \right)^2, }{CV(h) =  \sum (Y[ij] - \mu^{-i}(x[j];h))^2 / (n*p)}
where \eqn{Y_{ij}}{Y[ij]} is the measurement of the \eqn{i}-th curve at location \eqn{x_j}{x[j]} for \eqn{i=1,\ldots,n} and \eqn{j=1,\ldots,p}, and \eqn{\hat{\mu}^{-(i)}(x_j;h)}{\mu^{-i}(x[j];h)} is the local polynomial estimator with bandwidth \eqn{h} based on all curves except the \eqn{i}-th.

If the \code{x} values are not equally spaced, the data are first smoothed and evaluated on a grid of length \code{gridsize} spanning the range of \code{x}. The smoothed data are then interpolated back to \code{x}. 
}

\value{the cross-validation score.}

\references{
Rice, J. A. and Silverman, B. W. (1991). Estimating the mean and covariance structure nonparametrically when the data are curves. \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, \bold{53}, 233--243.
}



\seealso{
\code{\link{cv.select}}
}
\examples{
## Artificial example 
x <- seq(0, 1, len = 100)
mu <- x + .2 * sin(2 * pi * x)
y <- matrix(mu + rnorm(2000, sd = .25), 20, 100, byrow = TRUE)
h <- c(.005, .01, .02, .05, .1, .15)
cv <- numeric()
for (i in 1:length(h)) cv[i] <- cv.score(h[i], x, y, 1)
plot(h, cv, type = "l")

## Plasma citrate data
## Compare cross-validation scores and bandwidths  
## for local linear and local quadratic smoothing
\dontrun{
data(plasma)
time <- 8:21
h1 <- seq(.5, 1.3, .05)
h2 <- seq(.75, 2, .05)
cv1 <- sapply(h1, cv.score, x = time, y = plasma, degree = 1)
cv2 <- sapply(h2, cv.score, x = time, y = plasma, degree = 2)
plot(h1, cv1, type = "l", xlim = range(c(h1,h2)), ylim = range(c(cv1, cv2)), 
  xlab = "Bandwidth (hour)", ylab = "CV score", 
  main = "Cross validation for local polynomial estimation")
lines(h2, cv2, col = 2)
legend("topleft", legend = c("Linear", "Quadratic"), lty = 1, 
  col = 1:2, cex = .9)

## Note: using local linear (resp. quadratic) smoothing 
## with a bandwidth smaller than .5 (resp. .75) can result 
## in non-definiteness or numerical instability of the estimator. 
}
}
