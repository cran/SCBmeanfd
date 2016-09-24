\name{plugin.select}
\alias{plugin.select}

\title{
Plug-in Bandwidth Selection for Local Polynomial Estimation}

\description{
Select the plug-in bandwidth described in Benhenni and Degras (2011) for the local polynomial estimation of a mean function and its first derivative based on functional data. 
}

\usage{
plugin.select(x, y, drv = 0L, degree = drv+1L, gridsize = length(x), ...) 
}

\arguments{
  \item{x}{
numeric vector of x data. This observation grid must be uniform and missing values are not accepted.}
  \item{y}{
matrix or data frame with functional observations (= curves) stored in rows. The number of columns of \code{y} must match the length of \code{x}. Missing values are not accepted.}  
\item{drv}{order of the derivative to estimate. Must be 0 or 1. 
}
  \item{degree}{degree of local polynomial used. Must equal \code{drv} or \code{drv+1}.
}
 \item{gridsize}{
the size of the grid over which the mean function is to be estimated.  
Defaults to \code{length(x)}.
}
\item{...}{additional arguments to pass to the optimizer of the CV score.}
}

\details{
The plug-in method should not be used with small data sets, since it is based on asymptotic considerations and requires reasonably accurate estimates of derivatives of the mean and covariance functions. Both the number of observed curves and observation points should be moderate to large. The plug-in bandwidth is designed to minimize the asymptotic mean integrated  squared estimation error \deqn{AMISE(h) = \int (\mu(t) - \hat{\mu}(t; h))^2 dt,}{AMISE(h) = integral of (\mu(t) - \mu.hat(t;h))^2,} 
where \eqn{\mu(x)}{\mu(t)} is the mean function and \eqn{\hat{\mu}(t;h)}{mu.hat(t;h)} is a local polynomial estimator with kernel bandwidth \eqn{h}. The expression of the plug-in bandwidth can be found in Benhenni and Degras (2011).}

\value{
the plug-in bandwidth.}

\references{
Benhenni, K. and Degras, D. (2011). Local polynomial estimation of the average growth curve with functional data. \url{http://arxiv.org/abs/1107.4058}}

\seealso{
\code{\link{cv.select}}
}

\examples{
\dontrun{
## Phoneme data
data(phoneme)
classes <- phoneme[,151]
phoneme <- phoneme[,-151]
freq    <- 1:150
plugin.bandwidth  <- numeric(5) 
cv.bandwidth  <- numeric(5)  # compare with cross-validation
for (i in 1:5) {
  plugin.bandwidth[i] <- plugin.select(x = freq, y = phoneme[classes == i, ],
                          drv = 0, degree = 1) 
  cv.bandwidth[i]     <- cv.select(x = freq, y = phoneme[classes == i, ], 
                          degree = 1)
}

round(cbind(plugin.bandwidth, cv.bandwidth), 4)

}
}