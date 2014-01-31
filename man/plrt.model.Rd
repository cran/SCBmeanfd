\name{plrt.model}
\alias{plrt.model}

\title{
Pseudo-Likelihood Ratio Test for Models of the Mean Function  
}

\description{
Implement the Pseudo-Likelihood Ratio Test (PLRT) of Azzalini and Bowman (1993) with functional data. The test is used to assess whether a mean function belongs to a given finite-dimensional function space.}

\usage{ plrt.model(x, y, model, verbose = FALSE) }

\arguments{
  \item{x}{a numeric vector of x data. \code{x} must be a uniform grid; missing values are not accepted.
}
  \item{y}{
a matrix or data frame with functional observations (= curves) stored in rows. The number of columns of \code{y} must match the length of \code{x}. Missing values are not accepted.}

  \item{model}{an integer specifying the degree of a polynomial basis, or a data frame/matrix containing the basis functions stored in columns. In the latter case, the basis functions must be evaluated at \code{x}.
}

  \item{verbose}{
logical; if \code{TRUE}, information on the candidate model, kernel bandwidth, test statistic, and p value is printed; otherwise no print.}
}

\value{
the p value of the PLRT.
}

\references{
Azzalini, A. and Bowman, A. (1993). On the use of nonparametric regression for checking linear relationships. \emph{Journal of the Royal Statistical Society. Series B} \bold{55}, 549-557.
}


\seealso{
\code{\link{scb.model}}
}

\examples{
## Example: Gaussian process with mean = linear function + bump 
## and Onstein-Uhlenbeck covariance. The bump is high in the y   
## direction and narrow in the x direction. The SCB and PLRT 
## tests are compared.

# The departure from linearity in the mean function is strong  
# in the supremum norm (SCB test) but mild in the euclidean norm 
# (PLRT). With either n = 20 or n = 100 curves, the SCB test 
# strongly rejects the incorrect linear model for the mean 
# function while the PLRT retains it. 

p    <- 100    # number of observation points
x 	  <- seq(0, 1, len = p)
mu <- -1 + 1.5 * x + 0.2 * dnorm(x, .6, .02) 
plot(x, mu, type = "l")
R 	  <- (.25)^2 * exp(20 * log(.9) * abs(outer(x,x,"-"))) # covariance 
eigR <- eigen(R, symmetric = TRUE)  	
simR <- eigR$vectors \%*\% diag(sqrt(eigR$values)) 	 

n  <- 20
set.seed(100)
y  <- mu + simR \%*\% matrix(rnorm(n*p), p, n) 		
y  <- t(y) 	
points(x, colMeans(y))									
h  <- cv.select(x, y, 1)	
scb.model(x, y, 1, bandwidth = h)		 # p value: <1e-16
plrt.model(x, y, 1, verbose = TRUE)	# p value: .442
n  <- 100
y  <- mu + simR \%*\% matrix(rnorm(n*p), p, n) 		
y  <- t(y) 	
h  <- cv.select(x, y, 1)	
scb.model(x, y, 1, bandwidth = h)		 # p value: <1e-16
plrt.model(x, y, 1, verbose = TRUE)	# p value: .456
}

