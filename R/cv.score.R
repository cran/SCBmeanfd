cv.score <- function(bandwidth, x, y, degree=1, gridsize = length(x))
{
	n <- NROW(y)
	N <- NCOL(y)
	stopifnot(n > 1)
	stopifnot(length(x) == N)
	xgrid <- seq(min(x),max(x),len=gridsize)
	interpolation.flag <- (gridsize != N && !isTRUE(all.equal(x,xgrid)))
	
	y.hat <- if (interpolation.flag) {
		apply(y, 1, function(z) approx(xgrid, locpoly(x = x, y = z, 
			bandwidth = bandwidth, gridsize = gridsize)$y, x)$y)
	} else {
		apply(y, 1, function(z) locpoly(x = x, y = z, 
			bandwidth = bandwidth, gridsize = gridsize)$y)
	}

	mu.hat <- rowMeans(y.hat)
	residuals <- (n/(n-1)) * mu.hat - y.hat / (n-1) - t(y)
	mean(residuals^2)
}
