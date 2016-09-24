plugin.select <- function(x, y, drv = 0L, degree = drv+1L, gridsize = length(x), ...) 
{	

	stopifnot(NCOL(y) == length(x))
	if (!is.matrix(y)) y <- as.matrix(y)
	n <- nrow(y)
	N <- ncol(y)
	stopifnot(drv %in% c(0,1)) 
	stopifnot((degree-drv) %in% c(0,1)) 
	spacing <- diff(x)
	if (any(spacing < 0) || !isTRUE(all.equal(min(spacing),max(spacing))))
		stop("x must be a uniform grid") 
	range.x <- x[N] - x[1]
	int.K1 	<- 1.1283
	int.K2 	<- 0.5641
	ybar  <- colMeans(y)
	I.alpha <- mean(apply(y[,-1]-y[,-N], 1, crossprod))
	h.cv 	<- cv.select(x, y, degree, gridsize = gridsize, ...)
	mu.cv 	<- locpoly(x, ybar, degree = degree, bandwidth = h.cv, gridsize = N)$y
	Dmu.cv  <- diff(mu.cv, differences = drv + 2) * (N-1)^(drv+2)
	L2 		<- sum(Dmu.cv^2) / (N-1)
	if (drv == 0) h.plug <- max((I.alpha * int.K1 / 2 / L2 / n)^(1/3), .5/(N-1))
	if (drv == 1) h.plug <- max((I.alpha * int.K2 / 2 / L2 / n)^(1/5), .5/(N-1))	
	return(h.plug * range.x)
}
