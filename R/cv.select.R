cv.select <- function(x, y, degree = 1, interval = NULL, gridsize = length(x), ...) 
{ 		
	if (is.null(interval)) {
		rangex   <- diff(range(x))
		meshx    <- max(diff(sort(x)))
		interval <- c(ifelse(degree<2, meshx/2, meshx), rangex/2)
	}
	
	optimize(cv.score, interval, x=x, y=y, degree=degree, 
		gridsize=gridsize, ...)$minimum 
}
