scb.equal <- function(x, y, bandwidth, level = .05, degree = 1, 
	scbtype = c("normal","bootstrap","both","no"), gridsize = NULL, 
	keep.y = TRUE, nrep = 2e4, nboot = 1e4, parallel = c("no","multicore","snow"), 
	ncpus = getOption("boot.ncpus",1L), cl = NULL)
{
	caLL <- match.call()
	if (any(is.na(unlist(x))) || any(is.na(unlist(y))))
		stop("x and/or y contain missing values")
	
	if (is.list(x)) {
		x1 <- x[[1]] ; x2 <- x[[2]]
	} else if (NCOL(x) == 2) {
		x1 <- x[,1] ; x2 <- x[,2]
	} else if (NCOL(x) == 1) {
		x1 <- x2 <- x 
	}		
	if (x1[1] != x2[1] || x1[length(x1)] != x2[length(x2)])
		stop("The range of x is not the same in the two samples")

	y1 <- y[[1]]
	y2 <- y[[2]]
	if (length(x1) != ncol(y1) || length(x2) != ncol(y2))
		stop("The dimensions of x and y do not match in at least one sample") 
	if (is.null(gridsize))
		gridsize <- min(length(x1),length(x2))
		
	y1.hat <- apply(y1, 1, function(z) locpoly(x1, z, degree = degree, 
			  		bandwidth = bandwidth[1], gridsize = gridsize)$y)
	y2.hat <- apply(y2, 1, function(z) locpoly(x2, z, degree = degree, 
			 		bandwidth = bandwidth[2], gridsize = gridsize)$y)
	
	mu1.hat <- rowMeans(y1.hat)
	mu2.hat <- rowMeans(y2.hat)
	n1 <- nrow(y1)
	n2 <- nrow(y2)		
	R1.hat <- cov(t(y1.hat)) / n1
	R2.hat <- cov(t(y2.hat)) / n2
	R.hat <- R1.hat + R2.hat
	se <- sqrt(diag(R.hat))
	test.stat <- max(abs(mu1.hat - mu2.hat) / se)

	p.norm = p.boot = NULL
	q.norm = q.boot = NULL
	lb.norm = ub.norm = lb.boot = ub.boot = NULL
	scbtype <- match.arg(scbtype)	

	if (scbtype %in% c("normal","both")) {
		eigcoR 	<- eigen(cov2cor(R.hat), TRUE)
		ncomp 	<- which.max(cumsum(eigcoR$values) > .99 * sum(eigcoR$values))
		vars 	<- matrix(rnorm(ncomp * nrep), ncomp, nrep)
		M <- eigcoR$vectors[,1:ncomp] %*% 
			diag(sqrt(eigcoR$values[1:ncomp]),ncomp,ncomp)
		supnorm <- apply(abs(M %*% vars), 2, max)
		p.norm	<- 1 - ecdf(supnorm)(test.stat)
		q.norm	<- as.numeric(quantile(supnorm, 1-level)) 
		lb.norm	<- mu1.hat - q.norm * se
		ub.norm	<- mu1.hat + q.norm * se	
	}

	if (scbtype %in% c("bootstrap","both")) {
		r <- cbind(y1.hat - mu1.hat, y2.hat - mu2.hat)
		ix1 <- 1:n1
		ix2 <- (n1+1):(n1+n2)
		boot.stat <- function(mat, ix) {	
			mat1 <- mat[ix[ix1],]
			mat2 <- mat[ix[ix2],]
			mu1.boot <- colMeans(mat1)
			mu2.boot <- colMeans(mat2)
			var.mu12.boot <- (colSums(mat1^2)-n1*mu1.boot^2)/(n1*(n1-1)) + 
				(colSums(mat2^2)-n2*mu2.boot^2)/(n2*(n2-1))
			sqrt(max((mu1.boot-mu2.boot)^2 / var.mu12.boot))
		}
		supnorm <- boot(t(r), boot.stat, nboot, strata = rep(1:2, c(n1, n2)), 
			parallel = parallel, ncpus = ncpus, cl = cl)$t
		p.boot 	<- 1 - ecdf(supnorm)(test.stat)
		q.boot 	<- as.numeric(quantile(supnorm, 1 - level))
		lb.boot <- mu1.hat - q.boot * se
		ub.boot <- mu1.hat + q.boot * se	
	}

	result <- list( x = list(x1, x2), y = if(keep.y) y else NULL, call = caLL, model = NULL, par = NULL, nonpar = cbind(mu1.hat, mu2.hat), bandwidth = bandwidth, degree = degree, level = level, scbtype = scbtype, teststat = test.stat, pnorm = p.norm, pboot = p.boot, qnorm = q.norm, qboot = q.boot, normscb = cbind(lb.norm, ub.norm), bootscb = cbind(lb.boot, ub.boot), gridsize = gridsize, nrep = nrep, nboot = nboot )

	class(result) <- "SCBand"
	return(result)	
		
}
