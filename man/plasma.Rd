\name{plasma}
\alias{plasma}
\docType{data}
\title{
Plasma citrate data}
\description{
Plasma citrate concentrations measured on 10 human subjects on the same day. The measurements on an individual were taken each hour from 8:00AM to 9:00PM. A possible statistical analysis is to estimate the population mean plasma citrate concentration as a function of time of day.}

\usage{data(plasma)}
\format{
A matrix with 10 rows (corresponding to subjects) and 14 columns (corresponding to hours).
}

\source{
Andersen, A. H., Jensen, E. B. and Schou, G. (1981). Two-way analysis of variance with correlated errors. \emph{International Statistical Review} \bold{49}, 153--157.

Hart, T. D. and Wehrly, T. E. (1986). Kernel regression estimation using repeated measurements data. \emph{Journal of the American Statistical Association} \bold{81}, 1080--1088.
}

\examples{
\dontrun{
data(plasma)
matplot(x = 8:21, y = t(plasma), cex = .75, type = "l", col = 1, lty = 1,
  lwd = .5, xlab = "Time of day (hour)", ylab = "Concentration", 
  main = "Plasma citrate data")
lines(8:21, colMeans(plasma), col = 2, lwd = 1.5)
legend("topright", col = 2, lty = 1, lwd = 1.5, legend = "Mean")}
}
\keyword{datasets}
