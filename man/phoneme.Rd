\name{phoneme}
\alias{phoneme}
\docType{data}
\title{
Phoneme data}
\description{
Log-periodograms of five phonemes ("sh" as in "she", "dcl" as in "dark", "iy" as the vowel in "she", "aa" as the vowel in "dark", and "ao" as the first vowel in "water") spoken by a sample of ~50 male speakers. For each phoneme, 400 log-periodograms are observed on a uniform grid of 150 frequencies.}

\usage{data(phoneme)}
\format{
A data frame with 2000 rows and 151 columns.
Each row contains a discretized log-periodogram (150 first columns) followed by a phoneme indicator (last column) coded as follows:
\tabular{cc}{
Code \tab Phoneme \cr
1 \tab "sh" \cr
2 \tab "iy" \cr
3 \tab "dcl" \cr
4 \tab "aa" \cr
5 \tab "ao"}
}



\source{
This dataset was created by the STAPH group in Toulouse, France\cr
\url{http://www.math.univ-toulouse.fr/staph/npfda/}.\cr
The larger original dataset can be found at\cr
\url{http://www-stat.stanford.edu/ElemStatLearn}.
}

\examples{
\dontrun{
data(phoneme)
freq <- 1:150
classes <- phoneme[,151]
phoneme <- phoneme[,-151]
classnames <- c("sh", "iy", "dcl", "aa", "ao")

## Local linear fit to the mean log-periodogram for each phoneme 
llfit <- mapply(locpoly, y = by(phoneme, classes, colMeans), 
          MoreArgs = list(x = freq, bandwidth = 2, degree = 1, 
          gridsize = length(freq)))
llfit.y <- matrix(unlist(llfit["y",]), 150, 5)  
matplot(freq, llfit.y, type = "l", lty = 1, xlab = "Frequency (scaled)", 
	ylab = "Log-intensity", 
	main = "Local linear estimation\nof the population mean log-periodogram")
legend("topright", legend = classnames, col = 1:5, lty = 1)
}
}

\keyword{datasets}
