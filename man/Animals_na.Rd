% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VIM-package.R
\docType{data}
\name{Animals_na}
\alias{Animals_na}
\title{Animals_na}
\format{
A data frame with 28 observations on the following 2 variables.
\describe{
\item{lbody}{log body weight}
\item{lbrain}{log brain weight}
}
}
\source{
P. J. Rousseeuw and A. M. Leroy (1987) Robust Regression and Outlier Detection. Wiley, p. 57.
}
\description{
Average log brain and log body weights for 28 Species
}
\details{
The original data can be found in package MASS.
10 values on brain weight are set to be missing.
}
\examples{

data(Animals_na)
aggr(Animals_na)

}
\references{
Venables, W. N. and Ripley, B. D. (1999) Modern Applied Statistics with S-PLUS. Third Edition. Springer.

Templ, M. (2022) Visualization and Imputation of Missing Values. Springer Publishing. Upcoming book.
}
\keyword{datasets}
