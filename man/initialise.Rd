\name{initialise}
\alias{initialise}
\title{Initialization of missing values}
\description{
  Rough estimation of missing values in a vector according to its type.
}
\usage{
initialise(x,mixed,method="kNN")
}
\arguments{
  \item{x}{a vector.}
  \item{mixed}{a character vector containing the names of variables of type mixed (semi-continous).}
  \item{method}{Method used for Initialization (median or kNN)} 
}
\details{
  Missing values are imputed with the mean for vectors of class 
  \code{"numeric"}, with the median for vectors of class \code{"integer"}, 
  and with the mode for vectors of class \code{"factor"}.  Hence, \code{x} 
  should be prepared in the following way: assign class \code{"numeric"} 
  to numeric vectors, assign class \code{"integer"} to ordinal  vectors, 
  and assign class \code{"factor"} to nominal or binary  vectors. 
}
\value{
  the initialized vector.
}
\author{Matthias Templ, modifications by Andreas Alfons}
\note{  
  The function is used internally by some imputation algorithms.
}
\keyword{manip}
