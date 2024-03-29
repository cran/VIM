% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medianSamp.R
\name{medianSamp}
\alias{medianSamp}
\title{Aggregation function for a ordinal variable}
\usage{
medianSamp(x, weights = NULL)
}
\arguments{
\item{x}{ordered factor vector}

\item{weights}{numeric vector providing weights for the observations in x}
}
\description{
The function medianSamp chooses the level as the median or randomly between
two levels.
}
\seealso{
Other imputation methods: 
\code{\link{hotdeck}()},
\code{\link{impPCA}()},
\code{\link{irmi}()},
\code{\link{kNN}()},
\code{\link{matchImpute}()},
\code{\link{rangerImpute}()},
\code{\link{regressionImp}()},
\code{\link{sampleCat}()}
}
\concept{imputation methods}
