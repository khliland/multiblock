% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{explvar}
\alias{explvar}
\title{Explained predictor variance}
\usage{
explvar(object)
}
\arguments{
\item{object}{An object fitted using a method from the multiblock package}
}
\value{
A \code{vector} of component-wise explained variances for predictors.
}
\description{
Extraction and/or computation of explained variances for various
object classes in the \code{multiblock} package.
}
\examples{
data(potato)
so <- sopls(Sensory ~ Chemical + Compression, data=potato, ncomp=c(10,10), 
            max_comps=10)
explvar(so)
}
