% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supervised.R
\name{supervised}
\alias{supervised}
\title{Supervised Multiblock Methods}
\description{
Collection of supervised multiblock methods:
\itemize{
\item MB-PLS - Multiblock Partial Least Squares (\code{\link{mbpls}})
\item sMB-PLS - Sparse Multiblock Partial Least Squares (\code{\link{smbpls}})
\item SO-PLS - Sequential and Orthogonalized PLS (\code{\link{sopls}})
\item PO-PLS - Parallel and Orthogonalized PLS (\code{\link{popls}})
\item ROSA - Response Oriented Sequential Alternation (\code{\link{rosa}})
\item mbRDA - Multiblock Redundancy Analysis (\code{\link{mbrda}})
}
}
\examples{
data(potato)
mb <- mbpls(Sensory ~ Chemical + Compression, data=potato, ncomp = 5)
print(mb)

# Convert data.frame with AsIs objects to list of matrices
potatoList <- lapply(potato, unclass)
mbr <- mbrda(Sensory ~ Chemical + Compression, data=potatoList, ncomp = 10)
print(mbr)
scoreplot(mbr, labels="names")

}
\seealso{
Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
}
