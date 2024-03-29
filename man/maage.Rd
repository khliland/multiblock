% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maage.R
\name{maage}
\alias{maage}
\alias{maageSeq}
\title{Måge plot}
\usage{
maage(
  object,
  expl_var = TRUE,
  pure.trace = FALSE,
  pch = 20,
  xlab = "# components",
  ylab = ifelse(expl_var, "Explained variance (\%)", "RMSECV"),
  xlim = NULL,
  ylim = NULL,
  cex.text = 0.8,
  ...
)

maageSeq(
  object,
  compSeq = TRUE,
  expl_var = TRUE,
  pch = 20,
  xlab = "# components",
  ylab = ifelse(expl_var, "Explained variance (\%)", "RMSECV"),
  xlim = NULL,
  ylim = NULL,
  cex.text = 0.8,
  col = "gray",
  col.block = c("red", "blue", "darkgreen", "purple", "black", "red", "blue",
    "darkgreen"),
  ...
)
}
\arguments{
\item{object}{An SO-PLS model (\code{sopls} object)}

\item{expl_var}{Logical indicating if explained variance (default) or RMSECV should be displayed.}

\item{pure.trace}{Logical indicating if single block solutions should be traced in the plot.}

\item{pch}{Scalar or symbol giving plot symbol.}

\item{xlab}{Label for x-axis.}

\item{ylab}{Label for y-axis.}

\item{xlim}{Plot limits for x-axis (numeric vector).}

\item{ylim}{Plot limits for y-axis (numeric vector).}

\item{cex.text}{Text scaling (scalar) for better readability of plots.}

\item{...}{Additional arguments to \code{plot}.}

\item{compSeq}{Integer vector giving the sequence of previous components chosen for \code{maageSeq} (see example).}

\item{col}{Line colour in plot.}

\item{col.block}{Line colours for blocks (default = c('red','blue','darkgreen','purple','black'))}
}
\value{
The \code{maage} plot has no return.
}
\description{
Måge plot for SO-PLS (\code{\link{sopls}}) cross-validation visualisation.
}
\details{
This function can either be used
for global optimisation across blocks or sequential optimisation, using \code{maageSeq}.
The examples below show typical usage.
}
\examples{
data(wine)
ncomp <- unlist(lapply(wine, ncol))[-5]
so.wine <- sopls(`Global quality` ~ ., data=wine, ncomp=ncomp, 
            max_comps=10, validation="CV", segments=10)
maage(so.wine)

# Sequential search for optimal number of components per block
old.par <- par(mfrow=c(2,2), mar=c(3,3,0.5,1), mgp=c(2,0.7,0))
maageSeq(so.wine)
maageSeq(so.wine, 2)
maageSeq(so.wine, c(2,1))
maageSeq(so.wine, c(2,1,1))
par(old.par)
}
\seealso{
Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
}
