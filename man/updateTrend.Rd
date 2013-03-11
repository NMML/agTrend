\name{updateTrend}
\alias{updateTrend}
\title{Recalculate posterior predictive trend with a different time scale.}
\usage{
  updateTrend(start, end, x)
}
\arguments{
  \item{start}{A new start value for the time span}

  \item{end}{A new end value for the time span}

  \item{x}{An mcmc augmentation object produced by a call
  to \code{\link{mcmc.aggregate}} or an element of the list
  produced by a call to \code{\link{newAggregation}}.}
}
\description{
  Recalculate posterior predictive trend with a different
  time scale.
}

