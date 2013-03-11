\name{newAggregation}
\alias{newAggregation}
\title{Recalculate new site aggregations from a previous MCMC aggregation}
\usage{
  newAggregation(fit, aggregation.data)
}
\arguments{
  \item{fit}{The output list from a previous call to
  \code{\link{mcmc.aggregate}}. In order to use this
  function, \code{keep.site.abund = TRUE} had to be used in
  the original creation of \code{fit}. Else, there is
  nothing to be aggregated!}

  \item{aggregation.data}{A data frame with the sites in
  one column (with the same name as \code{site.name} used
  in the original call to create \code{fit}). The other
  columns are factor variables defining other site
  aggregations.}
}
\value{
  A named list with names equal to the variables in
  \code{aggregation.data}. Each element of the list is
  another list with elements: \item{aggregated.abund}{The
  MCMC sample of the new aggregation}
  \item{aggregation.summary}{A summary of the aggregation
  MCMC}
}
\description{
  If the site abundence sample was retained in a call to
  \code{\link{mcmc.aggregate}} the sites can be
  re-aggregated according to different region
  specifications.
}

