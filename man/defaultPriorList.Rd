\name{defaultPriorList}
\alias{defaultPriorList}
\title{Create list of prior parameters for \code{mcmc.aggregate(...)} trend estimation.}
\usage{
defaultPriorList(trend.limit = NULL, model.data, gamma.mean, gamma.prec)
}
\arguments{
  \item{trend.limit}{Bounds on the individual site trends.}

  \item{model.data}{A data set which provides the
  individual site models used in
  \code{\link{mcmc.aggregate}}.}

  \item{gamma.mean}{Prior mean of the (multivariate) normal
  prior distribution for gamma}

  \item{gamma.prec}{Precision (inverse of variance) of the
  gamma prior distribution}
}
\value{
A named list with entries for each of the parameters for
which prior specification is necessary.
}
\description{
This function takes user input to create the list of prior
paremters necessary for Bayesian inference of aggregated
trends. The prior parameter list is a necessary argument
for the \code{\link{mcmc.aggregate}} function.
}
\details{
Using the site model data set \code{model.data} and a soft
bound on the individual site trends (\code{trend.limit})
this function creates a sensible default for the list of
prior parameters necessary for Bayesian inference of
aggregated trends. The value of \code{trend.limit} should
be in terms of % growth (e.g., \code{trend.limit=0.2}
implies an upper limit of 20 percent growth at each
individual site). The lower bound is calculated as
1/(1+trend.limit). This is the ower bound that is symmetric
about zero on the log-scale. For \code{trend.limit=0.2},
the lower bound is approximately -0.17.

Specifically, the prior distribution for the linear trend
parameters (beta) in log abundance at each site is set to a
bivariate normal with mean \code{m.i=c(0,0)} and precision
matrix \code{Q.i=diag(c(0,q))}, where \code{q} is chosen
such that Pr[1/(1+trend.limit) < exp(beta[2])-1 <
trend.limit] = 0.95.
}

