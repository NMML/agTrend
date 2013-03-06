\name{mcmc.aggregate}
\alias{mcmc.aggregate}
\title{Posterior predictive sampling, aggregtion of abundance counts, and linear trend summary}
\usage{
  mcmc.aggregate(start, end, data, obs.formula = NULL,
    aggregation, model.data, abund.name, time.name,
    site.name, sig.abund, ln.adj = 0, upper = Inf,
    lower = -Inf, burn, iter, thin = 1, prior.list = NULL,
    keep.site.abund = FALSE, keep.site.param = FALSE,
    keep.obs.param = FALSE)
}
\arguments{
  \item{start}{}

  \item{end}{}

  \item{data}{A \code{data.frame} that contains the
  abundance survey data.}

  \item{obs.formula}{}

  \item{aggregation}{}

  \item{model.data}{}

  \item{abund.name}{}

  \item{time.name}{}

  \item{site.name}{}

  \item{sig.abund}{}

  \item{ln.adj}{}

  \item{upper}{}

  \item{lower}{}

  \item{burn}{}

  \item{iter}{}

  \item{thin}{}

  \item{prior.list}{}

  \item{keep.site.abund}{}

  \item{keep.site.param}{}

  \item{keep.obs.param}{}
}
\description{
  Function for sampling from the posterior predictive
  distribution of abundance (counts) at individual sites.
  Then aggregating the counts over the specified
  aggregation variable.
}
\details{
  This function is the workhorse for the \code{agTrend}
  package. It performs MCMC sampling of the posterior
  predictive distribution of the abundance at each site at
  each time, \eqn{N_{st}}{N_st}. The abundance at each site
  is modeled, in its most general form, with a
  zero-inflated, nonparameteric, log-linear model,
  \deqn{z_{st} = \beta_{s0} + \beta_{s1}t + \eta_{st} +
  \delta_{st} \mbox{ if } N_{st}>0,}{z_st = beta_s0 +
  beta_s1 * t + eta_st + delta_st if N_st > 0,} where
  \eqn{\beta_{s0}+\beta_{s1}t}{beta_s0 + beta_s1 * t} is
  the linear trend, \eqn{\eta}{eta} is a random walk of
  order 2 (RW2), and

  Then, for each iteration of the MCMC sampler, the
  complete data is aggregated (i.e., summed) over all sites
  within a specified region. Thus, one can obtain a
  posterior predictive sample from the aggregated abundance
  for every time between the first time in the data set and
  the last. In its most general form, the model fit to each
  site \eqn{s}{s} at time \eqn{t}{t} is specified by
  \deqn{y_{st} = x_{st}'\gamma + z_{st} +
  \epsilon_{st},}{y_st = x_st gamma + + eps_st,} where
  \eqn{y_{st}}{y_st} is the raw observed log abundance,
  \eqn{z_{st}}{z_st} is the 'corrected' standardized log
  abundance, \eqn{x_{st}}{x_st} is a vector of covariates
  used to standardize the raw abundance,
  \eqn{\gamma}{gamma} is a vector of coefficients, and
  \eqn{[\epsilon_{st}]=N(0,\sigma^2_{st})}{[eps_st]=N(0,sigma_st^2)}.
  Currently, \eqn{\sigma_{st}}{sigma_st} is considered to
  be known and is specified as a column in \code{data} by
  the \code{sig.abund} argument. The abundance
  standardization parameter is added to help account for
  changing survey methodology. The general

  Using the aggregations, the average linear trend is
  calculated between \code{start} and \code{end}. By using
  the posterior predictive distribution, we can account for
  parameter uncertainty at each site as well as sampling
  variability of the survey. even though we are using the
  Bayesian inference paradigm, we still capture the essence
  of frequentist inference by accounting for the
  'replication' of the survey effort by using the
  predictive distribution even for times and places where
  we have survey data.

  The data is specified by the
}

