<!-- README.md is generated from README.Rmd. Please edit that file -->
[![DOI](https://zenodo.org/badge/8592075.svg)](https://zenodo.org/badge/latestdoi/8592075)

Fit regional trends to site-specific abundence data
---------------------------------------------------

This package fits a log-linear trend models to regions aggregated over sites. The sites may contain missing surveys that are not temporally aligned with the missing data at other sites, making direct aggregation impossible. The functions within the package model the indivdual sites with a semi-parametric (possibly, zero-inflated) model to interpolate missing data from which regional aggregations can be made. By using Markov Chain Monte Carlo, on can sample from the posterior predictive distribution of the regional aggregations Then calculate the log-linear trend over the time period of interest as a derived parameter. Using the posterior predictive distribution allows incorporation of both parameter uncertainty as well as uncertainty due to sampling the local abundance processes.

### Disclaimer

*This software package is developed and maintained by scientists at the NOAA Fisheries Alaska Fisheries Science Center and should be considered a fundamental research communication. The reccomendations and conclusions presented here are those of the authors and this software should not be construed as official communication by NMFS, NOAA, or the U.S. Dept. of Commerce. In addition, reference to trade names does not imply endorsement by the National Marine Fisheries Service, NOAA. While the best efforts have been made to insure the highest quality, tools such as this are under constant development and are subject to change.*

### Example

Load packages for this example

``` r
library(agTrend)
#> Loading required package: coda
#> Loading required package: Matrix
#> Loading required package: mgcv
#> Loading required package: nlme
#> This is mgcv 1.8-23. For overview type 'help("mgcv-package")'.
#> Loading required package: truncnorm
#> 
#>  agTrend 0.17.7 (2017-03-23) 
#>  A demo is available at https://github.com/NMML/agTrend
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:nlme':
#> 
#>     collapse
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

Now we can load the data that is included with the `agTrend` package and filter it to the data we want for this example (i.e., 1990-2016).

``` r
data(wdpsNonpups)
wdpsNonpups %>% filter(YEAR>=1990) %>% droplevels() -> wdpsNonpups
```

Now, we count the number of positive counts at each site so that we can remove sites that had only 1 positive count

``` r
wdpsNonpups %>% group_by(SITE) %>% 
  summarise(num_counts=sum(COUNT>0)) %>% ungroup() -> nz_counts
wdpsNonpups %>% left_join(nz_counts) %>% filter(num_counts>1) %>% select(-num_counts) %>% 
  mutate(SITE=factor(SITE)) -> wdpsNonpups
#> Joining, by = "SITE"
```

Now we'll add a photo method covariate to data (oblique photos prior to 2004 surveys = 1)

``` r
wdpsNonpups %>% mutate(obl = as.integer(YEAR<2004)) %>% as.data.frame() -> wdpsNonpups
```

The next step is to create prediction and availability models for each site based on the number of surveys. The trend models will be

-   0-5 nonzero counts = constant trend signal
-   6-10 nonzero counts = linear trend signal
-   10 nonzero counts = RW2 (i.e., spline) trend signal

The zero inflation (avail) models are

-   0-5 surveys = constant inflation effect
-   5 surveys = linear inflation effect

-   All surveys have nonzero counts = no availability model

``` r
wdpsModels = wdpsNonpups %>% group_by(SITE) %>% 
  summarize(
    num_surv=n(),
    nz_count=sum(COUNT>0)) %>% 
  ungroup() %>%  
  mutate(
    trend = as.character(cut(nz_count, c(0,5,10,30), labels=c("const","lin","RW"))),
    avail = as.character(cut(num_surv, c(0,5,30), labels=c("const","lin")))
  ) %>% 
  mutate(avail = if_else(num_surv==nz_count, "none", avail)) %>% 
  select(-num_surv, -nz_count) %>% as.data.frame()

head(wdpsModels)
#>                 SITE trend avail
#> 1               ADAK    RW  none
#> 2 ADAK/ARGONNE POINT   lin   lin
#> 3  ADAK/CRONE ISLAND   lin   lin
#> 4             ADUGAK    RW  none
#> 5 AFOGNAK/TONKI CAPE   lin   lin
#> 6             AGATTU    RW  none
```

The next step involves creating a prior distribution list for MCMC site updating. The prior distribtuions for the trend parameters will enforce the assumption that site trends unlikely to be greater than 20% or less that about -17%. More exactly, for *r* = 0.2, *P**r**e**x**p*(−*l**n*(1 + *r*)) − 1 &lt; trend &lt; *r* = 0.95. In addition, we create an informative gamma prior for the survey methodology correction.

``` r
data("photoCorrection")
photoCorrection %>% mutate(log_ratio = OBLIQUE/VERTICAL) -> photoCorrection
gamma.0 = photoCorrection %>% summarize(mean(log_ratio)) %>% as.double()
gamma.se = photoCorrection %>% summarize(sd(log_ratio)/sqrt(n())) %>% as.double()

prior.list=defaultPriorList(trend.limit=0.2, model.data=wdpsModels,
                            gamma.mean=gamma.0, gamma.prec=1/gamma.se^2)
```

Finally, before we run the MCMC, the last step is to create a data set of upper bounds for predictive counts. For this bound, we use 3 times the maximum count observed at the site.

``` r
upper = wdpsNonpups %>% group_by(SITE) %>% summarize(upper=3*max(COUNT)) %>% ungroup()
```

Now we begin the MCMC sampling using the `mcmc.aggregate` function. This function performs the site augmentation and samples from the posterior predictive distribution of the count data. **Note: Only a small number of MCMC iterations are shown here. For a more robust analysis change to burn=1000 and iter=5000.**

``` r
set.seed(123) 
fit <- mcmc.aggregate(start=1990, end=2026, data=wdpsNonpups, obs.formula=~obl-1, model.data=wdpsModels, 
                      rw.order=list(omega=2), aggregation="REGION",
                      abund.name="COUNT", time.name="YEAR", site.name="SITE", forecast = TRUE,
                      burn=50, iter=100, thin=5, prior.list=prior.list, upper=upper, 
                      keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)
```

Let's look at the results

``` r
fitdat <- fit$aggregation.pred.summary
```

Compute trends for just 2000-2016

``` r
trend2006 = updateTrend(fit, 2003, 2016, "pred")
```

Change to percent growth form

``` r
growth2006 = mcmc(100*(exp(trend2006[,7:12])-1))
summary(growth2006)
#> 
#> Iterations = 1:100
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 100 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                  Mean     SD Naive SE Time-series SE
#> C ALEU (Trend)  2.180 0.9093  0.09093        0.23691
#> C GULF (Trend)  6.563 0.6461  0.06461        0.08083
#> E ALEU (Trend)  3.601 0.7187  0.07187        0.09984
#> E GULF (Trend)  7.453 1.5409  0.15409        0.20988
#> W ALEU (Trend) -2.663 1.7765  0.17765        0.17765
#> W GULF (Trend)  5.381 0.7341  0.07341        0.09077
#> 
#> 2. Quantiles for each variable:
#> 
#>                   2.5%    25%    50%    75%   97.5%
#> C ALEU (Trend)  0.4052  1.597  2.160  2.650  4.2380
#> C GULF (Trend)  5.3525  6.097  6.518  6.952  7.9334
#> E ALEU (Trend)  2.4247  3.076  3.547  4.007  5.1194
#> E GULF (Trend)  4.6746  6.385  7.363  8.413 10.6180
#> W ALEU (Trend) -5.8824 -3.685 -2.682 -1.758  0.7968
#> W GULF (Trend)  3.8947  4.904  5.420  5.872  6.6947
# Obtain posterior median % growth and 90% credible interval
print(
  data.frame(
    post.median=round(apply(growth2006, 2, median),2),
    HPD.90=round(HPDinterval(growth2006, 0.95),2)
  )
)
#>                post.median HPD.90.lower HPD.90.upper
#> C ALEU (Trend)        2.16         0.59         4.41
#> C GULF (Trend)        6.52         5.50         8.06
#> E ALEU (Trend)        3.55         2.42         5.16
#> E GULF (Trend)        7.36         4.49        10.52
#> W ALEU (Trend)       -2.68        -6.15         0.61
#> W GULF (Trend)        5.42         3.73         6.72
```
