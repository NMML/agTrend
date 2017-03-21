<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/NMML/crawl.svg?branch=devel)](https://travis-ci.org/NMML/crawl)

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
#> This is mgcv 1.8-17. For overview type 'help("mgcv-package")'.
#> Loading required package: truncnorm
#> 
#> 
#>  agTrend 0.17.7 (2017-03-20) 
#>  A demo is available at
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
wdpsNonpups %>% left_join(nz_counts) %>% filter(num_counts>1) %>% select(-num_counts) -> wdpsNonpups
#> Joining, by = "SITE"
```

Now we'll add a photo method covariate to data (oblique photos prior to 2004 surveys = 1)

``` r
wdpsNonpups %>% mutate(obl = as.integer(YEAR<2004)) -> wdpsNonpups
```

The next step is to create prediction and availability models for each site based on the number of surveys. The trend models will be

-   0-5 nonzero counts = constant trend signal
-   6-10 nonzero counts = linear trend signal
-   10 nonzero counts = RW2 (i.e., spline) trend signal

The zero inflation (avail) models are \* 0-5 surveys = constant inflation effect \* &gt;5 surveys = linear inflation effect \* All surveys have nonzero counts = no availability model

``` r
wdpsModels = wdpsNonpups %>% group_by(SITE) %>% 
  summarize(
    num_surv=n(),
    nz_count=sum(COUNT>0)) %>% 
  ungroup() %>%  
  mutate(
    trend = cut(num_surv, c(0,5,10,30), labels=c("const","lin","RW")),
    avail = cut(num_surv, c(0,5,30), labels=c("const","lin"))
  ) %>% 
  mutate(avail = ifelse(num_surv==nz_count, "none", avail))

head(wdpsModels)
#> # A tibble: 6 × 5
#>                 SITE num_surv nz_count  trend avail
#>                <chr>    <int>    <int> <fctr> <chr>
#> 1               ADAK       16       16     RW  none
#> 2 ADAK/ARGONNE POINT       13        8     RW     2
#> 3  ADAK/CRONE ISLAND       11        9     RW     2
#> 4             ADUGAK       18       18     RW  none
#> 5 AFOGNAK/TONKI CAPE       16        8     RW     2
#> 6             AGATTU       16       16     RW  none
```

### Create prior distribution list for MCMC site updating

### Assumes site trends unlikely to be greater than 20% or less that about -17%

### or, more exactly Pr{exp(-ln(1+0.2))-1 &lt; trend &lt; 0.2} = 0.95

Create informative gamma prior
==============================

x = log(photoCorrection*X*2000*O**B**L*/*p**h**o**t**o**C**o**r**r**e**c**t**i**o**n*X2000VERT) gamma.0 = mean(x) gamma.se = sd(x)/sqrt(nrow(photoCorrection))

prior.list=defaultPriorList(trend.limit=0.2, model.data=wdpsModels, gamma.mean=gamma.0, gamma.prec=1/gamma.se^2)

### Create upper bounds for predictive counts (= 3 x max(count) = 1.2^6)

upper &lt;- aggregate(wdpsNonpups*c**o**u**n**t*, *l**i**s**t*(*w**d**p**s**N**o**n**p**u**p**s*site), function(x){3\*max(x)}) colnames(upper) &lt;- c("site", "upper")

### Perform site augmentation and obtain posterior predictive distribution

### Note: Only a small number of MCMC iterations are shown here. For a more robust

### analysis change to burn=1000 and iter=5000.

set.seed(123) fit &lt;- mcmc.aggregate(start=1990, end=2012, data=wdpsNonpups, obs.formula=~obl-1, model.data=wdpsModels, rw.order=list(omega=2), aggregation="Region", abund.name="count", time.name="year", site.name="site", burn=50, iter=100, thin=5, prior.list=prior.list, upper=upper, keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)

### Look at the results

fitdat &lt;- fit$aggregation.pred.summary \# Compute trends for just 2000-2012 trend2000 &lt;- updateTrend(fit, 2000, 2012, "pred") \# Change to % growth form growth2000 &lt;- mcmc(100\*(exp(trend2000\[,7:12\])-1)) summary(growth2000) \# Obtain posterior median % growth and 90% credible interval print( data.frame( post.median=round(apply(growth2000, 2, median),2), HPD.90=round(HPDinterval(growth2000, 0.90),2) ) )

Add fitted 2000-2012 trend to aggregation summary
=================================================

b &lt;- apply(trend2000, 2, median) X &lt;- model.matrix(~(Region-1) + (Region-1):(year), data=fitdat) fitdat$trend2000 &lt;- apply(  apply(as.matrix(trend2000), 1, FUN=function(b,Mat){as.vector(exp(Mat%\*%b))}, Mat=X),  1, median ) fitdat$trend2000\[fitdat$year&lt;2000\] &lt;- NA

Make a plot of the results (requires ggplot2 package)
=====================================================

library(ggplot2)

envCol = "\#2b83ba" lnCol = "\#d7191c"

surv.yrs &lt;- unique(fit*o**r**i**g**i**n**a**l*.*d**a**t**a*year) ag.sum.data &lt;- fit*a**g**g**r**e**g**a**t**i**o**n*.*p**r**e**d*.*s**u**m**m**a**r**y**r**e**a**l*.*d**a**t* &lt; −*f**i**t*aggregation.real.summary real.dat &lt;- real.dat\[real.dat$year %in% surv.yrs,\] colnames(real.dat)\[3:5\] &lt;- paste(colnames(real.dat)\[3:5\], "REAL", sep="") ag.sum.data &lt;- merge(ag.sum.data, real.dat, all=TRUE) ag.nm &lt;- "Region" ag.sum.data\[,ag.nm\] &lt;- factor(ag.sum.data\[,ag.nm\]) b &lt;- apply(trend2000, 2, median) X &lt;- model.matrix(~(ag.sum.data\[,ag.nm\]-1) + (ag.sum.data\[,ag.nm\]-1):(year), data=ag.sum.data) ag.sum.data$trend2000 &lt;- apply(apply(as.matrix(trend2000), 1, FUN=function(b,Mat){as.vector(exp(Mat%\*%b))}, Mat=X), 1, median) ag.sum.data$trend2000\[ag.sum.data$year&lt;2000\] &lt;- NA ag.sum.data*R**e**g**i**o**n* = *f**a**c**t**o**r*(*a**g*.*s**u**m*.*d**a**t**a*Region, levels=c("W ALEU", "C ALEU", "E ALEU", "W GULF", "C GULF", "E GULF")) fig1 &lt;- ggplot(ag.sum.data, aes(x=year, y=post.median.abund)) + facet\_wrap(~Region, ncol=2) + geom\_line() + geom\_ribbon(aes(ymin=low.hpd, ymax=hi.hpd), alpha=0.4, fill=envCol) + geom\_line(aes(y=trend2000), color=lnCol, lwd=1.5, data=ag.sum.data\[!is.na(ag.sum.data$trend2000),\]) + geom\_pointrange(aes(y=post.median.abundREAL, ymin=low.hpdREAL, ymax=hi.hpdREAL), data=ag.sum.data\[!is.na(ag.sum.data$post.median.abundREAL),\]) + xlab("") + ylab("Aggregated count") + theme\_bw() + theme(panel.grid=element\_blank(), text=element\_text(size=14))

print(fig1) \# ggsave(fig1, file="figure/figure1.pdf", width=6.5, height=8)

suppressMessages(library(gridExtra))

site.pred &lt;- fit*m**c**m**c*.*s**a**m**p**l**e*pred.site.abund yr.site &lt;- expand.grid(c(1990:2012), levels(wdpsNonpups$site)) colnames(yr.site) &lt;- c("year","site") yr.site &lt;- merge(yr.site, wdpsModels, by="site")

GLACIER in the E GULF
=====================

glacier.pred &lt;- mcmc(site.pred\[,yr.site$site=="GLACIER"\]) glacier.dat &lt;- fit*o**r**i**g**i**n**a**l*.*d**a**t**a*\[*f**i**t*original.data$site=="GLACIER",\] glacier.plot &lt;- ggplot() + geom\_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacier.pred)), alpha=0.4, fill=envCol) + geom\_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacier.pred, 0.5)), alpha=0.4, fill=envCol) + geom\_line(aes(x=c(1990:2012), y=apply(glacier.pred, 2, median))) + geom\_point(aes(y=count, x=year), data=glacier.dat, size=3) + xlab("Year") + ylab("Survey count") + ggtitle("(a) Counts at Glacier") + theme\_bw() + theme(panel.grid=element\_blank(), text=element\_text(size=12), plot.title=element\_text(size=12))

Zero inflation process for GLACIER
==================================

glacierAV &lt;- mcmc(fit*m**c**m**c*.*s**a**m**p**l**e*prob.avail\[,yr.site*s**i**t**e*\[*y**r*.*s**i**t**e*avail!="none"\]=="GLACIER"\]) glacier.av &lt;- ggplot() + geom\_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacierAV)), alpha=0.4, fill=envCol) + geom\_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacierAV, 0.5)), alpha=0.4, fill=envCol) + geom\_line(aes(x=c(1990:2012), y=apply(glacierAV, 2, median))) + geom\_point(aes(y=1.0\*c(glacier.dat$count&gt;0), x=year), data=glacier.dat, size=3) + xlab("Year") + ylab("Probability survey count &gt; 0") + ggtitle("(b) Availability at Glacier") + theme\_bw() + theme(panel.grid=element\_blank(), text=element\_text(size=12), plot.title=element\_text(size=12))

MARMOT in the C ALEU
====================

marmot.pred &lt;- mcmc(site.pred\[,yr.site$site=="MARMOT"\]) marmot.dat &lt;- fit*o**r**i**g**i**n**a**l*.*d**a**t**a*\[*f**i**t*original.data$site=="MARMOT",\] marmot.plot &lt;- ggplot() + geom\_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred)), alpha=0.4, fill=envCol) + geom\_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred, 0.5)), alpha=0.4, fill=envCol) + geom\_line(aes(x=c(1990:2012), y=apply(marmot.pred, 2, median))) + geom\_point(aes(y=count, x=year), data=marmot.dat, size=3) + xlab("Year") + ylab("Survey count") + ggtitle("(c) Counts at Marmot") + theme\_bw() + theme(panel.grid=element\_blank(), text=element\_text(size=12), plot.title=element\_text(size=12))

fig2 &lt;- arrangeGrob(glacier.plot, glacier.av, marmot.plot, ncol=2) print(fig2) \# ggsave(fig2, file="figure/figure2.pdf", width=6.5, height=6.5)

Save the results-- uncomment to save
====================================

save(list=ls(), file="wdpsNonpupsDemoResults.RData", compress=TRUE)
===================================================================
