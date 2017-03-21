<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/NMML/crawl.svg?branch=devel)](https://travis-ci.org/NMML/crawl)

Fit regional trends to site-specific abundence data
---------------------------------------------------

This package fits a log-linear trend models to regions aggregated over sites. The sites may contain missing surveys that are not temporally aligned with the missing data at other sites, making direct aggregation impossible. The functions within the package model the indivdual sites with a semi-parametric (possibly, zero-inflated) model to interpolate missing data from which regional aggregations can be made. By using Markov Chain Monte Carlo, on can sample from the posterior predictive distribution of the regional aggregations Then calculate the log-linear trend over the time period of interest as a derived parameter. Using the posterior predictive distribution allows incorporation of both parameter uncertainty as well as uncertainty due to sampling the local abundance processes.

### Disclaimer

*This software package is developed and maintained by scientists at the NOAA Fisheries Alaska Fisheries Science Center and should be considered a fundamental research communication. The reccomendations and conclusions presented here are those of the authors and this software should not be construed as official communication by NMFS, NOAA, or the U.S. Dept. of Commerce. In addition, reference to trade names does not imply endorsement by the National Marine Fisheries Service, NOAA. While the best efforts have been made to insure the highest quality, tools such as this are under constant development and are subject to change.*

### Example