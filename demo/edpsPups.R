data(edpsPups)

# Multiply pup counts by 4.5 to get estimate of total abundance
edpsPups$abund <- edpsPups$count*4.5

# Fit the site models and estimate aggregated trend summary
fit <- mcmc.aggregate(start=1979, end=2010, data=edpsPups, rw.order=2,
                      abund.name="count", time.name="year", site.name="region", 
                      burn=5000, iter=10000, thin=1, prior.list=NULL, ln.adj=0, keep.site.param=TRUE, 
                      keep.site.abund=TRUE)

fitdat <- fit$aggregation.summary

# Make a plot of the results (requires ggplot2 package)
library(ggplot2)
ggplot(fitdat, aes(x=year, y=4.5*count, color=region)) +
  geom_point(size=3) +
  geom_line(aes(y=4.5*post.mean.abund)) +
  geom_line(aes(y=fit, x=year), color="black", lwd=2, data=data.frame(year=1979:2010, fit=4.5*exp(fit$trend.summary$statistics[1,1]+fit$trend.summary$statistics[2,1]*c(1979:2010)))) +
  geom_ribbon(aes(ymin=4.5*low90.hpd, ymax=4.5*hi90.hpd, fill=region), alpha=0.25) +
  xlab("Year") + ylab("eDPS estimated abundance")

# Examine individual site trends
summary(fit$mcmc.sample$site.param$beta)

# Estimate posterior predictive mode and HPD credible interval growth rate (percentage form) for EDPS stock
dd <- density(exp(fit$mcmc.sample$pred.trend[,2])-1)
dd$x[dd$y==max(dd$y)]
HPDinterval(exp(fit$mcmc.sample$pred.trend[,2])-1, 0.9)

