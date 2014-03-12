data(edpsPups)

# Multiply pup counts by 4.5 to get estimate of total abundance
edpsPups$abund <- edpsPups$count*4.5

# Define models to update within each region
edpsModels <- data.frame(region=levels(edpsPups$region), trend=c("RW","lin","lin","RW"))

# Define priors for regional trend parameters
prior.list=defaultPriorList(trend.limit=0.2, model.data=edpsModels)

### Create upper bounds for predictive counts (= 3 x max(count) = 1.2^6)
upper <- aggregate(edpsPups$abund, list(edpsPups$region), function(x){3*max(x)})
colnames(upper) <- c("region", "upper")

### Perform site augmentation and obtain posterior predictive distribution
### Note: Only a small number of MCMC iterations are shown here. For a more robust 
### analysis change to burn=1000 and iter=5000.
set.seed(123)
fit <- mcmc.aggregate(start=1979, end=2010, data=edpsPups, model.data=edpsModels, rw.order=list(eta=2),
                      abund.name="abund", time.name="year", site.name="region", 
                      burn=10, iter=50, thin=5, prior.list=prior.list, upper=upper, 
                      keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)

# Extract MCMC summaries
fitdat <- fit$aggregation.pred.summary
colnames(fitdat)[2] <- "region"
fitdat <- rbind(fitdat, fit$site.summary)
fitdat <- merge(fitdat, fit$original.data, all.x=TRUE)
fitdat <- fitdat[order(fitdat$region, fitdat$year),]

# Make a plot of the results (requires ggplot2 package)
library(ggplot2)
plt <- ggplot(aes(x=year, y=abund, color=region), data=fitdat) + 
  geom_point(size=3) +
  geom_line(aes(y=post.median.abund)) +
  geom_ribbon(aes(ymin=low90.hpd, ymax=hi90.hpd, fill=region), alpha=0.25) +
  xlab("\nYear") + ylab("Estimated SSL abundance (4.5 x pup count)\n") + ggtitle("eDPS 30 Year Trend\n")
suppressWarnings(print(plt))
# Uncomment to save figure:
# ggsave("edpsTotalTrend.pdf", plt)

# Examine individual site trends
summary(fit$mcmc.sample$site.param$beta)

# Estimate posterior predictive mode and HPD credible interval growth rate (percentage form) for EDPS stock
trend.mcmc = exp(fit$mcmc.sample$pred.trend[,2])-1
dd <- density(trend.mcmc)
trend.mode = dd$x[dd$y==max(dd$y)]
trend.HPD = HPDinterval(trend.mcmc, 0.95)
cat(paste("\nTrend estimate:", 100*round(trend.mode,3), 
          "(", 100*round(trend.HPD[,1],3),",",100*round(trend.HPD[,2],3),")","\n"))

# Save the results-- uncomment to save
# save(fit, file="edpsPupsDemoResults.rda", compress=TRUE)

