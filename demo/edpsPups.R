data(edpsPups)

# Multiply pup counts by 4.5 to get estimate of total abundance
edpsPups$abund <- edpsPups$count*4.5

# Define models to update within each region
edpsModels <- data.frame(region=levels(edpsPups$region), trend=c("RW2","lin","lin","RW2"))

# Define priors for regional trend parameters
library(Matrix)
n.sites <- length(unique(edpsPups$region))
prior.list <- list( 
  # Set Q.beta such that average trend not likely to exceed +/- 20%/year from 1979-2010
  beta=list(
    beta.0=rep(0,sum(edpsModels$trend=="const")+(sum(edpsModels$trend!="const")*2)), 
    Q.beta=Matrix(diag(
      unlist(lapply(edpsModels$trend, 
                    function(x){
                      if(x!="const") return(c(0,200))
                      else return(c(0))
                    }
      ))
    ))
  )
)

### Create upper bounds for predictive counts (= 3 x max(count) = 1.2^6)
upper <- aggregate(edpsPups$abund, list(edpsPups$region), function(x){3*max(x)})
colnames(upper) <- c("region", "upper")


# Fit the site models and estimate aggregated trend summary
set.seed(123)
fit <- mcmc.aggregate(start=1979, end=2010, data=edpsPups, model.data=edpsModels,
                      abund.name="abund", time.name="year", site.name="region", 
                      burn=1000, iter=5000, thin=5, prior.list=prior.list, upper=upper, 
                      keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)

# Extract MCMC summaries
fitdat <- fit$aggregation.pred.summary
colnames(fitdat)[2] <- "region"
fitdat <- rbind(fitdat, fit$site.summary)
fitdat <- merge(fitdat, fit$original.data, all.x=TRUE)
fitdat <- fitdat[order(fitdat$region, fitdat$year),]

# Make a plot of the results (requires ggplot2 package)
library(ggplot2)
plt <- ggplot(aes(x=year, y=abund, color=region), data=fitdat) + #[fitdat$region%in%c("Total","CA","OR"),]) +
  geom_point(size=3) +
  geom_line(aes(y=post.median.abund)) +
  geom_ribbon(aes(ymin=low90.hpd, ymax=hi90.hpd, fill=region), alpha=0.25) +
  xlab("\nYear") + ylab("Estimated SSL abundance (4.5 x pup count)\n") + ggtitle("eDPS 30 Year Trend\n")
ggsave("edpsTotalTrend.pdf", plt)

# Examine individual site trends
summary(fit$mcmc.sample$site.param$beta)

# Estimate posterior predictive mode and HPD credible interval growth rate (percentage form) for EDPS stock
dd <- density(exp(fit$mcmc.sample$pred.trend[,2])-1)
dd$x[dd$y==max(dd$y)]
HPDinterval(exp(fit$mcmc.sample$pred.trend[,2])-1, 0.95)

# Save the results-- uncomment to save
save(fit, file="edpsPupsDemoResults.rda", compress=TRUE)

