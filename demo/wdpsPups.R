data(wdpsPups)

### Subset data to post-1989
wdpsPups <- wdpsPups[wdpsPups$year>=1990 , ]
wdpsPups <- droplevels(wdpsPups)

### Remove sites that had only 1 positive count
nz.counts <- aggregate(wdpsPups$count[wdpsPups$count>0], by=list(wdpsPups$site[wdpsPups$count>0]), FUN=length)
wdpsPups <- wdpsPups[wdpsPups$site%in%as.character(nz.counts[nz.counts[,2]>1,1]),]
wdpsPups <- droplevels(wdpsPups)

### Create prediction and zero-inflation (ZI) models for each site 
# Count number of surveys
num.surv <- aggregate(wdpsPups$count, by=list(wdpsPups$site), FUN=length)
# Count number of surveys with non-zero counts
nz.counts <- aggregate(wdpsPups$count[wdpsPups$count>0], by=list(wdpsPups$site[wdpsPups$count>0]), FUN=length)
# Set up model data.frame for 'mcmc.aggregate(...)' function
wdpsModels <- data.frame(site=num.surv[,1])
# Trend models:
#   0-5 nonzero counts = constant trend signal
#   6-10 nonzero counts = linear trend signal
#   >10 nonzero counts = RW2 (i.e., spline) trend signal
wdpsModels$trend <- cut(nz.counts[,2], c(0,5,10,30), labels=c("const","lin","RW2"))
# Zero inflation models:
#   0-5 surveys = constant inflation effect
#   >5 surveys = linear inflation effect
#   All surveys have nonzero counts = no ZI model
wdpsModels$zero.infl <- cut(nz.counts[,2], c(0,5,30), labels=c("const","lin"))
levels(wdpsModels$zero.infl) <- c(levels(wdpsModels$zero.infl), "none")
wdpsModels$zero.infl[nz.counts[,2]==num.surv[,2]] <- "none"

head(wdpsModels)


### Create prior distribution list for MCMC site updating 
library(Matrix)
n.sites <- length(unique(wdpsPups$site))
prior.list <- list( 
  # Set Q.beta such that average trend not likely to exceed +/- 20%/year from 1990-2012
  beta=list(
    beta.0=rep(0,sum(wdpsModels$trend=="const")+(sum(wdpsModels$trend!="const")*2)), 
    Q.beta=Matrix(diag(
      unlist(lapply(wdpsModels$trend, 
        function(x){
          if(x!="const") return(c(0,200))
          else return(c(0))
        }
      ))
    ))
  )
)

### Create upper bounds for predictive counts (= 3 x max(count) = 1.2^6)
upper <- aggregate(wdpsPups$count, list(wdpsPups$site), function(x){3*max(x)})
colnames(upper) <- c("site", "upper")

### Perform site augmentation and obtain posterior predictive distribution
fit <- mcmc.aggregate(start=1990, end=2012, data=wdpsPups, model.data=wdpsModels, #aggregation="Region",
                      abund.name="count", time.name="year", site.name="site", 
                      burn=100, iter=1000, thin=5, prior.list=prior.list, upper=upper, 
                      keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)


### Look at the results
fitdat <- fit$aggregation.summary
# Compute trends for just 2000-2012
trend2000 <- updateTrend(2000, 2012, fit)
# Change to % growth form
growth2000 <- mcmc(100*(exp(trend2000[,2])-1))
summary(growth2000)
# Obtain posterior median % growth and 90% credible interval
data.frame(
  post.median=median(growth2000),
  HPD.90=round(HPDinterval(growth2000, 0.90),2)
  )

# Add fitted 2000-2012 trendd to aggregation summary
b <- apply(trend2000, 2, median)
X <- model.matrix(~year, data=fitdat)
fitdat$trend2000 <- exp(X%*%b)
fitdat$trend2000[fitdat$year<2000] <- NA

# Make a plot of the results (requires ggplot2 package)
library(ggplot2)
wdpsfig <- ggplot(fitdat, aes(x=year, y=post.median.abund)) +
  geom_line(aes(y=post.median.abund)) +
  geom_ribbon(aes(ymin=low90.hpd, ymax=hi90.hpd), alpha=0.15) +
  geom_line(aes(y=trend2000),lwd=3, data=fitdat[fitdat$year>=2000,], color="blue") + 
  xlab("Year") + ylab("WDPS estimated pup production")
ggsave("wdpstrends_pup.png", wdpsfig, width=6.5, height=6.5, units="in")


# MARMOT in the C ALEU
marmot.pred <- mcmc(site.pred[,yr.site$site=="MARMOT"])
marmot.dat <- fit$original.data[fit$original.data$site=="MARMOT",]
marmot.plot <- ggplot() + 
  geom_line(aes(x=c(1990:2012), y=apply(marmot.pred, 2, median))) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred)), alpha=0.25) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred, 0.5)), alpha=0.25) +
  geom_point(aes(y=count, x=year), data=marmot.dat, size=3) +
  xlab("Year") + ylab("Abundance")
ggsave("marmotPred_pup.png", marmot.plot, dpi=300)
