data(wdpsNonpups)
data(photoCorrection)

### Subset data to post-1989
wdpsNonpups <- wdpsNonpups[wdpsNonpups$year>=1990 , ]
wdpsNonpups <- droplevels(wdpsNonpups)

### Remove sites that had only 1 positive count
nz.counts <- aggregate(wdpsNonpups$count[wdpsNonpups$count>0], by=list(wdpsNonpups$site[wdpsNonpups$count>0]), FUN=length)
wdpsNonpups <- wdpsNonpups[wdpsNonpups$site%in%as.character(nz.counts[nz.counts[,2]>1,1]),]
wdpsNonpups <- droplevels(wdpsNonpups)

### Add photo method covariate to data (oblique photos prior to 2004 surveys = 1)
wdpsNonpups$obl <- as.numeric(wdpsNonpups$year<2004)


### Create prediction and zero-inflation (ZI) models for each site 
# Count number of surveys
num.surv <- aggregate(wdpsNonpups$count, by=list(wdpsNonpups$site), FUN=length)
# Count number of surveys with non-zero counts
nz.counts <- aggregate(wdpsNonpups$count[wdpsNonpups$count>0], by=list(wdpsNonpups$site[wdpsNonpups$count>0]), FUN=length)
# Set up model data.frame for 'mcmc.aggregate(...)' function
wdpsModels <- data.frame(site=num.surv[,1])
# Trend models:
#   0-5 nonzero counts = constant trend signal
#   6-10 nonzero counts = linear trend signal
#   >10 nonzero counts = RW2 (i.e., spline) trend signal
wdpsModels$trend <- cut(nz.counts[,2], c(0,5,10,30), labels=c("const","lin","RW"))
# Zero inflation models:
#   0-5 surveys = constant inflation effect
#   >5 surveys = linear inflation effect
#   All surveys have nonzero counts = no ZI model
wdpsModels$zero.infl <- cut(nz.counts[,2], c(0,5,30), labels=c("const","lin"))
levels(wdpsModels$zero.infl) <- c(levels(wdpsModels$zero.infl), "none")
wdpsModels$zero.infl[nz.counts[,2]==num.surv[,2]] <- "none"

head(wdpsModels)


### Create prior distribution list for MCMC site updating 
### Assumes site trends unlikely to be greater than 20% or less that about -17%
### or, more exactly Pr{exp(-ln(1+0.2))-1 < trend < 0.2} = 0.95
# Create informative gamma prior
x = log(photoCorrection$X2000OBL/photoCorrection$X2000VERT)
gamma.0 = mean(x)
gamma.se = sd(x)/sqrt(nrow(photoCorrection))

prior.list=defaultPriorList(trend.limit=0.2, model.data=wdpsModels,
                            gamma.mean=gamma.0, gamma.prec=1/gamma.se^2)

### Create upper bounds for predictive counts (= 3 x max(count) = 1.2^6)
upper <- aggregate(wdpsNonpups$count, list(wdpsNonpups$site), function(x){3*max(x)})
colnames(upper) <- c("site", "upper")

### Perform site augmentation and obtain posterior predictive distribution
### Note: Only a small number of MCMC iterations are shown here. For a more robust 
### analysis change to burn=1000 and iter=5000.
set.seed(123) 
fit <- mcmc.aggregate(start=1990, end=2012, data=wdpsNonpups, obs.formula=~obl-1, model.data=wdpsModels, 
                      rw.order=list(eta=2), aggregation="Region",
                      abund.name="count", time.name="year", site.name="site", 
                      burn=10, iter=50, thin=5, prior.list=prior.list, upper=upper, 
                      keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)


### Look at the results
fitdat <- fit$aggregation.pred.summary
# Compute trends for just 2000-2012
trend2000 <- updateTrend(fit, 2000, 2012, "pred")
# Change to % growth form
growth2000 <- mcmc(100*(exp(trend2000[,7:12])-1))
summary(growth2000)
# Obtain posterior median % growth and 90% credible interval
print(
  data.frame(
    post.median=round(apply(growth2000, 2, median),2),
    HPD.90=round(HPDinterval(growth2000, 0.90),2)
  )
)

# Add fitted 2000-2012 trend to aggregation summary
b <- apply(trend2000, 2, median)
X <- model.matrix(~(Region-1) + (Region-1):(year), data=fitdat)
fitdat$trend2000 <- apply(
  apply(as.matrix(trend2000), 1, FUN=function(b,Mat){as.vector(exp(Mat%*%b))}, Mat=X), 
  1, median
  )
fitdat$trend2000[fitdat$year<2000] <- NA

# Make a plot of the results (requires ggplot2 package)
library(ggplot2)
surv.yrs <- unique(fit$original.data$year)
ag.sum.data <- fit$aggregation.pred.summary
rel.dat <- fit$aggregation.rel.summary
rel.dat <- rel.dat[rel.dat$year %in% surv.yrs,]
colnames(rel.dat)[3:5] <- paste(colnames(rel.dat)[3:5], "REL", sep="")
ag.sum.data <- merge(ag.sum.data, rel.dat, all=TRUE)
ag.nm <- "Region"
ag.sum.data[,ag.nm] <- factor(ag.sum.data[,ag.nm])
b <- apply(trend2000, 2, median)
X <- model.matrix(~(ag.sum.data[,ag.nm]-1) + (ag.sum.data[,ag.nm]-1):(year), data=ag.sum.data)
ag.sum.data$trend2000 <- apply(apply(as.matrix(trend2000), 1, FUN=function(b,Mat){as.vector(exp(Mat%*%b))}, Mat=X), 1, median)
ag.sum.data$trend2000[ag.sum.data$year<2000] <- NA
ag.sum.data$Region = factor(ag.sum.data$Region, 
                            levels=c("W ALEU", "C ALEU", "E ALEU", "W GULF", "C GULF", "E GULF"))
wdpsfig <- ggplot(ag.sum.data, aes(x=year, y=post.median.abund)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=low90.hpd, ymax=hi90.hpd), alpha=0.4, na.rm=TRUE) + 
  geom_line(aes(y=trend2000), color="blue", lwd=1.5, na.rm=TRUE) + 
  geom_pointrange(aes(y=post.median.abundREL, ymin=low90.hpdREL, ymax=hi90.hpdREL), na.rm=TRUE) + 
  facet_wrap(~Region, ncol=2) +
  xlab("\nYear") + ylab("Aggregated count\n")

# Examine a couple of specific sites
site.pred <- fit$mcmc.sample$pred.site.abund
yr.site <- expand.grid(c(1990:2012), levels(wdpsNonpups$site))
colnames(yr.site) <- c("year","site")
yr.site <- merge(yr.site, wdpsModels, by="site")
# GLACIER in the E GULF
glacier.pred <- mcmc(site.pred[,yr.site$site=="GLACIER"])
glacier.dat <- fit$original.data[fit$original.data$site=="GLACIER",]
glacier.plot <- ggplot() + 
  geom_line(aes(x=c(1990:2012), y=apply(glacier.pred, 2, median))) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacier.pred)), alpha=0.25) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacier.pred, 0.5)), alpha=0.25) +
  geom_point(aes(y=count, x=year), data=glacier.dat, size=3) +
  xlab("Year") + ylab("Abundance")

# Zero inflation process for GLACIER
glacierZI <- mcmc(fit$mcmc.sample$prob.zero.infl[,yr.site$site[yr.site$zero.infl!="none"]=="GLACIER"])
glacier.zi <- ggplot() + 
  geom_line(aes(x=c(1990:2012), y=apply(glacierZI, 2, median))) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacierZI)), alpha=0.25) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacierZI, 0.5)), alpha=0.25) +
  geom_point(aes(y=1.0*c(glacier.dat$count>0), x=year), data=glacier.dat, size=3) +
  xlab("Year") + ylab("Probability survey count > 0")

# MARMOT in the C ALEU
marmot.pred <- mcmc(site.pred[,yr.site$site=="MARMOT"])
marmot.dat <- fit$original.data[fit$original.data$site=="MARMOT",]
marmot.plot <- ggplot() + 
  geom_line(aes(x=c(1990:2012), y=apply(marmot.pred, 2, median))) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred)), alpha=0.25) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred, 0.5)), alpha=0.25) +
  geom_point(aes(y=count, x=year), data=marmot.dat, size=3) +
  xlab("Year") + ylab("Abundance")

# Print the figs-- uncomment 'ggsave' lines to create PDFs
suppressWarnings(print(wdpsfig)) # ggplot2 printing throws warnings for the NAs in the plots
print(glacier.plot)
print(glacier.zi)
print(marmot.plot)
# ggsave("wdpstrends.pdf", wdpsfig)
# ggsave("glacierPred.pdf", glacier.plot)
# ggsave("glacierZI.pdf", glacier.zi)
# ggsave("marmotPred.pdf", marmot.plot)

# Save the results-- uncomment to save
## save(fit, file="wdpsNonpupsDemoResults.RData", compress=TRUE)
