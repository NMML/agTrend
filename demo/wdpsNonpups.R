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


### Create prediction and availability models for each site 
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
#   All surveys have nonzero counts = no availability model
wdpsModels$avail <- cut(nz.counts[,2], c(0,5,30), labels=c("const","lin"))
levels(wdpsModels$avail) <- c(levels(wdpsModels$avail), "none")
wdpsModels$avail[nz.counts[,2]==num.surv[,2]] <- "none"

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
                      rw.order=list(omega=2), aggregation="Region",
                      abund.name="count", time.name="year", site.name="site", 
                      burn=50, iter=100, thin=5, prior.list=prior.list, upper=upper, 
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

envCol = "#2b83ba"
lnCol = "#d7191c"



surv.yrs <- unique(fit$original.data$year)
ag.sum.data <- fit$aggregation.pred.summary
real.dat <- fit$aggregation.real.summary
real.dat <- real.dat[real.dat$year %in% surv.yrs,]
colnames(real.dat)[3:5] <- paste(colnames(real.dat)[3:5], "REAL", sep="")
ag.sum.data <- merge(ag.sum.data, real.dat, all=TRUE)
ag.nm <- "Region"
ag.sum.data[,ag.nm] <- factor(ag.sum.data[,ag.nm])
b <- apply(trend2000, 2, median)
X <- model.matrix(~(ag.sum.data[,ag.nm]-1) + (ag.sum.data[,ag.nm]-1):(year), data=ag.sum.data)
ag.sum.data$trend2000 <- apply(apply(as.matrix(trend2000), 1, FUN=function(b,Mat){as.vector(exp(Mat%*%b))}, Mat=X), 1, median)
ag.sum.data$trend2000[ag.sum.data$year<2000] <- NA
ag.sum.data$Region = factor(ag.sum.data$Region, 
                            levels=c("W ALEU", "C ALEU", "E ALEU", "W GULF", "C GULF", "E GULF"))
fig1 <- ggplot(ag.sum.data, aes(x=year, y=post.median.abund)) + 
  facet_wrap(~Region, ncol=2) +
  geom_line() + 
  geom_ribbon(aes(ymin=low.hpd, ymax=hi.hpd), alpha=0.4, fill=envCol) + 
  geom_line(aes(y=trend2000), color=lnCol, lwd=1.5, data=ag.sum.data[!is.na(ag.sum.data$trend2000),]) + 
  geom_pointrange(aes(y=post.median.abundREAL, ymin=low.hpdREAL, ymax=hi.hpdREAL), data=ag.sum.data[!is.na(ag.sum.data$post.median.abundREAL),]) + 
  xlab("\nYear") + ylab("Aggregated count\n") + 
  theme_bw() + theme(panel.grid=element_blank(), text=element_text(size=14)) 

print(fig1)
# ggsave(fig1, file="figure/figure1.pdf", width=6.5, height=8)

suppressMessages(library(gridExtra))

site.pred <- fit$mcmc.sample$pred.site.abund
yr.site <- expand.grid(c(1990:2012), levels(wdpsNonpups$site))
colnames(yr.site) <- c("year","site")
yr.site <- merge(yr.site, wdpsModels, by="site")

# GLACIER in the E GULF
glacier.pred <- mcmc(site.pred[,yr.site$site=="GLACIER"])
glacier.dat <- fit$original.data[fit$original.data$site=="GLACIER",]
glacier.plot <- ggplot() + 
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacier.pred)), alpha=0.4, fill=envCol) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacier.pred, 0.5)), alpha=0.4, fill=envCol) +
  geom_line(aes(x=c(1990:2012), y=apply(glacier.pred, 2, median))) +
  geom_point(aes(y=count, x=year), data=glacier.dat, size=3) +
  xlab("Year") + ylab("Survey count") + ggtitle("(a) Counts at Glacier\n") +
  theme_bw() + theme(panel.grid=element_blank(), text=element_text(size=12), plot.title=element_text(size=12))

# Zero inflation process for GLACIER
glacierAV <- mcmc(fit$mcmc.sample$prob.avail[,yr.site$site[yr.site$avail!="none"]=="GLACIER"])
glacier.av <- ggplot() + 
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacierAV)), alpha=0.4, fill=envCol) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacierAV, 0.5)), alpha=0.4, fill=envCol) +
  geom_line(aes(x=c(1990:2012), y=apply(glacierAV, 2, median))) +
  geom_point(aes(y=1.0*c(glacier.dat$count>0), x=year), data=glacier.dat, size=3) +
  xlab("Year") + ylab("Probability survey count > 0") + ggtitle("(b) Availability at Glacier\n") +
  theme_bw() + theme(panel.grid=element_blank(), text=element_text(size=12), plot.title=element_text(size=12))

# MARMOT in the C ALEU
marmot.pred <- mcmc(site.pred[,yr.site$site=="MARMOT"])
marmot.dat <- fit$original.data[fit$original.data$site=="MARMOT",]
marmot.plot <- ggplot() + 
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred)), alpha=0.4, fill=envCol) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred, 0.5)), alpha=0.4, fill=envCol) +
  geom_line(aes(x=c(1990:2012), y=apply(marmot.pred, 2, median))) +
  geom_point(aes(y=count, x=year), data=marmot.dat, size=3) +
  xlab("Year") + ylab("Survey count") + ggtitle("(c) Counts at Marmot\n") + 
  theme_bw() + theme(panel.grid=element_blank(), text=element_text(size=12), plot.title=element_text(size=12))

fig2 <- arrangeGrob(glacier.plot, glacier.av, marmot.plot, ncol=2)
print(fig2)
# ggsave(fig2, file="figure/figure2.pdf", width=6.5, height=6.5)


# Save the results-- uncomment to save
# save(list=ls(), file="wdpsNonpupsDemoResults.RData", compress=TRUE)
