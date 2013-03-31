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
n.sites <- length(unique(wdpsNonpups$site))
prior.list <- list(
  # Use sampling distribution for oblique/med. format pilot study
  gamma=list(gamma.0=mean(log(photoCorrection$X2000OBL) - log(photoCorrection$X2000VERT)),
             Q.gamma=(nrow(photoCorrection)-1)/var(log(photoCorrection$X2000OBL) - log(photoCorrection$X2000VERT))
             ), 
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
upper <- aggregate(wdpsNonpups$count, list(wdpsNonpups$site), function(x){3*max(x)})
colnames(upper) <- c("site", "upper")

### Perform site augmentation and obtain posterior predictive distribution
set.seed(123) 
fit <- mcmc.aggregate(start=1990, end=2012, data=wdpsNonpups, obs.formula=~obl-1, model.data=wdpsModels, aggregation="Region",
                      abund.name="count", time.name="year", site.name="site", 
                      burn=5000, iter=10000, thin=5, prior.list=prior.list, upper=upper, 
                      keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)


### Look at the results
fitdat <- fit$aggregation.pred.summary
# Compute trends for just 2000-2012
trend2000 <- updateTrend(fit, 2000, 2012, "pred")
# Change to % growth form
growth2000 <- mcmc(100*(exp(trend2000[,7:12])-1))
summary(growth2000)
# Obtain posterior median % growth and 90% credible interval
data.frame(
  post.median=round(apply(growth2000, 2, median),2),
  HPD.90=round(HPDinterval(growth2000, 0.90),2)
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
wdpsfig <- ggplot(fitdat, aes(x=year, y=post.median.abund, color=Region)) +
  geom_line(aes(y=post.median.abund)) +
  geom_ribbon(aes(ymin=low90.hpd, ymax=hi90.hpd, fill=Region), alpha=0.15) +
  geom_line(aes(y=trend2000,color=Region),lwd=3, data=fitdat[fitdat$year>=2000,]) + 
  xlab("Year") + ylab("WDPS estimated abundance")
ggsave("wdpstrends.pdf", wdpsfig)

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
ggsave("glacierPred.pdf", glacier.plot)

# Zero inflation process for GLACIER
glacierZI <- mcmc(fit$mcmc.sample$prob.zero.infl[,yr.site$site[yr.site$zero.infl!="none"]=="GLACIER"])
glacier.plot <- ggplot() + 
  geom_line(aes(x=c(1990:2012), y=apply(glacierZI, 2, median))) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacierZI)), alpha=0.25) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(glacierZI, 0.5)), alpha=0.25) +
  geom_point(aes(y=1.0*c(glacier.dat$count>0), x=year), data=glacier.dat, size=3) +
  xlab("Year") + ylab("Probability survey count > 0")
ggsave("glacierZI.pdf", glacier.plot)

# MARMOT in the C ALEU
marmot.pred <- mcmc(site.pred[,yr.site$site=="MARMOT"])
marmot.dat <- fit$original.data[fit$original.data$site=="MARMOT",]
marmot.plot <- ggplot() + 
  geom_line(aes(x=c(1990:2012), y=apply(marmot.pred, 2, median))) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred)), alpha=0.25) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=c(1990:2012)), data=data.frame(HPDinterval(marmot.pred, 0.5)), alpha=0.25) +
  geom_point(aes(y=count, x=year), data=marmot.dat, size=3) +
  xlab("Year") + ylab("Abundance")
ggsave("marmotPred.pdf", marmot.plot)

# Save the results-- uncomment to save
save(fit, file="wdpsNonpupsDemoResults.rda", compress=TRUE)
