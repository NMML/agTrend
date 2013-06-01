library(agTrend)
data(wdpsNonpups)
data(photoCorrection)

### Subset data to post-1989
wdpsNonpups <- wdpsNonpups[wdpsNonpups$year>=1990, ]
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
upper <- aggregate(wdpsNonpups$count, list(wdpsNonpups$site), function(x){(1.05^50)*max(x)})
colnames(upper) <- c("site", "upper")

### Perform site augmentation and obtain posterior predictive distribution
set.seed(123) 
fit <- mcmc.aggregate(start=1990, end=2062, data=wdpsNonpups, obs.formula=~obl-1, model.data=wdpsModels, aggregation="Region",
                      abund.name="count", time.name="year", site.name="site", forecast=TRUE,
                      burn=1000, iter=15000, thin=5, prior.list=prior.list, upper=upper, 
                      keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)


### Look at the results
ag.sum.data <- fit$aggregation.pred.summary
surv.yrs <- unique(fit$original.data$year)

rel.dat <- fit$aggregation.rel.summary
rel.dat <- rel.dat[rel.dat$year %in% surv.yrs,]
colnames(rel.dat)[3:5] <- paste(colnames(rel.dat)[3:5], "REL", sep="")
ag.sum.data <- merge(ag.sum.data, rel.dat, all=TRUE)
ag.sum.data[,"Region"] <- factor(ag.sum.data[,"Region"])


# Make a plot of the results (requires ggplot2 package)
library(ggplot2)
ggfig <- ggplot(ag.sum.data, aes(x=year, y=post.median.abund)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=low90.hpd, ymax=hi90.hpd), alpha=0.4) + 
  geom_pointrange(aes(y=post.median.abundREL, ymin=low90.hpdREL, ymax=hi90.hpdREL)) + 
  facet_wrap(~ Region, ncol=3) +
  xlab("\nYear") + ylab("Forecast count\n") + ggtitle("wDPS Region\n")
ggsave("wdpstrendsForecast.pdf", ggfig)

pe <- data.frame(Region=c("W ALEU", "C ALEU", "E ALEU", "W GULF", "C GULF", "E GULF", "Total"),
                 rookeries = c(4, 12, 7, 5, 6, 3, 37),
                 pe.abund=round(4743*c(4, 12, 7, 5, 6, 3, 37)/sum(c(4, 12, 7, 5, 6, 3)),0),
                 pe.cnt=round(0.5*4743*c(4, 12, 7, 5, 6, 3, 37)/sum(c(4, 12, 7, 5, 6, 3)),0),
                 Extinction.Prob.20=0,
                 Extinction.Prob.50=0
                 )
extinct <- function(x,val,year=73){1.0*any(x[1:year] < val)}

region.mcmc <- fit$mcmc.sample$aggregated.pred.abund

for(i in 1:6){
  idx <- grep(pe$Region[i], colnames(region.mcmc))
  Xmat <- region.mcmc[,idx]
  pe$Extinction.Prob.20[i] <- mean(apply(Xmat, 1, extinct, val=pe$pe.cnt[i], year=43))
  pe$Extinction.Prob.50[i] <- mean(apply(Xmat, 1, extinct, val=pe$pe.cnt[i], year=73))
}


# look at WDPS as a whole
new.agg.data <- data.frame(site=levels(fit$original.data$site))
new.agg.data$total <- "Total"
new.agg.list.pred <- newAggregation(fit, new.agg.data, "pred")
total.mcmc <- new.agg.list.pred[[1]]$mcmc.sample
ag.sum.data <- new.agg.list.pred[[1]]$aggregation.pred.summary
ggfig <- ggplot(ag.sum.data, aes(x=year, y=post.median.abund)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=low90.hpd, ymax=hi90.hpd), alpha=0.4) + 
  xlab("\nYear") + ylab("Forecast count\n") + ggtitle("wDPS Total\n")
ggsave("wdpstrendsForecast_total.pdf", ggfig)

Xmat <- total.mcmc$aggregated.pred.abund
pe$Extinction.Prob.20[7] <- mean(apply(Xmat, 1, extinct, val=pe$pe.cnt[7], year=43))
pe$Extinction.Prob.50[7] <- mean(apply(Xmat, 1, extinct, val=pe$pe.cnt[7], year=73))

# Save the results-- uncomment to save
write.csv(pe, file="wdpsPEresults.csv", row.names=FALSE)
save(fit, pe, new.agg.list.pred, file="wdpsNonpupsForecastDemoResults.rda", compress=TRUE)


