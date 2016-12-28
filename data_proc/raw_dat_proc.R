library(reshape2)

### WDPS Nonpups
dat <- read.csv("data_proc/wdps_ssl_nonpups.csv")
nsites <- nrow(dat)
dat <- melt(dat, id.vars=1:10, variable.name="year")
dat$year <- as.numeric(sapply(strsplit(as.character(dat$year), "X"), function(x){x[[2]]}))
dat <- dat[!is.na(dat$value),]
dat <- dat[order(dat$Site, dat$year),]
names(dat)[ncol(dat)] <- "count"
wdpsNonpups <- dat
save(wdpsNonpups, file="data/wdpsNonpups.rda")

### Photo experiment
photoCorrection <-  read.csv("data_proc/photoCorrection.csv")
save(photoCorrection, file="data/photoCorrection.rda")

### WDPS pups
dat <- read.csv("data_proc/wdpsPups.csv")
nsites <- nrow(dat)
dat <- melt(dat, id.vars=c("Site","Region","RCA"), variable.name="year")
dat$year <- as.numeric(sapply(strsplit(as.character(dat$year), "X"), function(x){x[[2]]}))
dat <- dat[!is.na(dat$value),]
dat <- dat[order(dat$Site, dat$year),]
names(dat)[ncol(dat)] <- "count"
wdpsPups <- dat
save(wdpsPups, file="data/wdpsPups.rda")

dat <- read.csv("data_proc/edps_pups.csv")
edpsPups <- dat
save(edpsPups, file="data/edpsPups.rda")

dat <- read.csv("data_proc/edps_nonpups.csv")
edpsNonpups <- dat
save(edpsNonpups, file="data/edpsNonpups.rda")
