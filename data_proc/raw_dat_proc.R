library(tidyverse)
library(stringr)

### WDPS Nonpups
wdpsNonpups = read.csv("data_proc/wdps_ssl_nonpups.csv") %>% 
  gather(year, count, X1978:X2016) %>% 
  mutate(year = as.numeric(str_sub(year, -4)), count=as.integer(count)) %>% 
  filter(!is.na(count)) %>% arrange(Site)
save(wdpsNonpups, file="data/wdpsNonpups.rda")


### Photo experiment
photoCorrection = read.csv("data_proc/photoCorrection.csv") %>% rename(Site=site)
save(photoCorrection, file="data/photoCorrection.rda")

### WDPS pups
wdpsPups <- read.csv("data_proc/wdpsPups.csv") %>% gather(year, count, -c(1:3)) %>% 
  mutate(year=str_sub(year,-4), count=as.integer(count)) %>% 
  filter(!is.na(count)) %>% arrange(Site)
save(wdpsPups, file="data/wdpsPups.rda")

edpsPups <- read.csv("data_proc/edps_pups.csv")
save(edpsPups, file="data/edpsPups.rda")

edpsNonpups <- read.csv("data_proc/edps_nonpups.csv")
save(edpsNonpups, file="data/edpsNonpups.rda")
