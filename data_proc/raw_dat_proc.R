library(tidyverse)
library(stringr)
library(readxl)

### WDPS Nonpups
wdpsNonpups = read_excel("data_proc/ALLCOUNTS_v3.xlsx", "wdpsnp") %>% 
  select(.,which(colnames(.)!="")) %>% 
  mutate(SITE=toupper(SITE), REGION=toupper(REGION)) %>% 
  gather(YEAR, COUNT, -c(1:3)) %>% 
  mutate(
    YEAR = as.integer(YEAR), 
    COUNT=as.integer(COUNT),
    RCA=as.character(RCA)
  ) %>% 
  filter(!is.na(COUNT)) %>% arrange(SITE, YEAR)
save(wdpsNonpups, file="data/wdpsNonpups.rda")


### Photo experiment
photoCorrection = read.csv("data_proc/photoCorrection.csv") %>% 
  rename(SITE=site, REGION=Region, OBLIQUE=X2000OBL, VERTICAL=X2000VERT) 
save(photoCorrection, file="data/photoCorrection.rda")

### WDPS pups
wdpsPups = read_excel("data_proc/ALLCOUNTS_v3.xlsx", "wdpspup") %>% 
  select(.,which(colnames(.)!="")) %>% 
  mutate(SITE=toupper(SITE), REGION=toupper(REGION)) %>% 
  gather(YEAR, COUNT, -c(1:3)) %>% 
  mutate(
    YEAR = as.integer(YEAR), 
    COUNT=as.integer(COUNT),
    RCA=as.character(RCA)
  ) %>% 
  filter(!is.na(COUNT)) %>% arrange(SITE, YEAR)
save(wdpsPups, file="data/wdpsPups.rda")

edpsPups = read_excel("data_proc/ALLCOUNTS_v3.xlsx", "edpspup") %>% 
  select(.,which(colnames(.)!="")) %>% 
  mutate(SITE=toupper(SITE), REGION=toupper(REGION)) %>% 
  gather(YEAR, COUNT, -c(1:2)) %>% 
  filter(!is.na(COUNT)) %>% arrange(SITE, YEAR)
save(edpsPups, file="data/edpsPups.rda")

edpsNonpups = read_excel("data_proc/ALLCOUNTS_v3.xlsx", "edpsnp") %>% 
  select(.,which(colnames(.)!="")) %>% 
  mutate(SITE=toupper(SITE), REGION=toupper(REGION)) %>% 
  gather(YEAR, COUNT, -c(1:2)) %>% 
  filter(!is.na(COUNT)) %>% arrange(SITE, YEAR)
save(edpsNonpups, file="data/edpsNonpups.rda")
