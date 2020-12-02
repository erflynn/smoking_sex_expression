# 01_load_cel_files.R
# 11/22/2020
# 
# Code for loading CEL files and running RMA
# set up for GPL570 + GPL96 airway epithelium data
#
# TODO:
#  - expand to larger dataset
#  - move data files to intermediate location

library(tidyverse)
library(affy)
library(oligo)
library(limma)
library(hgu133plus2cdf)

ae_ds <- read_csv("data/ae_samples_studies_ds.csv")

setwd("data/cel_files/")
cel_files <- list.files()
cel_to_geo <- tibble("cel_file"=cel_files) %>% 
  mutate(geo_accession=str_extract(cel_file, "^[0-9A-Za-z]+"))

cel_w_ds <- cel_to_geo %>% left_join(ae_ds %>% select(geo_accession, ds))

disc_gpl570 <- cel_w_ds %>% filter(str_detect(ds, "gpl570")) %>% 
  pull(cel_file)
disc_gpl96 <- cel_w_ds %>% filter(str_detect(ds, "gpl96")) %>% 
  pull(cel_file)

# read in all the GPL570 data - except for the held-out datasets
data1 <- ReadAffy(filenames=disc_gpl570) 
save(data1, file="../data_disc1.RData")

eset1 <- affy::rma(data1)
save(eset1, file="../eset_disc1.RData")


# read in all the GPL96 data
data2 <- ReadAffy(filenames=disc_gpl96) 
save(data2, file="../data_disc2.RData")

eset2 <- affy::rma(data2)
save( eset2, file="../eset_disc2.RData")


