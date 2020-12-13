# 01b_load_cel_files.R
# 11/22/2020
# 
# Code for loading CEL files and running RMA


library(tidyverse)
library(affy)
library(oligo)
library(limma)
library(hgu133plus2cdf)

# ---- AFFY DATA ---- #
list_f <- read_csv("data/list_f_to_download.csv")
affy_f <- list_f %>% filter(manufacturer == "Affymetrix")  %>%
  select(-manufacturer, -supplementary_file) %>% distinct()

affy_plat <- affy_f %>% 
  group_split(platform)

read_f <- function(my_plat){
  setwd(sprintf("data/00_raw_data/Affymetrix/%s/", my_plat))
  cel_files <- list.files(recursive=TRUE)
  data1 <- ReadAffy(filenames=cel_files)
  # todo - RMA?
  
  # save the file
}




# ---- Agilent data ---- #





