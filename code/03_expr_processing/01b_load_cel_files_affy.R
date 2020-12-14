# 01b_load_cel_files.R
# 11/22/2020
# 
# Code for loading CEL files and running RMA


library(tidyverse)
library(affy)
library(oligo)
#library(limma)
#library(hgu133plus2cdf)

args <- commandArgs(trailingOnly=TRUE)
idx <- as.numeric(args[[1]])
# ---- AFFY DATA ---- #
list_f <- read_csv("data/list_f_to_download.csv")
affy_f <- list_f %>% filter(manufacturer == "Affymetrix")  %>%
  select(-manufacturer, -supplementary_file) %>% distinct()

affy_plat <- setdiff(affy_f %>% distinct(gpl) %>% pull(gpl), "GPL16384") # remove the miRNA platform
my_plat <- affy_plat[[idx]]
print(my_plat)

setwd(sprintf("data/00_raw_data/Affymetrix/%s/", my_plat))
cel_files <- list.files(recursive=TRUE)
data <- ReadAffy(filenames=cel_files)
print("read")

eset <- affy::rma(data)
# save the file
setwd("../../../../")
save(eset, file=sprintf("data/01_loaded_data/%s.RData", my_plat))



# --- ones that failed --- #
# GPL11532: "hugene11sttranscriptcluster.db"
# GPL16686: "hugene20sttranscriptcluster.db"

