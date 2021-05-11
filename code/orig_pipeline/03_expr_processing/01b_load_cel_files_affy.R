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

# // read in differently depending on whether the platform is different
if (my_plat %in% c("GPL6244", "GPL11532", "GPL16686")){
  data <- oligo::read.celfiles(filenames=cel_files)
  print("read")
  eset <- oligo::rma(data)
  print("rma run")
} else {
  data <- affy::ReadAffy(filenames=cel_files)
  print("read")
  eset <- affy::rma(data)
  print("rma run")
}

# save the file
setwd("../../../../")
save(eset, file=sprintf("data/01_loaded_data/%s.RData", my_plat))



# --- ones that failed --- #

#oligo::read.celfiles

# The affy package can process data from the Gene ST 1.x series of arrays,
# but you should consider using either the oligo or xps packages, which are specifically
# designed for these arrays.
# (then couldn't load CDF)
# GPL11532: "hugene11sttranscriptcluster.db"

#Error: 
#  The affy package is not designed for this array type.
#Please use either the oligo or xps package.
# GPL16686: "hugene20sttranscriptcluster.db"

# GPL6244

