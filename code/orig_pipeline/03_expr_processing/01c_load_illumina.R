library(tidyverse)
#library(lumi)
library(beadarray)
library(illuminaio)

# NOTE: 
#  - the tar files are the .bgx files for the platforms and not the data itself

list_f <- read_csv("data/list_f_to_download.csv")
illumina_f <- list_f %>% filter(manufacturer == "Illumina")  %>%
  select(-manufacturer) 

# read one per study
my_plat <- "GPL6102"
setwd(sprintf("data/00_raw_data/Illumina/", my_plat))
my_f <- list.files(recursive=TRUE)
my_f[!str_detect(my_f, ".tar|.bgx|control")]

example.lumi <- lumiR.batch(fileName, lib.mapping='lumiHumanIDMapping', 
sampleInfoFile='sam')

ctls_tab <- read_tsv("GSE34198/GSE34198_controlsTable.txt.gz")
nn_tab <- read_tsv("GSE34198/GSE34198_non-normalized.txt.gz", comment="#")