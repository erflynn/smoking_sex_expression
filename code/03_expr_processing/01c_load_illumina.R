library(tidyverse)
library(beadarray)


list_f <- read_csv("data/list_f_to_download.csv")
illumina_f <- list_f %>% filter(manufacturer == "Illumina")  %>%
  select(-manufacturer) 

# read one per study