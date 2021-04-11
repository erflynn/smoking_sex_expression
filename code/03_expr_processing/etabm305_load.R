# Code for loading ETABM305
# 
# - ran after downloading
# for i in {1..13};
# do
# wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/TABM/E-TABM-305/E-TABM-305.raw.${i}.zip;
# done



library('lumi')
library('tidyverse')
BLOOD.DIR <- "data/blood"
# fix file columns

rename_cols <- function(batch, id) {
  my_df <- read_tsv(sprintf("%s/batch%s/%s.raw.csv", BLOOD.DIR, batch, id), 
                    col_types="ccccccccc") %>% 
    select(-`CompositeSequence Identifier`)
  
  # - change null to zero and then make the detection pvalue 0
  my_df2 <- my_df %>% 
    mutate(across(contains("BeadStudio"), as.numeric)) %>%
    mutate(`BeadStudio:Detection`=ifelse(
      is.na(`BeadStudio:AVG_Signal`) |
        is.na(`BeadStudio:BEAD_STDEV`), 0, `BeadStudio:Detection`)) %>%
    mutate(across(everything(), ~replace_na(., 0)))
  
  # update column names
  colnames(my_df2) <- str_replace_all(colnames(my_df2) ,
                                      "BeadStudio:", sprintf("%s.", id))
  my_df2 %>% write_tsv(sprintf("%s/batch%s/renamed/%s.csv", BLOOD.DIR, batch, id))
}

args <- commandArgs(trailingOnly=TRUE)
batch <- args[1]
BATCH.DIR <- sprintf("%s/batch%s/", BLOOD.DIR, batch)
PROC.DIR <- sprintf("%s/renamed", BATCH.DIR)
dir.create(BATCH.DIR)
dir.create(PROC.DIR)
unzip(sprintf("%s/E-TABM-305.raw.%s.zip", BLOOD.DIR, batch), exdir=BATCH.DIR)

list_samples <- str_replace_all(
  list.files(path=BATCH.DIR, pattern="*.csv"), ".raw.csv", "")
renamed <- sapply(list_samples, function(id) rename_cols(batch, id))

fList <- sapply(list_samples, function(x) 
  sprintf("%s/%s.csv", PROC.DIR, x))
myL <- lumiR.batch(fList, 
                   columnNameGrepPattern=list(beadNum=NA), 
                   inputAnnotation=F)
save(myL, file=sprintf("%s/processed/proc%s.RData", BLOOD.DIR, batch))
