# 00_phenotyping_utils.R
# 11/20/2020
# 
# Helpful utilities for extracting covariate data
#
# TODO:
# - update to new phenotyping setup
# - make less specific?

library('tidyverse')

clean_col <- function(x){
  if(!str_detect(x, ":")){
    return(x)
  }
  y <- str_split(x, ":")[[1]]
  return(str_trim(y[[2]]))
}

combine_cols <- function(pData){
  # combine columns that might have different cases
  pData %>%
    unite("sex", contains("sex:ch1"), sep=";", na.rm=TRUE) %>%
    unite("ethnic", contains("ethnic group:ch1"), sep=";", na.rm=TRUE) %>%
    unite("smok", contains("smoking status:ch1"), sep=";", na.rm=TRUE) %>%
    unite("age", contains("age:ch1"), sep=";", na.rm=TRUE)
}

check_relation_fields <- function(pheno){
  unique(unique(pheno %>% select(contains("relation")) %>% pull(relation.1)),
         unique(pheno %>% select(contains("relation")) %>% pull(relation)),
         unique(pheno %>% select(contains("relation")) %>% pull(relation.2)))
}

clean_pData <- function(pData){
  pData %>%
    select(geo_accession, source_name_ch1, title, sex, ethnic, smok, age) %>%
    mutate(pkyrs=str_extract(smok,"[0-9|\\.]+"),
           smok=case_when(str_detect(smok, "non-smoker") ~ "NS",
                          str_detect(smok, "smoker") ~ "S",
                          TRUE ~ smok),
           dgm_id=str_extract(title,"[0-9]+")) %>%
    rename(source_name=source_name_ch1)
}

process_SB_pDat1 <- function(pheno_dat){
  pheno_dat %>%
    rename("age"="characteristics_ch1",
           "sex"="characteristics_ch1.1",
           "ethnic"="characteristics_ch1.2",
           "smok"="characteristics_ch1.3") %>%
    select(title, geo_accession, submission_date, source_name_ch1, 
           age, sex, smok, ethnic, 
           `age:ch1`:`Smoking status:ch1`) %>% # also grab the extra cols so we can confirm they are the same
    group_by(geo_accession) %>%
    mutate(across(age:ethnic, ~clean_col(.)))  %>%
    ungroup() %>%
    rename(age0=age, sex0=sex, ethnic0=ethnic, smok0=smok) %>%
    combine_cols() 
}
process_SB_pDat2 <- function(pheno1){
  pheno1 %>%
    clean_pData() %>%  
    group_by(geo_accession) %>%
    mutate(title=str_split(title, ",")[[1]][[1]]) %>%
    ungroup() %>%
    select(-source_name) %>%
    rename(source_name=title)
}
grab_copd_info <- function(pheno){
  pheno %>%
    mutate(copd=ifelse(str_detect(smok, "COPD"), "y", "n")) %>%
    mutate(copd_pheno=case_when(
      copd == "y" ~ str_squish(str_replace_all(str_extract(smok, "[A-z\\-, ]+" ), ",", "")),
      copd == "n" ~ copd
    )) %>%
    mutate(smok=case_when(
      copd == "y" ~ "S",
      copd == "n" ~ smok
    ))
}