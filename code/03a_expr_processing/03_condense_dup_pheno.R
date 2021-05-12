# 03_condense_dup_pheno.R
#
# Code for condensing the duplicated phenotype data

library('tidyverse')

load("data/eset_comb.RData") # --> eset_comb, dup_split

condense_cols <- function(x){
  paste(sort(unique(x[!is.na(x)])), collapse=";")
}

comb_metadata <- read_csv("data/comb_sample_epithelium_metadata.csv")
dup_dat <- dup_split %>% 
  mutate(grp=1:n()) %>%
  select(-n) %>%
  rename(geo_accession=sample) %>%
  separate_rows(geo_accession) %>% 
  left_join(comb_metadata, by=c("geo_accession")) %>%
  group_by(grp) %>%
  summarize(across(everything(), ~condense_cols(.))) 

# for duplicated columns, put together the phenotype data
dup_dat %>% filter(str_detect(metadata_sex,";") |
                     str_detect(race_ethnicity,";") |
                     str_detect(smok,";") | # one
                     str_detect(age,";")  | #one
                     str_detect(expr_sex, ";"))

dup_dat %>% filter(str_detect(tissue, ";") ) %>% group_by(tissue) %>% count()

dup_dat2 <- dup_dat %>%
  mutate(
    tissue=case_when(
      ! str_detect(tissue, ";")  ~ tissue,
      tissue=="airway epithelium;small airway epithelium" ~ "small airway epithelium",
      tissue=="trachea;trachea epithelium" ~ "trachea epithelium" 
    ), # take the more specific tissue
    age=case_when(
      ! str_detect(age, ";") ~ age,
      TRUE ~ as.character(mean(as.numeric(str_split(age, ";")[[1]])))
    ) # take the mean of the ages
  ) %>%
  filter(!str_detect(smok, ";")) # remove conflicting

# add the phenotype data
present_cols <- colnames(eset_comb)


dup_dat2_long <- dup_dat2 %>% 
  separate_rows(geo_accession, sep=";") %>% select(-grp) %>%
  mutate(across(c(age, pack_years, pred), as.numeric))
comb_metadata2 <- dup_dat2_long %>%
  bind_rows(comb_metadata %>% anti_join(dup_dat2_long, by=c("geo_accession"))  )

comb_metadata3 <- comb_metadata2 %>% 
  filter(geo_accession %in% present_cols) %>%
  arrange(geo_accession) %>%
  mutate(across(everything(), ~ifelse(.=="", NA, .))) %>%
  mutate(tissue=ifelse(tissue=="trachea", "trachea epithelium", tissue)) %>%
  mutate(across(c(smok, tissue, metadata_sex, race_ethnicity,
                  expr_sex), as.factor))

summary(comb_metadata3 %>% select(-id, -geo_accession, -study))
eset_comb2 <- eset_comb[,comb_metadata3$geo_accession]
# remove missing data
row_miss <- apply(eset_comb2, 1, function(x) sum(is.na(x)))
eset_comb3 <- eset_comb2[row_miss==0,]

save(eset_comb3, comb_metadata3, file="data/ae_eset_pheno.RData")