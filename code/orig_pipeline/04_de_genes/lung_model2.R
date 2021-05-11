
library(tidyverse)
library(limma)
library(MetaIntegrator)
library(sva)


# Load list of samples
s1 = read_csv("data/supp_tables/supp_table_1_annot.csv")
head(s1)
lung_studies = s1 %>%
  filter(str_detect(tissue, 'lung'),
         included=='yes') %>%
  distinct(study_acc) 
# 31!
present_lung <- lung_studies$study_acc[(sapply(lung_studies$study_acc, function(x) 
  file.exists(sprintf("data/pdata/%s.csv",x))))]
table(present_lung)
load("data/study_sample_info.RData") # --> all_stats
# to remove: (incorrect tissue) GSE61628, GSE110780, GSE68896, GSE40364
my_studies <- all_stats[setdiff(lung_studies$study_acc, 
                                c("GSE61628", "GSE110780", "GSE68896", "GSE40364"))]

length(my_studies)
my_studies[[1]]$cols
lung_phe_dat = do.call(rbind, lapply(my_studies, function(x) {
       x$df_clean %>%
        select(study_acc, sample_acc, smok, sex, race_ethnicity, pkyrs, copd)}))
lung_phe_dat %>%
  fct_summ()

# moderate copd - GSE37768
# occasional use
# some smoking
#lung_phe_dat %>% filter(smok== "2") # "1"=current, 2=ex, 3=never "SRP074349"


smok_map <- lung_phe_dat %>%
  group_by(smok) %>%
  count() %>%
  ungroup() %>%
  mutate(smok_norm=case_when(
    smok %in% c("smoked in the past","ever-smoker", "ex-smoker", "ever", "2") ~ "FS",
    str_detect(smok, "former") | 
      str_detect(smok, "quit smoking") | 
      str_detect(smok, "previous use") |
      str_detect(smok, "current reformed") ~ "FS",
    smok %in% c("1", "yes", "smoker", "y", "healthy smoker") ~ "S",
    !is.na(as.numeric(smok)) & ! smok %in% c("0","1","2","3") ~ "S",
    str_detect(smok, "current") |
      str_detect(smok, "from smokers") |
      str_detect(smok, "some smoking") |
      str_detect(smok, "occasional use") ~ "S", 
    str_detect(smok, "never") | 
      str_detect(smok, "non-smoker") ~ "NS",
    smok %in% c("0", "0.0", "nonsmoker", "no", "n", "3") ~ "NS",
    str_detect(smok, "not available") ~ "",
    smok %in% c("--", "unknown", "unable to determine", "n.a.", "not specified", "na", "NA") ~ "",
    is.na(smok) ~ "",
    TRUE ~ ""
  )) %>%
  arrange(desc(n))
smok_map %>% View()

# --> todo normalize other covariates <-- 
## sex + race-ethnicity should be easier
# lung_phe_dat %>%
#   group_by(race_ethnicity) %>%
#   count() %>%
#   arrange(desc(n))

## age??
## pkyrs will be complicated, much of this is in the smoke field

# combine together

# grab the correct smoking data
(counts_per_study=lung_phe_dat %>%
  left_join(smok_map) %>%
  select(study_acc,sample_acc, smok_norm) %>%
  filter(smok_norm  %in% c('S', 'NS')) %>%
  group_by(study_acc,smok_norm) %>%
  count() %>%
  pivot_wider(names_from=smok_norm, values_from=n, values_fill=0))
  
sufficient_samples=counts_per_study %>%
  filter(NS >= 5, S >=5)

#Which of these are included in the expression data I already have downloaded?

### --- 1: Extract extra covariate data --- ###



### --- 2: Adjust with combat and SVA --- ###


### --- 3: Try glmnet --- ###
# study stratifies for CV
# try with and without adjusted data

### --- 4: Try MetaIntegrator --- ###
sufficient_samples$study_acc

### --- 5: Try additional covariates --- ###