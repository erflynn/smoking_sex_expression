
# tissue, study_acc*, title*, kept, authors/location*, platform*, 
# notes, tissue_type, metadata_sex_labels, demo_info_avail, overlapping_studies,
# date*, PMID or paper_link*, paper_sex_breakdown, study_sex_lab_qc,
# num_samples*, healthy_smok, healthy_ns, copd, cancer, other_dz

library(tidyverse)
annot_tab <- read_csv("data/supp_tables/supp_table_1_annot.csv")
study_sl <- read_csv("data/smok_study_sex_lab.csv")
sample_metadata_filt <- read_csv("../drug_trt/data/sample_metadata_filt.csv",
                                 col_types="cccccdldcc")

incl_studies <- annot_tab %>% 
  filter(included=="yes" & study_type=="smoking history")

incl_studies2 <- incl_studies %>% 
  filter(num_samples >= 8,
         tissue != "prostate",
         tissue != "placenta or umbilical cord") %>%
  select(-num_samples, -study_type, -included, -treatment) %>%
  left_join(study_sl %>% filter(labeling_method=="expression"), by="study_acc") %>%
  select(-dataset, -labeling_method) 

incl_studies2 %>% filter(study_type == "male-only" | study_type=="female-only" | study_type=="unlabeled")

incl_studies2 %>% write_csv("data/incl_studies.csv")

list_cols <- read_tsv("list_cols_idx.tsv")
smok_sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
kept_samples <- smok_sl %>% filter(study_acc %in% incl_studies2$study_acc)
kept_dat <- list_cols %>% semi_join(kept_samples, by="sample_acc")

samples <- kept_dat %>% rename(acc=sample_acc) %>% select(acc, f_idx, idx) # distinct??
samples %>% write_csv("data/list_smok_samples.csv")
# pull date, PMID + paper link, platform 
geo_smok_studies <- read_csv("data/geo_smok_studies_1231.csv") # gpl, submission_date, pubmed_id
#sra_studies <- read_csv("data/sra_studies_1231.csv") # adds no more information
study_dates <- read_csv("../drug_trt/data/data_old/study_dates.csv") # --> submission date
# TODO: should be able to get platform from rb metadata
table(incl_studies2$study_type)

table(incl_studies2$source) # 114 GEO, 14 SRA

