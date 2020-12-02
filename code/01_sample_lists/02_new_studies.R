# see if we ID any new studies from the metadata

require('tidyverse')

setwd("../labeling/geo2drug/")
# orig data -- is this all of it?
smok_dat <- read_csv("data/smok_dat/smoking_data.csv")
# 220 -- many duplicated (really it's 199)
smoking_labels <- read_csv("data/smok_dat/smoking_labels_reform.csv")


# ---- microarray ---- #

microarray_data <- read.csv("../../drug_trt/data/01_metadata/human_experiment_metadata.csv", stringsAsFactors = FALSE)
microarray_data <- microarray_data %>% 
  group_by(study_acc) %>% 
  mutate(str=paste(c(title, description), collapse=" ")) %>%
  ungroup()

microarray_smok <- microarray_data %>%
  filter(grepl( "smok|nicotine|tobacco|cigarette", str, 
                ignore.case=TRUE))  # 192

microarray_smok2 <- microarray_smok %>% 
  anti_join(smoking_labels %>% 
              select(gse) %>% 
              unique(), by=c("study_acc"="gse")) %>% 
  select(study_acc, title, description) 
  # 61 (79 w/o full list) 

microarray_smok2 %>% View()


# ---- rnaseq ---- #
rnaseq_data <- read_csv("../../drug_trt/data/01_metadata/human_rnaseq_experiment_metadata.csv")

rnaseq_data <- rnaseq_data %>% 
  group_by(study_acc) %>% 
  mutate(str=paste(c(title, description), collapse=" ")) %>%
  ungroup()

rnaseq_smok <- rnaseq_data %>%
  filter(grepl( "smok|nicotine|tobacco|cigarette", str, ignore.case=TRUE))  # 37

rnaseq_smok2 <- rnaseq_smok %>% 
  anti_join(smok_dat %>% 
              select(gse) %>% 
              unique(), by=c("study_acc"="gse")) %>% 
  select(study_acc, title, description) %>%
  unique() # 37 (no overlap)

microarray_smok2 %>% head()
rnaseq_smok2 %>% head()
intersect(microarray_smok2$study_acc, rnaseq_smok2$study_acc)

new_studies <- microarray_smok2 %>% anti_join(rnaseq_smok2, by=c("study_acc")) %>%
  bind_rows(rnaseq_smok2)

new_studies %>% write_csv("data/smok_dat/new_cersi_studies.csv")

rnaseq_smok2 %>% View()

# write out the info for these new studies!
metadata_new <- rbind(microarray_smok %>% mutate(src="microarray"), 
                      rnaseq_smok %>% mutate(src="rnaseq"))

# NOW read in previous manual annotations

# join and write out
metadata_new %>% write_csv("data/smok_dat/refinebio_metadata_to_annot.csv")

# ---> now time to manually annot

# ---- what other data is available for these samples? ---- #
microarray_exp_to_sample <- read_csv("../../drug_trt/data/01_metadata/human_exp_to_sample.csv")
microarray_sample <- read.csv("../../drug_trt/data/01_metadata/human_metadata.csv", stringsAsFactors = FALSE)
microarray_smok_samp <- microarray_smok %>%
  left_join(microarray_exp_to_sample) %>%
  left_join(microarray_sample, by=c("sample_acc"="acc"))

m_sl <- microarray_smok_samp %>% 
  left_join(labels4 %>% select(gsm, consensus_sex), by=c("sample_acc"="gsm")) %>% 
  select(study_acc, sample_acc, sex,consensus_sex)  %>%
  rename(text_sex=sex, expr_sex=consensus_sex, gse=study_acc, gsm=sample_acc)

alluvPlotSample(m_sl %>% select(-gse) %>% 
                  mutate(text_sex=ifelse(!text_sex %in% c("female", "male"), "unknown", text_sex),
                                expr_sex=ifelse(is.na(expr_sex), "unknown", expr_sex)))


alluvPlotStudy(m_sl %>%
                  mutate(text_sex=ifelse(!text_sex %in% c("female", "male"), NA, text_sex),
                         expr_sex=ifelse(is.na(expr_sex), NA, expr_sex)))



rnaseq_exp_to_sample <- read_csv("../../drug_trt/data/01_metadata/human_exp_to_sample_counts.csv")
rnaseq_sample <- read.csv("../../drug_trt/data/01_metadata/human_rnaseq_sample_metadata.csv", stringsAsFactors = FALSE)
rnaseq_smok_samp <- rnaseq_smok2 %>% 
  left_join(rnaseq_exp_to_sample)  %>% 
  left_join(rnaseq_sample, by=c("sample_acc"="acc"))

# // TODO - clean RNA-seq metadata sex labels :(
# but seems like it's either m, f, or nothing here... so that's ok
# the majority are missing!!!

table(rnaseq_smok_samp$sex) # mostly unlabeled, 20% m/f
table(rnaseq_smok_samp$cl_line) # some beas-2b, no other repeats

# //TODO: read in the sex labels!! - TBD...
