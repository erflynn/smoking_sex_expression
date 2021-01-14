
# tissue, study_acc*, title*, kept, authors/location*, platform*, 
# notes, tissue_type, metadata_sex_labels, demo_info_avail, overlapping_studies,
# date*, PMID or paper_link*, paper_sex_breakdown, study_sex_lab_qc,
# num_samples*, healthy_smok, healthy_ns, copd, cancer, other_dz

library(tidyverse)
annot_tab <- read_csv("data/supp_tables/supp_table_1_annot.csv")
study_sl <- read_csv("data/smok_study_sex_lab.csv")
sample_metadata_filt <- read_csv("../drug_trt/data/sample_metadata_filt.csv",
                                 col_types="cccccdldcc")

fct_summ <- function(df) {summary(df %>% mutate(across(everything(), as.factor)))}
add_sl <- function(df) {df %>% left_join(sample_metadata_filt %>% filter(label_type=="expression") %>% 
            select(sample_acc, sex_lab))}

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
samples %>% write_csv("list_smok_samples.csv")
# pull date, PMID + paper link, platform 
geo_smok_studies <- read_csv("data/geo_smok_studies_1231.csv") # gpl, submission_date, pubmed_id
#sra_studies <- read_csv("data/sra_studies_1231.csv") # adds no more information
study_dates <- read_csv("../drug_trt/data/data_old/study_dates.csv") # --> submission date
# TODO: should be able to get platform from rb metadata
table(incl_studies2$study_type)

table(incl_studies2$source) # 114 GEO, 14 SRA
incl_studies2 %>% filter(source=="SRA")

# -- pull pheno data tables to get further demographic information -- #
# what info do we already have??
load("data/study_sample_info.RData") # --> all_stats
my_studies <- all_stats[incl_studies2$study_acc]

# go thru
# - remove if it's not:
#  >=4 in each category
#  all disease

# group by tissue (ish)
ae_studies <-incl_studies2 %>% filter(tissue %in% c("airway epithelium"))
head(ae_studies)
# looks like we have pheno data for >31< of these -- not very many
present_ae <- ae_studies$study_acc[(
  sapply(ae_studies$study_acc, function(x) file.exists(sprintf("data/pdata/%s.csv",tolower(x)))))]
# 7/35
present_ae # - we know we are using all of these!

all_f <- toupper(str_replace_all(list.files(sprintf("data/pdata/")), ".csv", ""))
present_s <- sample_metadata_filt %>% 
  separate_rows(study_acc, sep=";") %>%
  filter(study_acc %in% all_f) %>%
  distinct(sample_acc) %>% 
  pull(sample_acc)

prev_ae_samples <- smok_sl %>% filter(study_acc %in% present_ae)  # 485 but in 661 rows

missing_ae <- setdiff(ae_studies$study_acc, present_ae)
my_studies[[missing_ae[[1]]]]


missing_ae
# "SRP082973" - very few non-COPD smokers
fct_summ(my_studies[["SRP082973"]]$df)
ae1 <- my_studies[["SRP082973"]]$df %>% 
  filter(`smoking status` != "copd smoker") %>%
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
  #group_by(`smoking status`, sex_lab) %>%
  #count()

# "SRP067922"  - no female non-smokers
fct_summ(my_studies[["SRP067922"]]$df)
ae2 <- my_studies[["SRP067922"]]$df %>%
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
  #group_by(`smoking status`, sex_lab) %>%
  #count()


# "SRP024274" - looks same as above?
fct_summ(my_studies[["SRP024274"]]$df)
ae3 <- my_studies[["SRP024274"]]$df %>%
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
#  group_by(`smoking status`, sex_lab) %>%
#  count()

#"SRP012656" - all nsclc/adenocarcinoma, non-smokers, female - 6 samples/patient
#fct_summ(my_studies[["SRP012656"]]$df)
# my_studies[["SRP012656"]]$df %>% filter(tissue=="normal") %>%
#   group_by(patient, smoker, sex) %>%
#   count()

#"SRP006676" - looks like pooled samples
#fct_summ(my_studies[["SRP006676"]]$df)


# "GSE8545"  
fct_summ(my_studies[["GSE8545"]]$df)
ae4 <- my_studies[["GSE8545"]]$df %>%
  mutate(smok=case_when(
    `smoking status`== "non-smoker" ~ "NS",
    str_detect(`smoking status`, "copd") ~ "copd",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  filter(smok!="copd") %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
  #group_by(smok, sex) %>%
  #count()


# "GSE64614" , good m/f breakdown
#  trachea, large airway epithelium - all NS
#  SAE - sufficient of all
#  previously used?
fct_summ(my_studies[["GSE64614"]]$df)
#my_studies[["GSE64614"]]$df %>% nrow()
#my_studies[["GSE64614"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
ae5 <- my_studies[["GSE64614"]]$df %>%
  mutate(smok=case_when(
    `smoking status`=="smoker" ~ "s",
    `smoking status`=="nonsmoker" ~ "ns",
    TRUE ~ `smoking status`
  )) %>%
  mutate(tissue2=case_when(
    str_detect(source_name_ch1, "trachea") ~ "trachea epithleium",
    str_detect(source_name_ch1, "small") ~ "small airway epithelium",
    str_detect(source_name_ch1, "large") ~ "large airway epithelium"
  )) %>%
  filter(is.na(`copd status`), smok!="copd") %>%
  mutate(smok=toupper(smok)) %>%  
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
  #add_sl() %>%
  #group_by(tissue2, smok, sex_lab) %>%
  #count()

# "GSE63127"- 11 samples previously used, many others not but likely overlap
nrow(my_studies[["GSE63127"]]$df)
#my_studies[["GSE63127"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
fct_summ(my_studies[["GSE63127"]]$df)
ae6 <- my_studies[["GSE63127"]]$df %>%
  filter(`copd status`!="yes") %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
    add_sl() %>%
    select(study_acc, sample_acc, smok, sex_lab)
  
table(my_studies[["GSE63127"]]$df$`smoking status`)

# "GSE52237" - previously used
fct_summ(my_studies[["GSE52237" ]]$df)
my_studies[["GSE52237"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
ae7 <- my_studies[["GSE52237"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
    add_sl() %>%
    select(study_acc, sample_acc, smok, sex_lab)

# "GSE43939" - previously used
fct_summ(my_studies[["GSE43939" ]]$df)
my_studies[["GSE43939"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
ae8 <- my_studies[["GSE43939" ]]$df %>%
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "GSE43079" - previously used
fct_summ(my_studies[["GSE43079"]]$df)
my_studies[["GSE43079"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()

ae9 <- my_studies[["GSE43079"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "GSE19667" - 89 previously used out of 121
fct_summ(my_studies[["GSE19667"]]$df)
my_studies[["GSE19667"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
ae10 <- my_studies[["GSE19667"]]$df %>%
  mutate(smok=case_when(
  `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
  `smoking status` %in% c("s", "smoker") ~ "S",
  str_detect(`smoking status`, "^smoker") ~ "S",
  TRUE ~ `smoking status`
)) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
nrow(my_studies[["GSE19667"]]$df)

# "GSE16696"  - 24 previously used out of 50
fct_summ(my_studies[["GSE16696"]]$df)
ae11 <- my_studies[["GSE16696"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

my_studies[["GSE16696"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
nrow(my_studies[["GSE16696"]]$df)

# "GSE13933" - 62 previously used out of 87
fct_summ(my_studies[["GSE13933"]]$df)
my_studies[["GSE13933"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
nrow(my_studies[["GSE13933"]]$df)
ae12 <- my_studies[["GSE13933"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "GSE12815" - all former or current, lots dz (mostly dysplasia - and those do not have sex labels?)
# only three healthy
# fct_summ(my_studies[["GSE12815"]]$df)
# my_studies[["GSE12815"]]$df %>% filter(source_name_ch1=="bronchial epithelium") %>%
#   filter(!str_detect(characteristics_ch1, "post-treatment")) %>%
#   select(sample_acc,characteristics_ch1,title) %>%
#   mutate(dx=case_when(
#     str_detect(title, "dysplasia") ~ "dysplasia",
#     str_detect(title, "nocopd") & str_detect(title, "nolungcancer") ~ "healthy",
#     str_detect(title, "copd") & str_detect(title, "nolungcancer") ~ "copd",
#     str_detect(title, "nocopd") & str_detect(title, "lungcancer") ~ "lung cancer",
#     !str_detect(title, "copd") & str_detect(title, "nolungcancer") ~ "healthy",
#     str_detect(title, "copd") & str_detect(title, "lungcancer") ~ "copd;lung cancer",
#     str_detect(title, "lungcancer") ~ "lung cancer",
#     TRUE ~ "healthy"
#   ), smok=case_when(
#            str_detect(characteristics_ch1, "current") ~ "S",
#            str_detect(characteristics_ch1, "former") ~ "FS",
#            TRUE ~ characteristics_ch1)) %>%
#   add_sl() %>%
#   group_by(dx, smok, sex_lab) %>%
#   count()


# "GSE11952" - 66 of 83 present before
fct_summ(my_studies[["GSE11952"]]$df)
ae13 <- my_studies[["GSE11952"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
my_studies[["GSE11952"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
nrow(my_studies[["GSE11952"]]$df)

# "GSE11784" - 110 of 171 present before
fct_summ(my_studies[["GSE11784"]]$df)
my_studies[["GSE11784"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
nrow(my_studies[["GSE11784"]]$df)
ae14 <- my_studies[["GSE11784"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

#>> keep? "GSE108134"  - none present before
#  includes data with pm2.5 exposure, but similar means,etc
fct_summ(my_studies[["GSE108134"]]$df)
#my_studies[["GSE108134"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
#nrow(my_studies[["GSE108134"]]$df)
ae15 <- my_studies[["GSE108134"]]$df %>%
  filter(`smoking status` !="copd smoker",
         `time of serial bronchoscopy (month)`=="m0") %>%
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
  # group_by(`smoking status`, sex_lab) %>%
  # summarize(n=n(), mean_pm25=mean(as.numeric(`30-day mean pm2.5 exposure level`), na.rm=T),
  #           min_pm25=min(as.numeric(`30-day mean pm2.5 exposure level`), na.rm=T),
  #           max_pm25=max(as.numeric(`30-day mean pm2.5 exposure level`), na.rm=T))

# "GSE89809" - mostly former vs never, some asthma/allergic rhinitis
# only 25 healthy controls, and 4 are former smokers
fct_summ(my_studies[["GSE89809"]]$df)
table(my_studies[["GSE89809"]]$df$source_name_ch1) 
#my_studies[["GSE89809"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
#nrow(my_studies[["GSE89809"]]$df)
ae16 <- my_studies[["GSE89809"]]$df %>%
  filter(str_detect(source_name_ch1, "healthy control"), `allergic rhinitis`=="no") %>%
  mutate(smok=case_when(
    smoking=="current" ~ "S",
    smoking=="former" ~ "FS",
    smoking=="never" ~ "NS"
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
  #group_by(smoking, sex, `allergic rhinitis`) %>%
  #count()

# "GSE5372"  - repeated measures response to injury, very small n (4,4,2)
fct_summ(my_studies[["GSE5372"]]$df)
#my_studies[["GSE5372"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
#nrow(my_studies[["GSE5372"]]$df)
ae17 <- my_studies[["GSE5372"]]$df %>% 
  filter(time=="day 0") %>%
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
  #group_by(sex, `smoking status`) %>%
  #count()

# "GSE34450" - plenty of data, possibly duplicate but not repeated?
fct_summ(my_studies[["GSE34450"]]$df)
#my_studies[["GSE34450"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
#nrow(my_studies[["GSE34450"]]$df)
ae18 <- my_studies[["GSE34450"]]$df %>%
  add_sl() %>%
  mutate(tissue=ifelse(str_detect(source_name_ch1, "large"), "large airway", "small airway")) %>%
  filter(!is.na(`smoking status`)) %>%
  mutate(smok=case_when(
    `smoking status` == "nonsmoker" ~ "ns",
    `smoking status`=="smoker" ~ "s",
    TRUE ~ `smoking status`)) %>%
  mutate(smok=toupper(smok)) %>% 
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
  #group_by(smok, tissue, sex_lab) %>%
  #count()

# "GSE18637" - all but 1 is a non-smoker
fct_summ(my_studies[["GSE18637"]]$df)
my_studies[["GSE18637"]]$df %>% group_by(`smoking status`, sex) %>% count()
ae19 <- my_studies[["GSE18637"]]$df %>%
  mutate(smok=ifelse(`smoking status`=="non-smoker","NS", "S"))%>% 
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "GSE20257" - 110 of 135 previously present
# also contains copd
fct_summ(my_studies[["GSE20257"]]$df)
my_studies[["GSE20257"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
nrow(my_studies[["GSE20257"]]$df)
ae20 <- my_studies[["GSE20257"]]$df %>%
  filter(!str_detect(title, "copd")) %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
    add_sl() %>%
    select(study_acc, sample_acc, smok, sex_lab)

# "GSE76327" - not previously listed but likely overlapping
fct_summ(my_studies[["GSE76327"]]$df)
ae21 <- my_studies[["GSE76327"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

#my_studies[["GSE76327"]]$df %>% filter(sample_acc %in% present_s) %>% nrow()
#nrow(my_studies[["GSE76327"]]$df)

# "GSE77659" - looks v similar
fct_summ(my_studies[["GSE77659"]]$df)
ae22 <- my_studies[["GSE77659"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "GSE10135" - looks v similar
fct_summ(my_studies[["GSE10135"]]$df)

ae23 <- my_studies[["GSE10135"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "GSE92662" - waterpipe smoking vs nonsmoker
# --> only keep the nonsmokers?
fct_summ(my_studies[["GSE92662"]]$df)
ae24 <- my_studies[["GSE92662"]]$df %>%
  filter(`smoking status`=="nonsmoker") %>%
  mutate(smok="NS") %>%
  add_sl() %>% 
  select(study_acc, sample_acc, smok, sex_lab)
  #group_by(sex_lab) %>%
  #count()

# "GSE22047"  - looks v similar
fct_summ(my_studies[["GSE22047"]]$df)

ae25 <-my_studies[["GSE22047"]]$df %>%
  mutate(smok=case_when(
    `smoking status` %in% c("ns", "nonsmoker", "non-smoker") ~ "NS",
    `smoking status` %in% c("s", "smoker") ~ "S",
    str_detect(`smoking status`, "^smoker") ~ "S",
    TRUE ~ `smoking status`
  )) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)


ae_df <- bind_rows(list(ae1, ae2, ae3, ae4, ae5, ae6, ae7, ae8, ae9, 
               ae10,  ae11, ae12, ae13, ae14, ae15, ae16, 
               ae17, ae18, ae19, ae20, ae21, ae22, ae23,
               ae24, ae25))
# --------------- BLOOD STUDIES ------------------- #

blood_studies <- incl_studies2 %>% filter(str_detect(tissue, "blood"))
present_blood <- blood_studies$study_acc[(sapply(blood_studies$study_acc, 
                                                 function(x) file.exists(sprintf("data/pdata/%s.csv",x))))]
# 4/23
# END RESULT: "GSE55962", "GSE42057", "GSE20189" :/ 
# LARGER: exclude: GSE27272, GSE73408, GSE7055, GSE50011
# other tissue: ERP110816, GSE112260


#---> KEEP!!! "GSE55962"  - used in other pubs, COPD
my_studies[["GSE55962" ]]$df %>%  
          filter(`disease status`=="healthy") %>% 
  group_by(gender, age, `smoking status`) %>%
  count()

blood_s1 <- my_studies[["GSE55962" ]]$df %>%  
  filter(`disease status`=="healthy") %>%
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# --> KEEP! "GSE42057" copd, 42 controls, 94 smoking
# - though check that former smokers are not included ("smokcignow")
incl_studies2 %>% filter(study_acc=="GSE42057") %>% pull(description)
my_studies[["GSE42057"]]$df %>% 
  filter(source_name_ch1=="control subject pbmc") %>%
  group_by(gender, source_name_ch1, smokcignow) %>% count()
blood_s2 <- my_studies[["GSE42057"]]$df %>% 
  filter(source_name_ch1=="control subject pbmc", smokcignow!="na") %>%
  mutate(smok=case_when(smokcignow=="0" ~ "NS", smokcignow=="1" ~ "S")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

#"GSE103174"  copd, insufficient non-COPD  (3f, 3m S, 1m ns, 9f ns) 
#incl_studies2 %>% filter(study_acc=="GSE103174")
#my_studies[["GSE103174"]]$df %>% filter(copd=="no") %>% group_by(smoking, sex) %>% count()
blood_s3 <- my_studies[["GSE103174"]]$df %>% 
  group_by(smoking, sex) %>% 
  mutate(smok=ifelse(smoking=="yes", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)


# "GSE12585" - all smokers, 10 heavy + 13 light
#incl_studies2 %>%  filter(study_acc=="GSE12585") %>% pull(description)
blood_s4 <- my_studies[["GSE12585"]]$df %>% 
  mutate(smok="S") %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

missing_blood <- setdiff(blood_studies$study_acc, present_blood)
# GSE87072 - all male
blood_s5 <- my_studies[["GSE87072"]]$df %>%
  filter(`sample group`!="moist_snuff_user") %>%
  mutate(smok=ifelse(`sample group`=="non-tobacco_user", "NS", "S")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# GSE27272 - fetal

# GSE18723 - all female
#   table(my_studies[[missing_blood[[3]]]]$df[,c("smoking", "menopausal")])
blood_s6 <- my_studies[["GSE18723"]]$df %>%
  mutate(smok=ifelse(smoking=="yes", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)


# SRP045500 - no males with smoker/non-smoker data, mostly non-smokers
#  my_studies[[missing_blood[[4]]]]$df %>% 
#  group_by(smoker, gender) %>% 
#  count()
blood_s7 <- my_studies[["SRP045500"]]$df %>%
  mutate(smoker !="--") %>%
  mutate(smok=case_when(smoker=="y" ~ "S", smoker=="n" ~ "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)


# GSE94138 - only a couple smokers
#   - possibly repeated measures
blood_s8 <- my_studies[["GSE94138"]]$df_clean %>% 
  filter(treatment=="placebo", smok!="na") %>% 
  mutate(smok=ifelse(smok==0, "NS", "S")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# GSE73408
# - all disease samples (TB, latent TB, or pneumonia)
#  my_studies[["GSE73408"]]$df %>% filter(`clinical group`=="pna") %>% group_by(`current smoker`,sex) %>% count()

# "GSE71220"
# ... only former vs never when we filter for not COPD
#  regardless of statin status
my_studies[["GSE71220"]]$df %>% 
   filter(disease!="copd", `statin user (y/n)`=="n") %>%
   group_by(sex, `smoking status`) %>%
   count() 
blood_s9 <- my_studies[["GSE71220"]]$df %>% 
  filter(disease!="copd", `statin user (y/n)`=="n") %>%
  mutate(smok=case_when(`smoking status`=="never smoked" ~ "NS",
                        `smoking status`=="former smoker" ~ "FS"))  %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)



# "GSE7055" -- all pmale / prostate cancer
# df_wider <- my_studies[["GSE7055"]]$df %>%
#   separate_rows(race, sep=", ") %>%
#   separate(race, into=c("key", "value"), sep=": ") %>%
#   mutate(value=ifelse(is.na(value), key, value),
#          key=ifelse(key==value, "race", key)) %>%
#   distinct() %>%
#   group_by(sample_acc) %>%
#   pivot_wider(names_from=key, values_from=value, values_fn=condense_dedup)
#table(df_wider$`smoking status`)

# "GSE69683"
#  no healthy smokers :/ 
my_studies[["GSE69683"]]$df %>%
  separate(cohort, into=c("asthma", "smok"), sep=", ") %>%
  group_by(asthma, smok, gender) %>%
  count()
blood_s10 <- my_studies[["GSE69683"]]$df %>%
  separate(cohort, into=c("asthma", "smok"), sep=", ") %>%
  filter(asthma=="healthy") %>%
  mutate(smok=ifelse(smok=="non-smoking", "NS", "S")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "GSE54837": very few non-smokers (6)
#my_studies[["GSE54837"]]$df_clean %>% group_by(smok, sex) %>% count()
blood_s11 <- my_studies[["GSE54837"]]$df_clean %>%
  mutate(smok=ifelse(smok=="smoker controls", "S", "NS"))  %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "GSE21862" 
#  - only one female smoker :/ 
blood_s12 <- my_studies[["GSE21862"]]$df %>% 
  filter(benzene=="ctr") %>% 
  mutate(smok=ifelse(smoking=="yes", "S", "NS")) %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)
#  group_by(smoking, gender) %>%
#  count()

#!---> "GSE20189" <--- 
# KEEP :): 10 m never; 11 f never; 14 f current, 13 m current
my_studies[["GSE20189"]]$df %>% filter(`case status`=="control") %>%
 add_sl() %>%
  group_by(`smoking status`, sex_lab) %>%
  count()

blood_s13 <- my_studies[["GSE20189"]]$df %>% 
  filter(`case status`=="control") %>%
  mutate(smok=case_when(`smoking status`=="current" ~ "S",
                        `smoking status`=="former" ~ "FS",
                        `smoking status`=="never" ~ "NS"))  %>%
  add_sl() %>%
  select(study_acc, sample_acc, smok, sex_lab)

# "ERP110816" - SKIN
# only one smoker, repeated measures, actually skin, psioarisis, etanercept
#summary(my_studies[["ERP110816"]]$df %>% 
#  select(sample_acc, sex, `organism part`, `sampling site`, disease, `clinical history`, treatment,
#         `clinical information`, `isolate`, `individual`) %>%
#  mutate(across(everything(), as.factor)))


# "ERP110814" (related to above)
# blood, repeated measures, only one smoker, psioarisis, etanercept
#summary(my_studies[["ERP110814"]]$df %>% 
#         select(sample_acc, sex, `organism part`,  disease, `clinical history`, treatment,
#                 `clinical information`, `isolate`, `individual`) %>%
 #         mutate(across(everything(), as.factor)))

# "GSE112260"
# actually induced sputum macrophages
# only 4 controls (rest asthma or COPD)
#my_studies[["GSE112260"]]$df %>% filter(diagnosis=="control") %>% group_by(`smoking status`) %>% count()

# "GSE13255"   NSCLC
# 3 f smokers, 5 m smokers - mostly quit, some non-smokers
# unclear what the dx/class columns are 
# 12 stage 1 NSCLC, 15 controls
# 18 paired pre/post surgery, signature reduced after tumor removal
# 2 independent sets of hold-out
# exp_filt <- my_studies[["GSE13255"]]$df %>%
#   select(class, dx, gender, smoking, patient.id, description, stage ) %>%
#   separate(description, into=c("data_type", "experiment"), sep=";") %>%
#     select(-data_type) %>%
#   filter(experiment != "before vs after surgery experiment") 
# exp_filt %>% filter(stage=="") %>%
#   group_by(gender, smoking) %>%
#   count() 
# summary(exp_filt%>%
#   mutate(across(everything(), as.factor)))

# "GSE50011": male-only, all HIV+
#incl_studies2 %>% filter(study_acc=="GSE50011") %>% pull(description) 

# "GSE56768": COPD, multiple blood components
#  very few non-smokers
# incl_studies2 %>% filter(study_acc=="GSE56768") %>% pull(description)
 blood_s14 <- my_studies[["GSE56768"]]$df %>%
   mutate(smok=ifelse(diagnosis=="smoker", "S", "NS")) %>%
   add_sl() %>%
   select(study_acc, sample_acc, smok, sex_lab)
#   group_by(source_name_ch1, diagnosis, sex_lab) %>% count()


# "GSE87650": UC vs CD
# filtered for healthy controls but insufficient smokers
#incl_studies2 %>% filter(study_acc=="GSE87650" )
#summary(my_studies[["GSE87650"]]$df   %>% mutate(across(everything(), as.factor)))
#my_studies[["GSE87650"]]$df   %>%
#  filter(`disease state (simplediagnosis)`=="hc") %>%
#  filter(`smoking status` %in% c("never", "current", "0")) %>%
#  mutate(`smoking status`=ifelse(`smoking status`=="0", "never", `smoking status`)) %>%
#  group_by(sex, `smoking status`, `cell type`) %>%
#  count() %>% arrange(desc(n))

 
blood_kept <- bind_rows(list(blood1, blood_s2, blood_s3, blood_s4, blood_s5, blood_s6, blood_s7,
                           blood_s8, blood_s9, blood_s10, blood_s11, blood_s12, blood_s13, blood_s14))

blood_kept2 <- blood_kept %>% select(-smoking, -sex) 
length(unique(blood_kept2$sample_acc))

# --------------- LUNG STUDIES ------------------- #


lung_studies <- incl_studies2 %>% filter(str_detect(tissue, "lung"))
present_lung <- lung_studies$study_acc[(sapply(lung_studies$study_acc, function(x) file.exists(sprintf("data/pdata/%s.csv",x))))]
# 8/30
present_pdat <- lapply(present_lung, function(x) read_csv(sprintf("data/pdata/%s.csv",x)))
# "GSE63882" -- all cell lines :/
# there is enough adenocarcinoma data
#fct_summ(present_pdat[[1]])
#present_pdat[[1]] %>% filter(smok != "FS", cancer_type=="Adenocarcinoma") %>% group_by(metadata_sex, smok) %>% count()

# "GSE43458" 
# 80 adenocarcinomas + 30 normal lung tissues, all normal lung tissues are non-smokers
present_pdat[[2]] %>% filter(cancer=="n") %>% group_by(smok) %>% count()

lung_s1 <- present_pdat[[2]] %>% filter(cancer=="n")  %>%
  rename(sample_acc=gsm, study_acc=gse) %>%
  add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
#fct_summ(present_pdat[[2]])

# "GSE37768"  - copd
#fct_summ(present_pdat[[3]])
# only one female smoker non-copd
lung_s2 <- present_pdat[[3]] %>% filter(copd=="n") %>% 
  rename(sample_acc=gsm, study_acc=gse) %>% add_sl() %>%
 # group_by(sex_lab, smok) %>% count()
  select(sample_acc, study_acc, smok, sex_lab)

#  "GSE1650"  - all SMOKERS, severe vs mild emphysema
#incl_studies2 %>% filter(study_acc=="GSE1650") %>% pull(description)

# "GSE12667"  - all cancer
#fct_summ(present_pdat[[5]])

# "GSE103174" - copd, insufficient healthy per group -- also listed in blood but actually LUNG
#fct_summ(present_pdat[[6]])
present_pdat[[6]] %>% filter(copd=="no") %>% group_by(smok, metadata_sex) %>% count()
lung_s3 <- present_pdat[[6]] %>% filter(copd=="no") %>%
  rename(sample_acc=gsm, study_acc=gse) %>% add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)


#"GSE10072"  - adenocarcinoma
# maybe keep! >=4 f + m per cateogry in healthy
incl_studies2 %>% filter(study_acc=="GSE10072") %>% pull(description)
fct_summ(present_pdat[[7]])
lung_s4 <- present_pdat[[7]] %>% filter(cancer=="n") %>%
  rename(sample_acc=gsm, study_acc=gse) %>% add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
#group_by(smok, metadata_sex) %>% count()


# "GSE43580" -- all cancer, plenty of samples tho
#incl_studies2 %>% filter(study_acc=="GSE43580") %>% pull(description)
#fct_summ(present_pdat[[8]])

missing_lung <- setdiff(lung_studies$study_acc, present_lung)
incl_studies2 %>% filter(study_acc %in% missing_lung) %>% View()
# GSE68896 - fetal

# "SRP078419" - all cancer, SCC
#fct_summ(my_studies[["SRP078419"]]$df)

# "GSE74777"  - all cancer, only 3 never smokers + 11 females
#fct_summ(my_studies[["GSE74777" ]]$df) 

# GSE68793" - all cancer, SCC almost all smokers
#fct_summ(my_studies[["GSE68793" ]]$df)

# "GSE68465" - lung adenocarcinoma
#fct_summ(my_studies[["GSE68465" ]]$df)

# "GSE50081" - all primary NSCLC
#fct_summ(my_studies[["GSE50081" ]]$df)

# "GSE40791" - only 10 normal tissue, all never smokers
#fct_summ(my_studies[["GSE40791"  ]]$df)
lung_s4 <- my_studies[["GSE40791"  ]]$df %>% 
  filter(`disease state`=="normal lung tissue") %>%
  mutate(smok=ifelse(smokingstat=="never", "NS", smokingstat)) %>% add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
  #group_by(gender, smokingstat) %>% 
  #count()

#--> KEEP, sufficient! 
# "GSE40364"  - actually small airway epithelium
fct_summ(my_studies[["GSE40364" ]]$df) 
lung_s5 <- my_studies[["GSE40364"]]$df %>% 
  filter(is.na(`copd status`)) %>% 
  mutate(`smoking status`=ifelse(`smoking status`=="smoker", "s", `smoking status`)) %>%
  mutate(smok=toupper(`smoking status`)) %>%
  add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
  #group_by(`smoking status`, sex_lab)  %>%
  #count()

# "GSE31210"  -- ever vs never, but potentially enough!!! 
#  3 male never, but still is some
fct_summ(my_studies[["GSE31210"]]$df)
lung_s6 <- my_studies[["GSE31210"]]$df %>%
  filter(tissue=="normal lung") %>%
  mutate(smok=ifelse(`smoking status`=="ever-smoker", "S", "NS")) %>%
  add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
  #group_by(`smoking status`, sex_lab) %>%
  #count()


# "GSE29013" - all NSCLC
#fct_summ(my_studies[["GSE29013"]]$df)


# "GSE17475"  - all tumor, lung adenocarcinoma
#fct_summ(my_studies[["GSE17475"]]$df)


# "ERP112291" - all tumor (adenocarcinoma or SCC) 
#fct_summ(my_studies[["ERP112291"]]$df)

# "GSE83227"  -- all NS for normal tissue
# fct_summ(my_studies[["GSE83227"]]$df)
lung_s7 <- my_studies[["GSE83227"]]$df %>% 
   filter(diseasestate=="normal") %>%
   mutate(smok="NS") %>%
   add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
   #group_by(sex_lab, normalizedsmokinghistoryinpackyears) %>%
   #count()

# "GSE61628"  - only 4 samples, actually human bronchiolar epithelial cells 
# exposed to lentivirus
#fct_summ(my_studies[["GSE61628"]]$df)

#> "SRP074349" - all NSCLC tissue, but fair number of smokers/non-smokers + m/f
#fct_summ(my_studies[["SRP074349"]]$df)
#my_studies[["SRP074349"]]$df %>% group_by(gender, smoking) %>% count()


# "GSE102287" - very few non-smokers, tons of smokers
# fct_summ(my_studies[["GSE102287"]]$df)
lung_s8 <- my_studies[["GSE102287"]]$df %>%
  mutate(pkyrs=ifelse(is.na(`smoking pack year`), `smoking pack years`, `smoking pack year`),
         tissue2=ifelse(is.na(tissue), `tumor_normal status`, tissue)) %>%
  filter(tissue2 %in% c("normal lung tissue", "n")) %>%
  mutate(smok=ifelse(pkyrs==0, "NS", "S")) %>%
  add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
  #group_by(gender, smok) %>%
  #count()
  

# "GSE110780" - no female current smokers, but lots of former/never
# actually bronchiolalveolar LAVAGE, all patients diseased
#fct_summ(my_studies[["GSE110780"]]$df)
#my_studies[["GSE110780"]]$df %>% 
#  filter(treatment=="untreated") %>%
#  group_by(sex, `smoking status`) %>%
#  count()

# "GSE31852" - all cancer, very few never smokers, no current, lots of former
#fct_summ(my_studies[["GSE31852"]]$df )
#my_studies[["GSE31852"]]$df  %>% group_by(gender, smoking_status) %>% count()

# "GSE33072"  - same as previous
#fct_summ(my_studies[["GSE33072"]]$df )

#> maybe? "GSE31908" - few current, multiple former, no male never
fct_summ(my_studies[["GSE31908" ]]$df)
lung_s9 <- my_studies[["GSE31908" ]]$df %>%
  filter(source_name_ch1=="normal lung") %>%
  mutate(smoking=case_when(
    is.na(smoking) ~ `smoking history`,
    TRUE ~ smoking)) %>%
  mutate(smok=case_when(smoking=="current" ~"S",
                        smoking=="never" ~ "NS",
                        smoking=="former" ~ "FS")) %>%
  add_sl() %>%
    select(sample_acc, study_acc, smok, sex_lab)
  #group_by(smoking, sex) %>%
  #count() 

#> KEEP "GSE32539"  - interstitial pneumonias (IIPs)
# - note it says 50 controls? but there appear to be 100
fct_summ(my_studies[["GSE32539" ]]$df)
lung_s10 <- my_studies[["GSE32539" ]]$df %>%
  filter(`final diagnosis`=="control") %>%
  filter(`smoking status`!="--") %>%
  mutate(smok=case_when(`smoking status`=="current" ~ "S",
                        `smoking status`=="former" ~ "FS",
                        `smoking status`=="nonsmoker"~"NS")) %>%
  add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
  #group_by(`smoking status`, gender) %>%
  #count()

# KEEP "GSE32867" 
fct_summ(my_studies[["GSE32867"]]$df)
lung_s11 <- my_studies[["GSE32867"]]$df %>% 
  filter(tissue=="normal lung") %>%
  mutate(smok=case_when(`smoking status`=="current" ~ "S",
                        `smoking status`=="former" ~ "FS",
                        `smoking status`=="nonsmoker"~"NS")) %>%
  add_sl() %>%
  select(sample_acc, study_acc, smok, sex_lab)
  #group_by(`smoking status`, gender) %>%
  #count()

lung_kept <- bind_rows(list(lung_s1, lung_s2, lung_s3, lung_s4, lung_s5, lung_s6, lung_s7,
               lung_s8, lung_s9, lung_s10, lung_s11))

# --------------- OTHER STUDIES ------------------- #
bronchial <- incl_studies2 %>% filter(tissue %in% c("bronchial epithelium or brushing"))
intersect(bronchial$study_acc, all_f) # 5 previously done
missing_bronchial <- setdiff(bronchial$study_acc, all_f)

# "GSE4635"  - appears all male? double check study info 
fct_summ(my_studies[["GSE4635"]]$df)
my_studies[["GSE4635"]]$df %>%
  mutate(smok=case_when(str_detect(source_name_ch1, "current") ~ "S",
                   str_detect(source_name_ch1, "never") ~ "NS" )) %>%
  add_sl() %>%
  group_by(smok, sex, sex_lab) %>%
  count()
         

# "GSE10038" -- all HIV1+, all smokers, all early copd
#fct_summ(my_studies[["GSE10038"]]$df)

#M "GSE103888" - mostly SCC or adenocarcinoma
# only 4 controls, all non-smokers 
fct_summ(my_studies[["GSE103888"]]$df)
#incl_studies2 %>% filter(study_acc=="GSE103888") %>% pull(description)
#incl_studies2 %>% filter(study_acc=="GSE103888") %>% pull(title)
#my_studies[["GSE103888"]]$df %>% 
#  filter(source_name_ch1=="bal_normal control") %>% 
#  group_by(gender, `smoking history`) %>%
#  count()

# "GSE114489" - bronchial dysplasia that progressed to SCC
# ALMOST all are smokers
# can either define as "normal" based on initial score -or- group 4
#
# baseline biopsies from 32 bronchial sites that persisted/progressed to 31 regressive sites
# 1=normal, 2=reserve cellhyperplasia, 3=squamous metaplasia without atypia, 4=mild dysplasia, 5=moderate dysplassia, 6=severe dysplasia and 7=carcinoma-in-situ.  
#/tGroup 1 = Persistent bronchial dysplasia (BL >/= histology score 4, FU >/= histology score 4);\t
# Group 2 = Regressive bronchial dysplasia (BL >/= histology score 4, FU </= histology score 2);\t
# Group 3 = Progressive non-dysplasia (BL </= histology score 2, FU >/= histology score 4);\t
# Group 4 = Stable non-dysplasia (BL </= histology score 2, FU </= histology score 2)"
fct_summ(my_studies[["GSE114489"]]$df)
table(as.numeric(my_studies[["GSE114489"]]$df$`smoking pack/yr`)==0) #2 non-smokers
my_studies[["GSE114489"]]$df %>%
  mutate(group=str_extract(title, "^group [0-9]")) %>%
  filter(group=="group 4") %>%
  #filter(`bl ffpe dx` %in% c(1,2) |`bl frozen dx` %in% c(1,2)) %>%
  add_sl() %>%
  mutate(smok=ifelse(`smoking pack/yr`=="0", "ns", "s")) %>%
  group_by(sex_lab, smok) %>%
  count()


# "GSE48798"  - looks POOLED, only 4 samples
#fct_summ(my_studies[["GSE48798" ]]$df)
#my_studies[["GSE48798" ]]$df %>% filter(`lung cancer`=="no") %>%



oral_studies <- incl_studies2 %>% filter(str_detect(tissue, "oral"))
intersect(oral_studies$study_acc, all_f) # 3 already present!
setdiff(oral_studies$study_acc, all_f)

# M"GSE16149" - all female?, s + ns
fct_summ(my_studies[["GSE16149" ]]$df)
my_studies[["GSE16149" ]]$df %>%
  add_sl() %>%
  group_by(sex_lab, `smoking status`) %>%
  count()

#M "GSE42743"  (few S/NS female)
incl_studies2 %>% filter(study_acc=="GSE42743") %>% pull(description)
fct_summ(my_studies[["GSE42743" ]]$df)
my_studies[["GSE42743" ]]$df %>%
  filter(`tu-vs-nl`=="normal") %>%
  group_by(`smoking status`, `gender`) %>%
  count()



#? "GSE26549" - oral cancer
#  ALL come from lesions?? need to be careful about which treatment arm
fct_summ(my_studies[["GSE26549"]]$df %>% select(-description))
my_studies[["GSE26549"]]$df %>% 
  filter(#outcome=="no oral cancer development" &
           `time of biopsy`=="biopsy at baseline") %>%
  group_by(sex, `smoking habits`) %>%
  count()

#M "GSE29330" - ACTUALLY head & neck cancer
# only 5 controls, all NS -- normal mucosa??
incl_studies2 %>% filter(study_acc=="GSE29330") %>% pull(title)
fct_summ(my_studies[["GSE29330" ]]$df)
#my_studies[["GSE29330" ]]$df %>%
#  filter(str_detect(source_name_ch1,"normal")) %>%
 # group_by(smoker, gender) %>%
#  count()

incl_studies2 %>% filter(str_detect(tissue, "head/neck"))

# GSE53355 - all HNSCC tumor
#fct_summ(my_studies[["GSE53355"]]$df %>% select(-description))

nasal <- incl_studies2 %>% filter(str_detect(tissue, "nasal"))
intersect(all_f,nasal$study_acc) # "GSE8987"
setdiff(nasal$study_acc, all_f)

# "GSE73129" - looking at schizophrenia OE + lymphoblast 
# SHOULD BE IN BLOOD + OE?
# no female smokers??
# fct_summ(my_studies[["GSE73129"]]$df)
# my_studies[["GSE73129"]]$df %>%
#   filter(!str_detect(source_name_ch1, "schizophrenia")) %>%
#   mutate(tissue=ifelse(is.na(`cell type`), "OE", `cell type`)) %>%
#   group_by(tissue, gender, smoking) %>%
#   count()


#F "GSE119136" - former vs never
# fct_summ(my_studies[["GSE119136"]]$df)
# my_studies[["GSE119136"]]$df %>%
#   filter(disease == "control") %>%
#   filter(`immune_meds (recorded at time of sample collection)`=="none") %>%
#   group_by(gender, smoking_status) %>%
#   count()

#* "GSE2109" - all cancer
# incl_studies2 %>% filter(tissue %in% c("multiple"))
# fct_summ(my_studies[["GSE2109"]]$df)
# my_studies[["GSE2109"]]$df %>%
#   group_by(gender, `tobacco use `) %>%
#   count()


incl_studies2 %>% filter(str_detect(tissue, "esophagus"))

incl_studies2 %>% filter(tissue %in% c("bladder"))
incl_studies2 %>% filter(str_detect(tissue, "alveolar"))
incl_studies2 %>% filter(str_detect(tissue, "sputum")) # one is done previous

incl_studies2 %>% filter(tissue %in% c("brain"))
incl_studies2 %>% filter(tissue %in% c("colon"))
incl_studies2 %>% filter(tissue %in% c("LCL"))
incl_studies2 %>% filter(tissue %in% c("liver"))
incl_studies2 %>% filter(tissue %in% c("kidney"))
incl_studies2 %>% filter(tissue %in% c("muscle"))

# --- PC PLOT --- #
library(vroom)
compendia_expr <- vroom("compendia_smoking_expr.csv")
rownames(compendia_expr) <- compendia_expr$rid
compendia_expr2 <- compendia_expr %>% select(-rid) %>% as.matrix.data.frame()
rownames(compendia_expr2) <- compendia_expr$rid

head(kept_samples)
kept_s2 <- kept_samples %>% 
  inner_join(incl_studies2 %>% select(study_acc, tissue), by="study_acc") %>%
  filter(!is.na(expression)) %>%
  distinct(sample_acc, expression, tissue)

# ---- deduplicate and identify duplicates ---- #
compendia_expr3 <- compendia_expr2[,!duplicated(t(compendia_expr2))] # -->11599 to 11326

# remove missing
row_sums <- apply(compendia_expr3, 1, function(x) sum(is.na(x)))
col_sums <- apply(compendia_expr3, 2, function(x) sum(is.na(x)))
compendia_expr4 <- compendia_expr3[row_sums==0,]
# figure out which are duplicated...
dup_dat <- data.frame(t(compendia_expr3[1:8,order(compendia_expr3[1,])]))
dup_dat$sample <- rownames(dup_dat)
dup_split <- dup_dat %>% 
  as_tibble() %>% 
  group_by(ENSG00000000003, ENSG00000000005, ENSG00000000419, ENSG00000000457, ENSG00000000460,
           ENSG00000000938, ENSG00000000971, ENSG00000001036) %>% 
  summarize(n=n(), sample=paste(sample, collapse=";")) %>%
  filter(n>=2) %>%
  arrange(desc(n), sample) %>%
  ungroup() %>%
  dplyr::select(n, sample) %>%
  mutate(grp=1:n()) %>%
  select(-n)


meta_sex <- sample_metadata_filt %>% 
  filter(label_type=="metadata") %>%
  dplyr::select(sample_acc, sex_lab) %>%
  dplyr::rename(metadata_sex=sex_lab)
kept_s3 <- kept_s2 %>% filter(sample_acc %in% colnames(compendia_expr4))
kept_s3 %>% separate_rows(tissue, sep=";") %>% group_by(tissue) %>% count() %>%
  arrange(desc(n))

compendia_expr5 <- compendia_expr4[,kept_s3$sample_acc]
save(compendia_expr5, file="data/compendia_pre_pc.RData")


blood_kept3 <- blood_kept2 %>% filter(sample_acc %in% colnames(compendia_expr4))
blood_expr <- compendia_expr4[,blood_kept3$sample_acc]

blood_kept3 %>% 
  left_join(meta_sex) %>%
  group_by(smok, sex_lab) %>%
  count()

blood_kept3 %>%
  left_join(meta_sex) %>%
  filter(metadata_sex != sex_lab & metadata_sex != "unlabeled") %>%
  nrow() # 17 mismatch
# 


lung_kept3 <- lung_kept %>% filter(sample_acc %in% colnames(compendia_expr4))
lung_expr <- compendia_expr4[,lung_kept3$sample_acc]
lung_kept3 %>% left_join(meta_sex)


ae_kept <- ae_df %>% filter(sample_acc %in% colnames(compendia_expr4)) %>% distinct() %>%
  filter(!str_detect(smok, "copd"))
ae_kept2 <- ae_kept %>% 
  group_by(sample_acc) %>% 
  summarize(study_acc=paste(unique(study_acc), collapse=";"),
            sex_lab=paste(unique(sex_lab), collapse=";"),
            smok=paste(unique(smok), collapse=";")) 
ae_kept2 %>%
  left_join(meta_sex) %>%
  filter(metadata_sex != sex_lab,
         metadata_sex != "unlabeled") %>%
  nrow() # 9

ae_kept2 %>% left_join(meta_sex) %>% 
  filter(metadata_sex=="unlabeled") %>% nrow()

ae_kept2 %>%
  left_join(meta_sex) %>%
  filter(metadata_sex == sex_lab,
         metadata_sex != "unlabeled") %>%
  nrow() # 237
# 649
ae_kept2 %>% 
  filter(smok!="FS", sex_lab!="unlabeled") %>%
  left_join(meta_sex) %>% nrow()
length(unique(ae_kept2$study_acc))
  group_by(smok, sex_lab) %>% 
  dplyr::count() %>%
  filter(smok!="FS", sex_lab!="unlabeled")

blood_kept3 %>% 
  #group_by(smok, sex_lab) %>% 
  #dplyr::count() %>%
  filter(smok!="FS", sex_lab!="unlabeled") %>%
  left_join(meta_sex) %>% 
  filter(metadata_sex!="unlabeled" &
           metadata_sex!=sex_lab) %>% nrow()
length(unique(blood_kept3$study_acc))

lung_kept3 %>% 
  #group_by(smok, sex_lab) %>% 
  #dplyr::count() %>%
  filter(smok!="FS", sex_lab!="unlabeled") %>%
  left_join(meta_sex) %>% 
  filter(metadata_sex!="unlabeled" &
           metadata_sex!=sex_lab) %>% nrow()

length(unique(lung_kept3$study_acc))


ae_expr <- compendia_expr4[,ae_kept2$sample_acc]
library(bigpca)
start_time = Sys.time(); 
bp1 <- big.PCA(blood_expr, pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time

start_time = Sys.time(); 
lung_bp1 <- big.PCA(lung_expr, pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time

start_time = Sys.time(); 
ae_bp1 <- big.PCA(ae_expr, pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time

# run PCA
#pcs <- prcomp(t(compendia_expr4), rank.=4, scale=TRUE)
pcs2 <- data.frame(bp1$PCs)
pcs2$smoking <- blood_kept3$smok
pcs2$sex <- blood_kept3$sex_lab
pcs2$study <- blood_kept3$study_acc

ggplot(pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()


ggsave("figures/blood_pcs.pdf")
ggplot(pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=study))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/blood_pcs_study.pdf")

# lung
lung_pcs2 <- data.frame(lung_bp1$PCs)
lung_pcs2$smoking <- lung_kept3$smok
lung_pcs2$sex <- lung_kept3$sex_lab
lung_pcs2$study <- lung_kept3$study_acc

ggplot(lung_pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()


ggsave("figures/lung_pcs.pdf")
ggplot(lung_pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=study))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/lung_pcs_study.pdf")

#  AE
ae_pcs2 <- data.frame(ae_bp1$PCs) 
ae_pcs2$smoking <- ae_kept2$smok
ae_pcs2$sex <- ae_kept2$sex_lab
ae_pcs2$study <- ae_kept2$study_acc

ggplot(ae_pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/ae_pcs.pdf")

ggplot(ae_pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()+
  xlim(-0.03, 0.0)+
  ylim(-0.05, 0.05)
ggsave("figures/ae_pcs_zoom.pdf")



# ggplot(ae_pcs2 %>%
#          filter(sex!="unlabeled", smoking!="FS"), 
#        aes(x=PC1, y=PC2, col=study))+
#   geom_point(alpha=0.5)+
#   theme_bw()
# ggsave("figures/ae_pcs_study.pdf")
# 
# 
# ggplot(ae_pcs2 %>%
#          filter(sex!="unlabeled", smoking!="FS"), 
#        aes(x=PC1, y=PC2, col=study))+
#   geom_point(alpha=0.5)+
#   theme_bw()+
#   xlim(-0.02, 0.015)+
#   ylim(-0.02, 0.03)
# 
# ggsave("figures/ae_pcs_study_zoom.pdf")

# -- REMOVE SAMPLES THAT ARE IDENTICAL IN PC SPACE -- #
# todo - figure out PC space cutoff
ae_pcs2$sample_acc <- rownames(ae_pcs2)
ae_pcs3 <- ae_pcs2 %>% 
  #arrange(PC1, PC2, PC3) %>%
  mutate(across(c(PC1, PC2, PC3), ~round(., digits=4))) %>%
  group_by(PC1, PC2, PC3) %>%
  summarize(study=paste(unique(study), collapse=";"),
            sample_acc=paste(unique(sample_acc), collapse=";"),
            sex=paste(unique(sex), collapse=";"),
            smoking=paste(unique(smoking), collapse=";")) 
ae_pcs3 %>% filter(str_detect(sample_acc, ";")) 

ae_pcs3 %>% filter(str_detect(sample_acc, ";"), 
                   str_detect(sex, ";"), str_detect(smoking, ";")) 
# at a 4 digit PC cutoff rounding, we have no conflicting labels
# 54 pairs+ of samples are identical (!)

ae_pcs3.2 <- ae_pcs3 %>%
  filter(PC1 > -0.03 , 
         PC1 < -0.015,
         PC2 > -0.05, 
         PC2 < 0.05)
# --> 407

ae_pcs3.2$sample_acc <- sapply(ae_pcs3.2$sample_acc, function(x) strsplit(x, ";")[[1]][[1]])

# NOW TRY our model. 