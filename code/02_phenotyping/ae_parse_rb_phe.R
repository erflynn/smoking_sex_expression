

library(tidyverse)

load("data/study_sample_info.RData") # --> all_stats
incl_studies <- read_csv("data/incl_studies.csv")
smok_sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
my_studies <- all_stats[incl_studies2$study_acc]

# group by tissue (ish)
ae_studies <-incl_studies2 %>% filter(tissue %in% c("airway epithelium"))
# looks like we have pheno data for >31< of these -- not very many
present_ae <- ae_studies$study_acc[(
  sapply(ae_studies$study_acc, function(x) 
    file.exists(sprintf("data/pdata/%s.csv",tolower(x)))))]
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


ae_df %>% write_csv("data/ae_sample_labels.csv")

# TODO - add from this one
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

