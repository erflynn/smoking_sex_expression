
library(tidyverse)

load("data/study_sample_info.RData") # --> all_stats
incl_studies <- read_csv("data/incl_studies.csv")
smok_sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
my_studies <- all_stats[incl_studies2$study_acc]

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

blood_kept2 %>% write_csv("data/blood_sample_labels.csv")
