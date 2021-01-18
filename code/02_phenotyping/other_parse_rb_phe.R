

library(tidyverse)
source("code/00_utils.R")

load("data/study_sample_info.RData") # --> all_stats
incl_studies2 <- read_csv("data/incl_studies.csv")
smok_sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
my_studies <- all_stats[incl_studies2$study_acc]



# ---- bronchial ---- #

bronchial <- incl_studies2 %>% filter(tissue %in% c("bronchial epithelium or brushing"))
intersect(bronchial$study_acc, all_f) # 5 previously done

bronchial %>% pull(study_acc)

# "GSE97010" - validation?
# acute cigarette exposure
#fct_summ(my_studies[["GSE97010"]]$df)
#my_studies[["GSE97010"]]$df %>%
#  select(-treatment_protocol_ch1, -description, -source_name_ch1) %>%
# fct_summ()


# "GSE8823"  - actually alveolar macrophages from BL
#fct_summ(my_studies[["GSE8823"]]$df)

# "GSE7895" -- KEEP
fct_summ(my_studies[["GSE7895"]]$df)
my_studies[["GSE7895"]]$df %>% select(-description) %>%
  mutate(tissue=source_name_ch1) %>%
  add_sl() %>%
  group_by(characteristics_ch1, sex_lab) %>%
  count()
# "GSE37147"
# "GSE19027" 




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


# add: "GSE8823"  - actually alveolar macrophages from BL
fct_summ(my_studies[["GSE8823"]]$df)


incl_studies2 %>% filter(str_detect(tissue, "sputum")) # one is done previous

incl_studies2 %>% filter(tissue %in% c("brain"))
incl_studies2 %>% filter(tissue %in% c("colon"))
incl_studies2 %>% filter(tissue %in% c("LCL"))
incl_studies2 %>% filter(tissue %in% c("liver"))
incl_studies2 %>% filter(tissue %in% c("kidney"))
incl_studies2 %>% filter(tissue %in% c("muscle"))


# from other:
# "ERP110816" - SKIN
# only one smoker, repeated measures, actually skin, psioarisis, etanercept
#summary(my_studies[["ERP110816"]]$df %>% 
#  select(sample_acc, sex, `organism part`, `sampling site`, disease, `clinical history`, treatment,
#         `clinical information`, `isolate`, `individual`) %>%
#  mutate(across(everything(), as.factor)))


# "GSE112260"
# actually induced sputum macrophages
# only 4 controls (rest asthma or COPD)
#my_studies[["GSE112260"]]$df %>% filter(diagnosis=="control") %>% group_by(`smoking status`) %>% count()

