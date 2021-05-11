

library(tidyverse)
source("code/00_utils.R")

load("data/study_sample_info.RData") # --> all_stats
incl_studies2 <- read_csv("data/incl_studies.csv")
smok_sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
my_studies <- all_stats[incl_studies2$study_acc]

incl_studies3 <- incl_studies %>% filter(!tissue %in% c("airway epithelium", "blood", "lung")) # 44

present_other <- incl_studies3$study_acc[(sapply(incl_studies3$study_acc, 
                                                 function(x) file.exists(sprintf("data/pdata/%s.csv",x))))] 
present_other # 15 of these I've already parsed

# ---- bronchial ---- #
# only keep 1!
bronchial <- incl_studies2 %>% filter(tissue %in% c("bronchial epithelium or brushing"))
intersect(bronchial$study_acc, present_other) # 5 previously done

# keep: "GSE7895" 
# other tissue:
#  "GSE8823" - alveolar macrophages - KEEP
# "GSE114489" - lung - do not keep

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
my_studies[["GSE7895"]]$df %>% dplyr::select(-description) %>%
  mutate(tissue=source_name_ch1) %>%
  add_sl() %>%
  group_by(characteristics_ch1, sex_lab) %>%
  count()
# "GSE37147"
fct_summ(my_studies[["GSE37147"]]$df)

# "GSE19027" 
fct_summ(my_studies[["GSE19027"]]$df)


missing_bronchial <- setdiff(bronchial$study_acc, all_f)

# "GSE4635"  - discard - only 8 samples
# fct_summ(my_studies[["GSE4635"]]$df)
# my_studies[["GSE4635"]]$df %>%
#   mutate(smok=case_when(str_detect(source_name_ch1, "current") ~ "S",
#                         str_detect(source_name_ch1, "never") ~ "NS" )) %>%
#   add_sl() %>%
#   group_by(smok, sex, sex_lab) %>%
#   count()


# "GSE10038" -- all HIV1+, all smokers, all early copd
#fct_summ(my_studies[["GSE10038"]]$df)

#M "GSE103888" - mostly SCC or adenocarcinoma
# only 4 controls, all non-smokers 
#fct_summ(my_studies[["GSE103888"]]$df)
#incl_studies2 %>% filter(study_acc=="GSE103888") %>% pull(description)
#incl_studies2 %>% filter(study_acc=="GSE103888") %>% pull(title)
#my_studies[["GSE103888"]]$df %>% 
#  filter(source_name_ch1=="bal_normal control") %>% 
#  group_by(gender, `smoking history`) %>%
#  count()

# "GSE114489" - bronchial dysplasia that progressed to SCC - ACTUALLY LUNG
# ALMOST all are smokers
# can either define as "normal" based on initial score -or- group 4
#
# baseline biopsies from 32 bronchial sites that persisted/progressed to 31 regressive sites
# 1=normal, 2=reserve cellhyperplasia, 3=squamous metaplasia without atypia, 4=mild dysplasia, 5=moderate dysplassia, 6=severe dysplasia and 7=carcinoma-in-situ.  
#/tGroup 1 = Persistent bronchial dysplasia (BL >/= histology score 4, FU >/= histology score 4);\t
# Group 2 = Regressive bronchial dysplasia (BL >/= histology score 4, FU </= histology score 2);\t
# Group 3 = Progressive non-dysplasia (BL </= histology score 2, FU >/= histology score 4);\t
# Group 4 = Stable non-dysplasia (BL </= histology score 2, FU </= histology score 2)"
# fct_summ(my_studies[["GSE114489"]]$df)
# table(as.numeric(my_studies[["GSE114489"]]$df$`smoking pack/yr`)==0) #2 non-smokers
# my_studies[["GSE114489"]]$df %>%
#   mutate(group=str_extract(title, "^group [0-9]")) %>%
#   filter(group=="group 4") %>%
#   mutate(smok=ifelse(`smoking pack/yr`=="0", "ns", "s")) %>%
#   group_by(smok) %>% count()
#   #filter(`bl ffpe dx` %in% c(1,2) |`bl frozen dx` %in% c(1,2)) %>%
#   add_sl() %>%
# 
#   group_by(sex_lab, smok) %>%
#   count()


# "GSE48798"  - looks POOLED, only 4 samples
#fct_summ(my_studies[["GSE48798" ]]$df)
#my_studies[["GSE48798" ]]$df %>% filter(`lung cancer`=="no") %>%



oral_studies <- incl_studies2 %>% filter(str_detect(tissue, "oral")) # 7
intersect(oral_studies$study_acc, present_other) # 3 already present!
setdiff(oral_studies$study_acc, all_f)

# "GSE8987" - buccal mucosa and nasal epithelium
fct_summ(my_studies[["GSE8987" ]]$df)

# "GSE40013" - stress, whole saliva
fct_summ(my_studies[["GSE40013" ]]$df)
my_studies[["GSE40013" ]]$df %>% group_by(`chronic stress level`, `smoking status`) %>% count()


# "GSE17913" - buccal mucosa, keep
fct_summ(my_studies[["GSE17913" ]]$df)

# M"GSE16149" - all female?, s + ns buccal mucosa
# keep: smoking study
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
fct_summ(my_studies[["GSE26549"]]$df %>% dplyr::select(-description))
my_studies[["GSE26549"]]$df %>% 
  filter(#outcome=="no oral cancer development" &
    `time of biopsy`=="biopsy at baseline") %>%
  group_by(sex, `smoking habits`) %>%
  count()

#M "GSE29330" - ACTUALLY head & neck cancer
# only 5 controls, all NS -- normal mucosa??
#incl_studies2 %>% filter(study_acc=="GSE29330") %>% pull(title)
#fct_summ(my_studies[["GSE29330" ]]$df)
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

# most barret's esophagus, otherwise former vs current (only 3 never)
#fct_summ(my_studies[["GSE77563"]]$df %>% filter(`clinical diagnosis`=="none"))

# all from tumors
#fct_summ(my_studies[["GSE42363"]]$df)

# all tumors? 
#fct_summ(my_studies[["GSE36223"]]$df %>% filter(dysplasia == "negative"))

incl_studies2 %>% filter(tissue %in% c("bladder"))
# only 4 per grp after filtering for non-malignant samples
#fct_summ(my_studies[["GSE21142"]]$df %>% filter(type=="non-malignant"))

# all tumor
#fct_summ(my_studies[["GSE93527"]]$df)

# all tumor
#fct_summ(my_studies[["GSE31684"]]$df)


incl_studies2 %>% filter(str_detect(tissue, "alveolar"))

# part AE, part AM
fct_summ(my_studies[["GSE13896"]]$df)

# all AM
fct_summ(my_studies[["GSE2125"]]$df)


# add: "GSE8823"  - actually alveolar macrophages from BL
fct_summ(my_studies[["GSE8823"]]$df)


incl_studies2 %>% filter(str_detect(tissue, "sputum")) # one is done previous

#fct_summ(my_studies[["GSE46903"]]$df) # actually blood, too many types of cells

#fct_summ(my_studies[["GSE54837"]]$df) # also blood, all COPD

incl_studies2 %>% filter(tissue %in% c("brain"))
#fct_summ(my_studies[["SRP115956"]]$df %>% filter(phenotype=="ctrl")) # multiple different brain parts

fct_summ(my_studies[["GSE44456"]]$df %>% filter(phenotype=="control"))

incl_studies2 %>% filter(tissue %in% c("colon"))

fct_summ(my_studies[["SRP135671"]]$df) # look up?

#fct_summ(my_studies[["GSE36807"]]$df) # no - all uc/cd

incl_studies2 %>% filter(tissue %in% c("LCL"))
fct_summ(my_studies[["GSE36868"]]$df %>% filter(str_detect(source_name_ch1, "control")))


incl_studies2 %>% filter(tissue %in% c("liver"))
#fct_summ(my_studies[["ERP109255"]]$df) # all some sort of liver disease


#fct_summ(my_studies[["GSE9166"]]$df) # incubated with something


incl_studies2 %>% filter(tissue %in% c("kidney"))
fct_summ(my_studies[["GSE46699"]]$df %>% filter(tissue =="normal"))

incl_studies2 %>% filter(tissue %in% c("muscle"))
fct_summ(my_studies[["GSE18732"]]$df %>% filter(str_detect(title, "normal")))

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





