library(tidyverse)

load("data/study_sample_info.RData") # --> all_stats
incl_studies <- read_csv("data/incl_studies.csv")
smok_sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
my_studies <- all_stats[incl_studies2$study_acc]

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

lung_kept %>% write_csv("data/lung_sample_labels.csv")
