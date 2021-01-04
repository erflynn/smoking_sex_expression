# Author: Annie Chang
### phenotype processing for lung data ###
# rename the columns/fields so we can start to put it together
# for lung specifically it is important to retain cancer information + subject IDs because many of the samples are paired cancer/non-cancer samples


# TODO: need to exclude some of these 


library(tidyverse)

lung_pheno <- load("data/lung_pheno_data.RData")

# GSE10072 
gse10072 <- lung_tabs["GSE10072"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`Cigarette Smoking Status`,
                'age'=`Age at Diagnosis`,
                'metadata_sex'=Gender,
                "cancer_stage"=Stage) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         cancer=ifelse(source_name_ch1=="Adenocarcinoma of the Lung", "y", "n"),
         cancer_type=ifelse(cancer=="y", "adenocarcinoma", NA),
         smok=case_when(
           smok=="Never Smoked" ~ "NS",
           smok=="Current Smoker" ~ "S",
           smok=="Former Smoker" ~ "FS"
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]][[2]]) %>%
  select(-source_name_ch1, -title, -description)
gse10072 %>% 
  mutate(across(c(smok, metadata_sex, cancer, cancer_stage, tissue, cancer_type,
                  tissue2), as.factor),
                    age=as.numeric(age)) %>%
  summary()

# todo - remove tissue2
gse10072.2 <- gse10072 %>% select(-tissue2)
write.csv(gse10072, 'data/pdata/gse10072.csv')

# GSE103174
gse103174 <- lung_tabs["GSE103174"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`smoking`,
                'age'=`age`,
                'lung_t_cells' = 'lung t cells',
                'lung_macrophages' = 'lung macrophages',
                'lung_monocytes' = 'lung monocytes',
                'pack_years' = 'pack years',
                'fev1_fvc' = 'fev1/fvc',
                'metadata_sex'=Sex) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         smok=case_when(
           smok=="no" ~ "NS",
           smok=="yes" ~ "S",
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]]) %>%
  select(-source_name_ch1, -title, -description)

gse103174 %>% mutate(across(c(smok, metadata_sex, 
                              copd, tissue, tissue2, batch, id, group), as.factor),
                    across(c(age, bmi, fev1_fvc, fev1, dlco, pack_years, lung_t_cells,
                             lung_macrophages,lung_monocytes, bmi), as.numeric))   %>%
  summary() 

gse103174 %>% ungroup() %>% distinct(group, smok) 
# make sure that FS --> FS, assuming S+CS are the same
# remvoe "group", "tissue2"
gse103174.2 <- gse103174 %>%
  mutate(smok=ifelse(group=="FS", "FS", smok)) %>%
  select(-group, -tissue2)

write.csv(gse103174.2, 'data/pdata/gse103174.csv')

# GSE12667
lung_tabs["GSE12667"][[1]] %>%
  mutate(tissue="lung")  %>% group_by(`Smoking_Status`) %>% count()
gse12667 <- lung_tabs["GSE12667"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`Smoking_Status`,
                'race'=`Race`,
                'cancer_type' = 'Pathology',
                'metadata_sex'=Gender) %>%
                #"BMI"=bmi) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         #cancer=ifelse(source_name_ch1=="Lung", "y", "n"), # none isnt cancer?
         metadata_sex = case_when(
           metadata_sex == "f" ~ "female",
           metadata_sex == "m" ~ "male",
         ),
         smok=case_when(
           smok=="Former" ~ "FS",
           smok=="Current" ~ "S",
           smok=="Never" ~ "NS" # not available for 24
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]]) %>%
  select(-source_name_ch1, -title, -description)

gse12667 %>% mutate(across(c(smok, metadata_sex, race, cancer), as.factor))
#gse12667 <- gse12667[, -c(7, 9, 10:83)] # make more robust by using contains
gse12667.2 <- gse12667 %>% select(-contains("TSP_Patient"))
gse12667 %>% pivot_longer(contains("TSP_Patient"), names_to="TSP_Patient", values_to="number") %>%
  select(gsm, TSP_Patient, number) %>%
  filter(!is.na(number)) # this doesn't seem to work --> have to redownload and check
geo_gse12667 <- getGEO("GSE12667")
additional_pheno <- pData(geo_gse12667[[1]])
patient_id <- additional_pheno %>% 
  select(geo_accession, characteristics_ch1.2) %>%
  dplyr::rename(patient_id=characteristics_ch1.2) %>%
  mutate(patient_id=str_extract(patient_id, "[0-9]+"))
gse12667.3 <- gse12667.2 %>% 
  left_join(patient_id, by=c("gsm"="geo_accession")) %>%
  mutate(race_ethnicity=case_when(
    race=="African American" ~ "black",
    race=="White" ~ "white"
  )) %>%
  select(-race, -tissue2)
gse12667.3 %>% mutate(across(everything(), as.factor)) %>% summary()
write.csv(gse12667.3, 'data/pdata/gse12667.csv')
  
# GSE1650 - note: very little phenotypic detail # Ask Emily what to do about this?
#  exclude this study: all smokers
gse1650 <- lung_tabs["GSE1650"][[1]] %>%
  #dplyr::rename('tissue' = `source_name_ch1`) %>%
  mutate(group=str_extract(title, "[A-Z]"),
         id=str_extract(title, "[0-9]+"),
         smok="S",  # appears all smokers?
         copd="y", # how did you get this?
         tissue="lung")
table(gse1650$group) # 12 w severe emphysema, 18 with mild/no
#L  N 
#18 12 
gse1650.2 <- gse1650 %>%
  mutate(severity=case_when(group=="L" ~ "mild",
                            group=="N" ~ "severe")) %>%
  select(gsm, gse, gpl, tissue, smok, severity, copd, id)

gse1650.2 %>% mutate(across(everything(), as.factor)) %>% summary()

#geo_gse1650 <- getGEO("GSE1650")
#pData(geo_gse1650[[1]]) # re-analyzed by GSE60486


write.csv(gse1650.2, 'data/pdata/gse1650.csv')


# GSE32537
gse32537 <- lung_tabs["GSE32537"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`smoking status`,
                'age'=`age`,
                'tissue_source' = `tissue source`,
                'years_since_quitting' = 'quit how many years ago',
                'pack_years' = 'pack years',
                'final_diagnosis' = 'final diagnosis',
                'fvc_prebronchodilator_perc_predicted' = 'fvc pre-bronchodilator % predicted',
                'dlco_perc_predicted' = 'dlco % predicted',
                'st_george_total_score' = "st. george's total score",
                'metadata_sex'=gender) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         smok=case_when(
           smok=="nonsmoker" ~ "NS",
           smok=="former" ~ "FS",
           smok=="current" ~ "S"
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]]) %>%
  select(-source_name_ch1, -title, -description, -tissue2)
gse32537.1 <- gse32537 %>% mutate(across(c(smok, metadata_sex, tissue, id, final_diagnosis,
                             repository, tissue_source, preservative), as.factor),
                    across(c(age, rin,  st_george_total_score), as.numeric)) # no warning
gse32537.2 <- gse32537.1 %>% mutate(across(c(pack_years, years_since_quitting, 
                               fvc_prebronchodilator_perc_predicted, dlco_perc_predicted), as.numeric))
# cause warnings -- I have checked, all are fine and should be converted to NAs
gse32537.2 %>%
  summary()

write.csv(gse32537.2, 'data/pdata/gse32537.csv')

# GSE37768
gse37768 <- lung_tabs["GSE37768"][[1]] %>%
  dplyr::rename('smok'=`phenotype`) %>%
  mutate(smok=case_when(
    smok=="Nonsmoker" ~ "NS",
    smok=="healthy smoker" ~ "S",
    smok=="moderate COPD" ~ "copd")) %>%
  mutate(copd=ifelse(smok=="copd", "y", "n")) %>%
  select(-title, -description, -source_name_ch1) %>%
  dplyr::rename(tissue_source=tissue,
                tissue=tissue2)   

gse37768 %>%
  mutate(across(everything(), as.factor)) %>%
  summary()
# advise against hardcoding
#gse37768$copd <- 'n'
#gse37768[21:38, 10] <- 'y'

write.csv(gse37768, 'data/pdata/gse37768.csv')

# GSE43458
gse43458 <- lung_tabs["GSE43458"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'="smoking status") %>%
  mutate(cancer=ifelse(str_detect(title, "adenocarcinoma"), "y", "n"),
         cancer_type=ifelse(cancer=="y", "adenocarcinoma", NA),
         smok=case_when(
           smok=="Never-smoker" ~ "NS",
           smok=="Smoker" ~ "S",
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]][[2]]) %>%
  select(-source_name_ch1, -title, -description, -tissue2, -histology)
gse43458 %>% mutate(across(c(smok, cancer, tissue, cancer_type, id), as.factor)) %>%
  summary()

write.csv(gse43458, 'data/pdata/gse43458.csv')

# GSE43580
# -- exclude: these are all cancer!
gse43580 <- lung_tabs["GSE43580"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`smoking status`,
                'age'= "age at excision (years)",
                'metadata_sex' = gender,
                'clinical_diagnosis_pt'=`clinical diagnosis patient`,
                'smok_dose_cigday' = "smoking dose (cigarettes/day)",
                'smok_dur_yrs' = "smoking duration (years)",
                'height_cm' = "height (cm)",
                "weight_kg" = "weight (kg)",
                "excision_year" = "excision year",
                "sample_recovery_type" = "sample recovery type", 
                "tnm_stage" = "tnm stage", 
                "ajcc_uicc_stage" = "ajcc uicc stage",
                "clinical_diagnosis_specimen" = "clinical diagnosis specimen",
                "sbv_challenge_gold_standard" = "sbv challenge gold standard",
                "biosample_confirmed_diagnosis" = "biosample confirmed diagnosis",
                "biosample_confirmed_subdiagnosis" = "biosample confirmed sub-diagnosis",
                "tumor_grade" = "tumor grade") %>%
  mutate(metadata_sex=tolower(metadata_sex),
         smok=case_when(
           smok=="Occasional Use" ~ "S",
           smok=="Previous Use" ~ "FS",
           smok=="Current Use" ~ "S",
           smok=="Never Used" ~ "NS"
         )) %>%
  group_by(gsm) %>%
  select(-source_name_ch1, -title, -description)
gse43580.1 <- gse43580 %>% mutate(across(c(smok, metadata_sex, clinical_diagnosis_pt,
                             sample_recovery_type, ethnicity, tumor_grade, tnm_stage, ajcc_uicc_stage,
                             clinical_diagnosis_specimen, sbv_challenge_gold_standard, biosample_confirmed_diagnosis,
                             biosample_confirmed_subdiagnosis), as.factor),
                    across(c(age, weight_kg, height_cm, bmi, smok_dur_yrs, smok_dose_cigday,
                             excision_year, rin), as.numeric)) %>%
  select(-tissue) %>%
  dplyr::rename(tissue=tissue2)
gse43580.1 %>%
  summary()

write.csv(gse43580.1, 'data/pdata/gse43580.csv')

# GSE63882
# ?? exclude: cell line + all cancer 
gse63882 <- lung_tabs["GSE63882"][[1]] %>%
  mutate(tissue="lung") %>%
  select(-tissue2, -description) %>%
  dplyr::rename('smok'=`smoker`,
                'metadata_sex' = "gender",
                cancer_type=histology,
                race_ethnicity=race,
                source=source_name_ch1,
                pack_years=`pack-years`) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         race_ethnicity=tolower(race_ethnicity),
         # cancer_type="adenocarcinoma", # not all are this
         metadata_sex = case_when(
           metadata_sex == "f" ~ "female",
           metadata_sex == "m" ~ "male",
         ),
         smok=case_when(
           smok=="Y" ~ "S",
           smok=="N" ~ "NS",
           smok=="Ex-smoker" ~ "FS"
         )) %>%
  group_by(gsm) %>%
 # mutate(id=str_split(title, "_")[[1]][[2]]) %>%
  select( -title)

gse63882.1 <- gse63882 %>% mutate(across(c(smok, metadata_sex, cancer_type,race_ethnicity, source, tissue, egfr, kras), as.factor),
                    across(c(age, pack_years), as.numeric))
gse63882.1 %>% summary()
#gse63882$cancer <- 'y'
#gse63882 <- gse63882[, -c(15)]
#colnames(gse63882)[11] <- c("pack_years")
#colnames(gse63882)[12] <- c("cancer_type")

write.csv(gse63882.1, 'data/pdata/gse63882.csv')


