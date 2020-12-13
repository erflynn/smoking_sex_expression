# Author: Annie Chang
### phenotype processing for lung data ###
# rename the columns/fields so we can start to put it together
# for lung specifically it is important to retain cancer information + subject IDs because many of the samples are paired cancer/non-cancer samples

lung_pheno <- load("data/lung_pheno_data.RData")

# GSE10072 
gse10072 <- lung_tabs["GSE10072"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`Cigarette Smoking Status`,
                'age'=`Age at Diagnosis`,
                'metadata_sex'=Gender,
                "cancer_stage"=Stage) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         cancer_type="adenocarcinoma",
         cancer=ifelse(source_name_ch1=="Adenocarcinoma of the Lung", "y", "n"),
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
                #"BMI"=bmi) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         cancer_type="lung",
         cancer=ifelse(cancer_type=="lung", "y", "n"),
         smok=case_when(
           smok=="no" ~ "NS",
           smok=="yes" ~ "S",
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]]) %>%
  select(-source_name_ch1, -title, -description)

gse103174 %>% mutate(across(c(smok, metadata_sex, contains("cancer"),
                              copd, tissue, tissue2, batch, id, group), as.factor),
                    across(c(age, bmi, fev1_fvc, fev1, dlco, pack_years, lung_t_cells,
                             lung_macrophages,lung_monocytes), as.numeric))   %>%
  summary() 

gse103174 %>% ungroup() %>% distinct(group, smok) 
# make sure that FS --> FS, assuming S+CS are the same
# remvoe "group", "tissue2"
gse103174.2 <- gse103174 %>%
  mutate(smok=ifelse(group=="FS", "FS", smok)) %>%
  select(-group, -tissue2)

write.csv(gse103174.2, 'data/pdata/gse103174.csv')

# GSE12667
gse12667 <- lung_tabs["GSE12667"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`Smoking_Status`,
                'race'=`Race`,
                'cancer_type' = 'Pathology',
                'metadata_sex'=Gender) %>%
                #"BMI"=bmi) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         cancer=ifelse(source_name_ch1=="Lung", "y", "n"),
         metadata_sex = case_when(
           metadata_sex == "f" ~ "female",
           metadata_sex == "m" ~ "male",
         ),
         smok=case_when(
           smok=="Former" ~ "FS",
           smok=="Current" ~ "S",
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]]) %>%
  select(-source_name_ch1, -title, -description)

gse12667 %>% mutate(across(c(smok, metadata_sex, race, cancer), as.factor))
gse12667 <- gse12667[, -c(7, 9, 10:83)]
write.csv(gse12667, 'data/pdata/gse12667.csv')
  
# GSE1650 - note: very little phenotypic detail # Ask Emily what to do about this?
gse1650 <- lung_tabs["GSE1650"][[1]] %>%
  dplyr::rename('tissue' = `source_name_ch1`)
gse1650$copd <- 'y'

gse1650 <- gse1650[, -c(3, 4, 5, 8)]
write.csv(gse1650, 'data/pdata/gse1650.csv')


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
         cancer_type="lung",
         cancer=ifelse(source_name_ch1=="Lung", "y", "n"),
         #tissue_source=case_when(
         #  tissue_source==
         #)
         smok=case_when(
           smok=="nonsmoker" ~ "NS",
           smok=="former" ~ "FS"
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]]) %>%
  select(-source_name_ch1, -title, -description)
gse32537 %>% mutate(across(c(smok, metadata_sex, cancer), as.factor),
                    age=as.numeric(age))

gse32537 <- gse32537[, -c(18, 20)]
write.csv(gse32537, 'data/pdata/gse32537.csv')

# GSE37768
gse37768 <- lung_tabs["GSE37768"][[1]] %>%
  dplyr::rename('smok'=`phenotype`) %>%
  mutate(smok=case_when(
    smok=="Nonsmoker" ~ "NS",
    smok=="healthy smoker" ~ "S",
    smok=="moderate COPD" ~ "copd"))

gse37768$copd <- 'n'
gse37768[21:38, 10] <- 'y'

gse37768 <- gse37768[, -c(3, 4)]
write.csv(gse37768, 'data/pdata/gse37768.csv')

# GSE43458
gse43458 <- lung_tabs["GSE43458"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'="smoking status") %>%
  mutate(cancer_type="adenocarcinoma",
         cancer=ifelse(str_detect(title, "adenocarcinoma"), "y", "n"),
         smok=case_when(
           smok=="Never-smoker" ~ "NS",
           smok=="Smoker" ~ "S",
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]][[2]]) %>%
  select(-source_name_ch1, -title, -description)
gse43458 %>% mutate(across(c(smok), as.factor))

gse43458 <- gse43458[, -c(5, 6)]
write.csv(gse43458, 'data/pdata/gse43458.csv')

# GSE43580
gse43580 <- lung_tabs["GSE43580"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`smoking status`,
                'age'= "age at excision (years)",
                'metadata_sex' = gender,
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
         cancer_type="lung",
         cancer=ifelse(source_name_ch1=="Lung", "y", "n"),
         #tissue_source=case_when(
         #  tissue_source==
         #)
         smok=case_when(
           smok=="Occasional Use" ~ "S",
           smok=="Previous Use" ~ "FS",
           smok=="Current Use" ~ "S"
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]][[2]]) %>%
  select(-source_name_ch1, -title, -description)
gse43580$cancer <- 'y'
gse43580 %>% mutate(across(c(smok, metadata_sex, cancer), as.factor),
                    age=as.numeric(age))

gse43580 <- gse43580[, -c(8)]
write.csv(gse43580, 'data/pdata/gse43580.csv')

# GSE63882
gse63882 <- lung_tabs["GSE63882"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`smoker`,
                'metadata_sex' = "gender",
                ) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         cancer_type="adenocarcinoma",
         metadata_sex = case_when(
           metadata_sex == "f" ~ "female",
           metadata_sex == "m" ~ "male",
         ),
         smok=case_when(
           smok=="Y" ~ "S",
           smok=="N" ~ "NS"
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]][[2]]) %>%
  select(-source_name_ch1, -title, -description)

gse63882 %>% mutate(across(c(smok, metadata_sex, cancer), as.factor),
                    age=as.numeric(age))

gse63882$cancer <- 'y'
gse63882 <- gse63882[, -c(15)]
colnames(gse63882)[11] <- c("pack_years")
colnames(gse63882)[12] <- c("cancer_type")

write.csv(gse63882, 'data/pdata/gse63882.csv')


