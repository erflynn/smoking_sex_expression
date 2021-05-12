
library(stringr)
library(tidyverse)

# phenotype processing for blood data
blood_pheno <- load("data/blood_pheno_data.rdata")

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

# GSE12585
gse12585 <- blood_tabs["GSE12585"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('smok_degree' = `title`,
                'num_cigpacks' = `source_name_ch1`,
                'metadata_sex'=gender) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         smok_degree=case_when(
           str_detect(smok_degree, "Heavy") ~ "heavy",
           str_detect(smok_degree, "Light") ~ "light"
         ),
         num_cigpacks=numextract(num_cigpacks)) %>%
  group_by(gsm) %>%
  select(-description, -tissue2)
gse12585 %>% mutate(across(c(smok_degree, metadata_sex), as.factor),
                    age=as.numeric(age))
write.csv(gse12585, 'data/pdata/gse12585.csv')

# GSE13985
gse13985 <- blood_tabs["GSE13985"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('age' = `Age`,
                'source' = `source_name_ch1`,
                'disease' = `Disease`,
                'drugs' = `Drugs`,
                'metadata_sex'=Gender) %>%
  mutate(metadata_sex=tolower(metadata_sex),
        age = numextract(age)) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -title, -Tissue)
gse13985 %>% mutate(across(c(disease, drugs, metadata_sex), as.factor),
                    age=as.numeric(age))
write.csv(gse13985, 'data/pdata/gse13985.csv')

# GSE34198
gse34198 <- blood_tabs["GSE34198"][[1]] %>%
  mutate(tissue="peripheral whole blood") %>%
  dplyr::rename('age' = `age [year]`,
                'diabetes_status' = `diabetes status`,
                'smoking_status' = `smoking status`,
                'weight_kg' = `weigth [kg]`,
                'height_cm' = `height [cm]`,
                'systolicbp_mmhg' = `sbp [mmhg]`,
                'ca_blockers' = `ca blockers`,
                'diastolicbp_mmhg' = `dbp [mmhg]`,
                'other_medication' = `other medication`,
                'technical_replicate' = `technical replicate`,
                'sentrixIDposition' = `sentrix id and position`,
                'metadata_sex'=Sex) %>%
  mutate(metadata_sex=tolower(metadata_sex)) %>%
  group_by(gsm) %>%
  select(-title, -source_name_ch1, -tissue2, -sentrixIDposition, -technical_replicate)
write.csv(gse34198, 'data/pdata/gse34198.csv')

# GSE42057 -- note on 12/17/20: still kind of messy, need to go back and improve
gse42057 <- blood_tabs["GSE42057"][[1]] %>%
  mutate(tissue="peripheral blood") %>%
  dplyr::rename('age' = `age_enroll`,
                'status' = `title`,
                'metadata_sex' = `gender`) %>%
  mutate(metadata_sex=case_when(
    metadata_sex==0 ~ "female",
    metadata_sex==1 ~ "male"),
    status=case_when(
           str_detect(status, "COPD") ~ "COPD",
           str_detect(status, "Control") ~ "control"
         )) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -source_name_ch1)
write.csv(gse42057, 'data/pdata/gse42057.csv')

# GSE45329
gse45329 <- blood_tabs["GSE45329"][[1]] %>%
  mutate(tissue="peripheral blood") %>%
  dplyr::rename('cell_type' = `cell type`) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -source_name_ch1, -cell_type)
write.csv(gse45329, 'data/pdata/gse45329.csv')

# GSE48424
gse48424 <- blood_tabs["GSE48424"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('disease_status' = `disease status`,
                'gestational_age' = `gestational age (weeks)`) %>%
  group_by(gsm) %>%
  select(-title, -source_name_ch1, -description, -tissue2)
write.csv(gse48424, 'data/pdata/gse48424.csv')

# GSE55962
gse55962 <- blood_tabs["GSE55962"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('smok' = `smoking status`,
                'disease_status' = `disease status`,
                'source' = `source_name_ch1`,
                'metadata_sex'=gender) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         source=tolower(source),
         disease_status=tolower(disease_status),
         smok=case_when(
           smok=="Non-smoker" ~ "NS",
           smok=="Smoker" ~ "S"
         )) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -title)
write.csv(gse55962, 'data/pdata/gse55962.csv')

# GSE60491
gse60491 <- blood_tabs["GSE60491"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('source' = `source_name_ch1`,
                'metadata_sex' = `male`) %>%
  mutate(metadata_sex=case_when(
         metadata_sex==0 ~ "female",
         metadata_sex==1 ~ "male"
       )) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -title, -negativeaffect, -neuroticism, -extraversion, -openness, -agreeableness, -conscientiousness)
write.csv(gse60491, 'data/pdata/gse60491.csv')

# GSE65213
gse65213 <- blood_tabs["GSE65213"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('income' = `lnincome`,
                'source' = `source_name_ch1`,
                'metadata_sex' = `female`) %>%
  mutate(metadata_sex=case_when(
       metadata_sex==0 ~ "male",
       metadata_sex==1 ~ "female"
     )) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -title)
write.csv(gse65213, 'data/pdata/gse65213.csv')

# GSE65298
gse65298 <- blood_tabs["GSE65298"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('income' = `lnincome`,
                'source' = `source_name_ch1`,
                'smok' = `smoke`,
                'metadata_sex' = `female`) %>%
  group_by(gsm) %>%
  mutate(metadata_sex=case_when(
     metadata_sex==0 ~ "male",
     metadata_sex==1 ~ "female"
   ),
    hispanic = case_when(
    hispanic==1 ~ "yes",
    hispanic==0 ~ "no"
  ),
  black = case_when(
    black==1 ~ "yes",
    black==0 ~ "no"
  ),
   smok=case_when(
     smok==0 ~ "NS",
     smok==1 ~ "S"
   )) %>%
  select(-description, -tissue2, -title)
write.csv(gse65298, 'data/pdata/gse65298.csv')

# GSE65317
gse65317 <- blood_tabs["GSE65317"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('income' = `lnincome`,
                'source' = `source_name_ch1`,
                'metadata_sex' = `female`) %>%
  mutate(metadata_sex=case_when(
    metadata_sex==0 ~ "male",
    metadata_sex==1 ~ "female"
  ),
  hispanic = case_when(
    hispanic==1 ~ "yes",
    hispanic==0 ~ "no"
  ),
  black = case_when(
    black==1 ~ "yes",
    black==0 ~ "no"
  )) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -title)
write.csv(gse65317, 'data/pdata/gse65317.csv')

# GSE65341
gse65341 <- blood_tabs["GSE65341"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('income' = `lnincome`,
                'source' = `source_name_ch1`,
                'metadata_sex' = `female`,
                'smok' = `smoke`) %>%
  mutate(metadata_sex=case_when(
  metadata_sex==0 ~ "male",
  metadata_sex==1 ~ "female"
  ),
  hispanic = case_when(
    hispanic==1 ~ "yes",
    hispanic==0 ~ "no"
  ),
  black = case_when(
    black==1 ~ "yes",
    black==0 ~ "no"
  ),
  married = case_when(
    married==1 ~ "yes",
    married==0 ~ "no"
  ),
   smok=case_when(
   smok==0 ~ "NS",
   smok==1 ~ "S"
   )) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -title)
write.csv(gse65341, 'data/pdata/gse65341.csv')

# GSE65403
gse65403 <- blood_tabs["GSE65403"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('income' = `lnincome`,
                'source' = `source_name_ch1`,
                'metadata_sex' = `female`,
                'smok' = `smoke`) %>%
  mutate(metadata_sex=case_when(
  metadata_sex==0 ~ "male",
  metadata_sex==1 ~ "female"
  ),
  hispanic = case_when(
    hispanic==1 ~ "yes",
    hispanic==0 ~ "no"
  ),
  black = case_when(
    black==1 ~ "yes",
    black==0 ~ "no"
  ),
  married = case_when(
    married==1 ~ "yes",
    married==0 ~ "no"
  ),
   smok=case_when(
     smok==0 ~ "NS",
     smok==1 ~ "S"
   )) %>%
  group_by(gsm) %>%
  select(-description, -tissue2, -title)
write.csv(gse65403, 'data/pdata/gse65403.csv')

# GSE75511
gse75511 <- blood_tabs["GSE75511"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('sample' = `sampletype`,
                'source' = `source_name_ch1`,
                'smok' = `smoke`,
                'metadata_sex' = `female`) %>%
  mutate(metadata_sex=case_when(
  metadata_sex==0 ~ "male",
  metadata_sex==1 ~ "female"
  ),
  black = case_when(
    black==1 ~ "yes",
    black==0 ~ "no"
  ),
  asian = case_when(
    asian==1 ~ "yes",
    asian==0 ~ "no"
  ),
  white = case_when(
    white==1 ~ "yes",
    white==0 ~ "no"
  ),
   smok=case_when(
     smok==0 ~ "NS",
     smok==1 ~ "S"
   )) %>%
  group_by(gsm) %>%
  select(-source, -description, -tissue2, -title)
write.csv(gse75511, 'data/pdata/gse75511.csv')

# GSE79092
gse79092 <- blood_tabs["GSE79092"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('source' = `source_name_ch1`,
                'smok' = `smoking`) %>%
  mutate(smok = case_when(
    smok==0 ~ 'NS',
    smok==1 ~ 'S'
  )) %>%
  group_by(gsm) %>%
  select(-source, -description, -tissue2, -title, -hedonic, -eudaimonic, -jobevaluation, -selfindependent, -selfinterdependent, -otherindependent, -otherinterdependent)
write.csv(gse79092, 'data/pdata/gse79092.csv')

# GSE87072
gse87072 <- blood_tabs["GSE87072"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('source' = `source_name_ch1`,
                'health_status' = `subject status`,
                'cell_type' = `cell type`,
                'metadata_sex' = `gender`) %>%
  mutate(smok = case_when(
    smok == "Cigarette_Smoker" ~ "S",
    str_detect(smok, "User") ~ "NS"
  )) %>%
  group_by(gsm) %>%
  select(-source, -description, -tissue2, -title, -cell_type)
gse87072$source <- 'PBMC'
write.csv(gse87072, 'data/pdata/gse87072.csv')

# GSE87656
gse87656 <- blood_tabs["GSE87656"][[1]] %>%
  mutate(tissue="blood") %>%
  dplyr::rename('source' = `source_name_ch1`,
                'impersonal_pronoun' = `impersonal pronoun`,
                'auxilary_verb' = `auxilary verb`,
                'black' = `africanamerican`,
                'white' = `caucasian`,
                'metadata_sex' = `male`,
                'smok' = `smoking`) %>%
  group_by(gsm) %>%
  mutate(metadata_sex=case_when(
    metadata_sex==0 ~ "female",
    metadata_sex==1 ~ "male",
  ),
  black = case_when(
    black==1 ~ "yes",
    black==0 ~ "no"
  ),
  white = case_when(
    white==1 ~ "yes",
    white==0 ~ "no"
  ),
  asian = case_when(
    asian==1 ~ "yes",
    asian==0 ~ "no"
  ),
  smok=case_when(
     smok==0 ~ "NS",
     smok==1 ~ "S"
   )) %>%
  select(-source, -description, -tissue2, -title, -i, -we, -you, -shehe, -they, -impersonal_pronoun, -article, -preposition, -auxilary_verb, -negation, -conjunction, -quantifier, -adverb, -cdi, -alone, -talking, -wordcount, -wordspersentence)
write.csv(gse87656, 'data/pdata/gse87656.csv')




