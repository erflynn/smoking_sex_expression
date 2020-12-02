


# LUNG
# GSE10072
# GSE103174
# GSE12667
# GSE1650
# GSE32537
# GSE37768
# GSE43458
# GSE63882

library('tidyverse')
library('GEOquery')
gse10 <- getGEO("GSE10072")
head(pData(gse10[[1]])) # NCI - Lee, adenocarcinoma study
# reanalyzed by GSE60486

pDat10 <- pData(gse10[[1]]) %>%
  select(geo_accession, title, source_name_ch1, contains("characteristics"),
         "Age at Diagnosis:ch1", "Cigarette Smoking Status:ch1", "Gender:ch1", "Stage:ch1") 

pDat10.2 <- pDat10 %>%
  group_by(geo_accession) %>%
  mutate(across(c(title, source_name_ch1),as.character)) %>%
  mutate(across(contains("characteristics_ch1"), ~clean_col(.))) %>%
  mutate(across(everything(), ~replace_na(., replace=""))) %>%
  ungroup() %>% 
  rename("sex"="characteristics_ch1", "age"="characteristics_ch1.1",
         "smok_status"="characteristics_ch1.2", "stage"="characteristics_ch1.3",
         "age2"=`Age at Diagnosis:ch1`, "smok_status2"=`Cigarette Smoking Status:ch1`, 
         "sex2"="Gender:ch1", "stage2"=`Stage:ch1`)
# all matches
pDat10.2 %>% 
  filter(sex != sex2 | age != age2 | stage != stage2 | smok_status != smok_status2) 
pDat10.3 <- pDat10.2 %>%
  select(-contains("2"), -source_name_ch1) %>%
  separate(title, into=c("cancer_status", "id"), sep="_") %>%
  mutate(sex=tolower(sex)) %>%
  mutate(smok=case_when(smok_status=="Never Smoked"~"NS", 
                        smok_status=="Current Smoker"~"S",
                        smok_status=="Former Smoker"~"FS")) %>%
  mutate(cancer=str_detect(cancer_status, "Tumor"))

# remove tumor samples???
pDat10.4 <- pDat10.3 %>% filter(!cancer, smok_status!="FS")


gse11 <- getGEO("GSE103174") # IDIBAPS Barcelona - Guillaume
pData11 <- pData(gse11[[1]]) %>% 
  select(geo_accession, title, contains("ch1")) %>% 
  select(-contains("protocol"), -contains("taxid"), -contains("label"), -contains("organism"))

gse12 <- getGEO("GSE12667")
gse13 <- getGEO("GSE1650") # no demo
gse14 <- getGEO("GSE32537")
gse15 <- getGEO("GSE37768") # no demo
gse16 <- getGEO("GSE43458") # no demo
gse17 <- getGEO("GSE63882")

