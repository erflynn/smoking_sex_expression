# 01b_ae_part2_extract.R
# 11/22/2020
#
# Code for extracting covariate data for airway epithelium studies.
#
# TODO:
#  - include the code for the previous AE studies?
#  - convert this process to the same as that for blood/lung
#  - reorganize data directory
#  - add missing? nasal epithelium (GSE8987), sputum macrophages (GSE46903)
#  - rescue data that does not have expr_sex labels


library('tidyverse')
library('GEOquery')

source("code/02_phenotyping/00_phenotyping_utils.R")

prev_data <- read_csv("data/sae_sl_mapped.csv") 
(prev_studies <- prev_data %>% 
    separate_rows(study, sep=";") %>% 
    distinct(study) %>% 
    pull(study))



# ------ GSE10006 ------ #
gse <- getGEO("GSE10006") # Cornell, Strulovici-Barel
pDat <- pData(gse$GSE10006_series_matrix.txt.gz)

# overlap note: 
#    GSE119087 [giant study reanalyzing a ton of microarray data, not important - used RMA on raw data]
#    GSE60486 [looked at COPD, interstudy variability in expression in normal + dz lung tissues]
pDat %>% distinct(relation) %>% pull(relation) 
pDat %>% distinct(relation.1) %>% pull(relation.1) 

# // TODO - add in COPD?





pDat1.2<- pDat %>% process_SB_pDat1() 
pDat1.4 <- pDat1.2 %>% process_SB_pDat2()

# look at counts by sex
pDat1.4 %>% 
  filter(smok %in% c("NS", "S")) %>%
  group_by(sex, smok) %>%
  count()


pDat1.4 %>% 
  filter(!smok %in% c("NS", "S")) %>% # check other data
  group_by(sex) %>% 
  count()

pDat1.4 %>% 
  group_by(source_name, sex) %>%
  count()

# add in a COPD column
# COPD are smokers
pDat1.5 <- pDat1.4 %>%
  grab_copd_info()


pDat1.5 %>%
  write_csv("data/pdata/gse10006.csv")

# ---- look at expression data ---- #
expData <- data.frame(exprs(gse$GSE10006_series_matrix.txt.gz))
# boxplot - is it appropriately transformed?
# PC space - what are the clusters?


# sex lab QC via clustering
# - convert to genes
# - plot
gpl570 <- getGEO("GPL570")
tab <- gpl570@dataTable@table %>% select("ID", `Gene Symbol`)
relevant_genes <- tab %>% filter(`Gene Symbol` %in% c("XIST", "KDM5D", "RPS4Y1"))
expData$Gene <- rownames(expData)
expData$Gene <- rownames(expData)

#--------- "GSE4302" ------ #

gse2 <- getGEO("GSE4302") # UCSF, Woodruff
pDat2 <- pData(gse2[[1]]) %>% select(geo_accession, source_name_ch1, characteristics_ch1, relation, `sample type:ch1`)
# reanalyzed by: GSE60486 

pDat2.2 <- pDat2 %>% 
  select(-relation, -characteristics_ch1) %>% 
  rename(smok=`sample type:ch1`) %>%
  filter(!str_detect(smok, "Asthmatic")) %>% 
  mutate(smok=case_when(smok=="Smoker" ~ "S",
                        smok=="Healthy control" ~ "NS")) %>%
  rename(source_name=source_name_ch1)

pDat2.2 %>% write_csv("data/pdata/gse4302.csv")

# -------- "GSE43939" ------ #

gse3 <- getGEO("GSE43939") # Cornell, Strulovici-Barel
pDat3 <- pData(gse3[[1]])  # no relation column
pDat3.1 <-pDat3%>% 
  select(geo_accession,  source_name_ch1, `smoking status:ch1`, title) %>%
  rename(smok=`smoking status:ch1`) %>%
  mutate(dgm_id=str_extract(title, "[0-9]+")) %>%
  select(-title) %>%
  mutate(smok=case_when(smok=="smoker"~ "S",
                        smok=="nonsmoker"~"NS")) %>%
  rename(source_name=source_name_ch1)
# cilia length info but ignoring
# NOTE: large -AND- small AE

pDat3.1 %>% write_csv("data/pdata/gse43949.csv")

# -------- "GSE994" ------ #

gse4 <- getGEO("GSE994") # BU, Spira
pDat4 <- pData(gse4[[1]]) 
pDat4.1 <- pDat4%>% 
  select(geo_accession,title, relation, contains("ch1")) %>%
  select(-contains("protocol"), -contains("taxid"), -contains("organism"), -contains("molecule"))
# reanalyzed by: GSE60486 -- double check


pDat4.2 <- pDat4.1 %>% 
  select(-relation) %>%
  group_by(geo_accession) %>%
  mutate(across(contains("characteristics_ch1"), ~clean_col(.))) %>%
  mutate(across(everything(), ~replace_na(., replace=""))) %>%
  rename("sex"=characteristics_ch1.1, "sex2"=`sex:ch1`,
         "age"="characteristics_ch1.3", "age2"="age:ch1",
         "race"="characteristics_ch1.2", "race2"="race:ch1", 
         "pkyrs"="characteristics_ch1.4", "pkyrs2"="pkyrs:ch1",
         "history"="characteristics_ch1.5", "history2"="history:ch1",
         "patient_id"="characteristics_ch1.7", "patient_id2"="patient_id:ch1",
         "status"="characteristics_ch1.6", "status2"="status:ch1",
         "tissue"="characteristics_ch1", "tissue2"="tissue:ch1")

# all match
pDat4.2 %>%
  filter(sex!=sex2 | age != age2 | race != race2 | pkyrs != pkyrs2 |
           history != history2 | patient_id != patient_id2 | status != status2 |
           tissue != tissue2)

pDat4.2 %>% pull(source_name_ch1)
pDat4.3 <- pDat4.2 %>% 
  ungroup() %>%
  select(-contains("2"), -source_name_ch1) %>% 
  mutate(title=tolower(as.character(title))) %>% 
  mutate(diff_formatted=str_detect(title, "sample")) %>%
  mutate(title=str_replace_all(title, "sample|\\_|without cancer", " ")) %>%
  mutate(id=str_extract(title,"[0-9]+")) %>%
  mutate(smok_status=str_replace_all(title,"[0-9]+| ", "")) %>%
  mutate(smok=case_when(smok_status=="currentsmoker"~"S",
                        smok_status=="neversmoker" ~ "NS")) %>%
  select(-title, -smok_status) %>%
  mutate(sex=tolower(sex),
         ethnic=case_when(race=="CAU" ~ "white",
                          race=="AFA"~"black",
                          race=="ASI"~"asian",
                          race=="" ~ "")) %>%
  select(-race)
# 9/75 have demographic info, 10 are diff formatted -- are these from a different study?
# //TODO: look at status, patient id, id, history

pDat4.3 %>% write_csv("data/pdata/gse994.csv")
```


```{r}
# -------- GSE16696" ------ #

gse5 <- getGEO("GSE16696") # Cornell - Strulovici-Barel

pDat5 <-pData(gse5[[1]]) 
pDat5.1 <- pDat5 %>% 
  select(geo_accession,title, source_name_ch1,relation, 
         contains("characteristics"), "Ethnic group:ch1", "Age:ch1", "Sex:ch1", "Smoking status:ch1")
# reanalyzed by: GSE119087

pDat5.2 <- pDat5.1 %>% select(-relation) %>%
  group_by(geo_accession) %>%
  mutate(across(contains("characteristics_ch1"), ~clean_col(.))) %>%
  mutate(across(everything(), ~replace_na(., replace=""))) %>%
  ungroup() %>%
  rename(age="characteristics_ch1", sex="characteristics_ch1.1", race="characteristics_ch1.2",
         smok_status="characteristics_ch1.3", race2="Ethnic group:ch1", age2="Age:ch1", 
         sex2="Sex:ch1", smok_status2="Smoking status:ch1")
pDat5.3 <- pDat5.2 %>% 
  select(-contains("2"), -source_name_ch1) %>%
  mutate(title=as.character(title)) %>%
  mutate(mas5=str_detect(title, "MAS5")) %>%
  separate(title, into=c("tissue", "smok1"), sep=",") %>%
  mutate(id=str_extract(smok1, "[0-9]+"),
         smok1=str_extract(smok1, "[A-z|\\-]+")) %>%
  separate(smok_status, into=c("smok", "pkyrs"), sep=",") %>%
  mutate(pkyrs=str_extract(pkyrs, "[0-9]+"),
         sex=case_when(sex=="M" ~ "male", sex=="F"~"female"),
         race=ifelse(race=="hispnaic","hispanic", race)) %>%
  rename(ethnic=race)

pDat5.3 %>% filter(smok != smok1) # this is messed up. let's stick w smok
pDat5.4 <- pDat5.3 %>% select(-smok1) %>%
  mutate(smok=case_when(smok=="non-smoker"~"NS",
                        smok=="smoker"~"S"))

pDat5.4 %>% write_csv("data/gse16696.csv")

gse6 <- getGEO("GSE11906") # Cornell - Strulovici-Barel
# reanalyzed by: GSE60486 , GSE119087, 
# GSE64985 - ICA paper :)

pDat6 <- pData(gse6[[1]])

pDat6 %>% check_relation_fields()

pDat6.2 <- pDat6 %>% process_SB_pDat1()
pDat6.2 %>% filter(smok0!=smok) %>% select(smok, smok0) # smok is incomplete
pDat6.2 %>% filter(age0!=age | sex0 != sex | ethnic0 != ethnic)  # none

pDat6.3 <- pDat6.2 %>% select(-smok) %>% rename(smok=smok0) %>%
  process_SB_pDat2()

pDat6.4 <- pDat6.3 %>% grab_copd_info()
summary(pDat6.4 %>% 
          mutate(across(c(source_name, sex, ethnic, smok, copd, copd_pheno), as.factor)) %>%
          mutate(across(c(age, pkyrs), as.numeric)))

pDat6.4 %>% write_csv("data/pdata/gse11906.csv")

# ------- GSE17905 ------ #
gse7 <- getGEO("GSE17905") # Cornell - Strulovici-Barel
# reanalyzed by: GSE64985 , GSE119087
pDat7 <- pData(gse7[[1]])
pDat7 %>% check_relation_fields()
pDat7.2 <- pDat7 %>% process_SB_pDat1()
pDat7.2 %>% filter(smok0!=smok) %>% select(smok, smok0) # smok is incomplete
pDat7.2 %>% filter(age0!=age | sex0 != sex | ethnic0 != ethnic)  # none

pDat7.3 <- pDat7.2 %>% select(-smok) %>% rename(smok=smok0) %>%
  process_SB_pDat2() 
# no COPD info
summary(pDat7.3 %>% 
          mutate(across(c(source_name, sex, ethnic, smok), as.factor)) %>%
          mutate(across(c(age, pkyrs), as.numeric)))

pDat7.3 %>% write_csv("data/pdata/gse17905.csv")

# ----- GSE19722 ----- #

gse8 <- getGEO("GSE19722") # Cornell - Strulovici-Barel
pDat8 <- pData(gse8[[1]])
# reanalyzed by: GSE119087
pDat8 %>% select(contains("relation"))
pDat8.2 <- pDat8 %>% process_SB_pDat1()
pDat8.2 %>% filter(smok0!=smok) %>% select(smok, smok0) # smok is incomplete
pDat8.2 %>% filter(age0!=age | sex0 != sex | ethnic0 != ethnic) 
pDat8.3 <- pDat8.2 %>% select(-smok) %>% rename(smok=smok0) %>%
  process_SB_pDat2() 
# no COPD info
summary(pDat8.3 %>% 
          mutate(across(c(source_name, sex, ethnic, smok), as.factor)) %>%
          mutate(across(c(age, pkyrs), as.numeric)))

pDat8.3 %>% write_csv("data/pdata/gse19722.csv")


# -------- GSE18385 -------- #
gse9 <- getGEO("GSE18385") # Cornell - Strulovici-Barel
pDat9 <- pData(gse9[[1]])
pDat9 %>% check_relation_fields()
pDat9.2 <- pDat9 %>% process_SB_pDat1()
pDat9.2 %>% filter(smok0!=smok) %>% select(smok, smok0) # smok is incomplete
pDat9.2 %>% filter(age0!=age | sex0 != sex | ethnic0 != ethnic)  # none

pDat9.3 <- pDat9.2 %>% select(-smok) %>% rename(smok=smok0) %>%
  process_SB_pDat2()  

# no COPD info
summary(pDat9.3 %>% 
          mutate(across(c(source_name, sex, ethnic, smok), as.factor)) %>%
          mutate(across(c(age, pkyrs), as.numeric)))

pDat9.3 %>% write_csv("data/pdata/gse18385.csv")



# --- bronchial epithelium --- #
gse19 <- getGEO("GSE7895") # Spira, BU
pDat19 <- pData(gse19[[1]]) # age
pDat19 %>% check_relation_fields() # lots of alternative to other GSMs... hmm
pDat19.1 <- pDat19 %>% 
  select(geo_accession, title, source_name_ch1, 
         contains("characteristics_ch1")) %>%
  rename(smok=characteristics_ch1,
         age=characteristics_ch1.1,
         source_name=source_name_ch1) %>%
  group_by(geo_accession) %>%
  mutate(age=clean_col(age)) %>%
  mutate(id=str_extract(title, "[0-9]+")) %>%
  mutate(smok0=str_split(title, ",")[[1]][[1]]) %>%
  ungroup() %>%
  mutate(smok=as.character(smok))
pDat19.1 %>% filter(smok!=smok0) # none
pDat19.2 <- pDat19.1 %>% 
  select(-smok0, -title) %>%
  mutate(smok=case_when(
    smok=="Never Smoker" ~ "NS",
    smok=="Current Smoker" ~ "S",
    smok=="Former Smoker"~ "FS",
    TRUE ~smok
  ))
pDat19.2 %>% write_csv("data/pdata/gse7895.csv")

# ---- GSE19027 ---- #

gse20 <- getGEO("GSE19027") # Spira - BU
pDat20 <- pData(gse20[[1]]) # many fields!
pDat20 %>% check_relation_fields()
# lots of alternative to other GSMs... hmm
# GSE60486
pDat20.1 <- pDat20 %>% 
  select(geo_accession, title, source_name_ch1, contains("characteristics_ch1"), contains(":ch1")) %>%
  group_by(geo_accession) %>%
  mutate(across(contains("characteristics_ch1"), clean_col)) %>%
  ungroup() %>%
  rename(source_name=source_name_ch1,
         tissue=characteristics_ch1,
         sex=characteristics_ch1.1,
         ethnic=characteristics_ch1.2,
         age=characteristics_ch1.3,
         pkyrs=characteristics_ch1.4,
         history=characteristics_ch1.5,
         status=characteristics_ch1.6,
         patient_id=characteristics_ch1.7)


pDat20.1 %>%
  filter(age!=`age:ch1` | sex!=`sex:ch1` | ethnic!=`race:ch1`|
           pkyrs!=`pkyrs:ch1`| history!=`history:ch1` | status!=`status:ch1` |
           patient_id!=`patient_id:ch1`| tissue!=`tissue:ch1`) # none conflict
# do I need to keep checking?

pDat20.2 <- pDat20.1 %>% select(-contains(":ch1")) %>%
  mutate(sex=tolower(sex),
         ethnic=case_when(
           ethnic=="ASI" ~ "asian",
           ethnic=="AFA"~"black",
           ethnic=="CAU" ~ "white",
           ethnic=="HIS"~"hispanic",
           ethnic=="OTH"~"other"),
         id = str_extract(title, "[0-9]+"),
         title =str_replace_all(title, "Sample\\_[0-9]+\\_", "")
  ) 
unique(pDat20.2$history)
pDat20.2 %>% 
  filter(title != source_name | title != history) 

pDat20.2 %>% 
  select(-source_name, -title) %>%
  filter(status=="NS" & history!="never smoker" |
           status=="SNC" & !str_detect(history, "without cancer") |
           status=="SC" & !str_detect(history, "with cancer")) # none conflict!

pDat20.3 <- pDat20.2 %>%
  mutate(smok=case_when(
    history=="never smoker" ~ "NS",
    str_detect(history, "current") ~ "S",
    str_detect(history, "former") ~ "FS"
  ),
  cancer=case_when(
    status=="SC" ~ "yes",
    status=="SNC" ~ "no",
    status=="NS" ~ "no")) %>% # TODO ok assumption?
  select(-history, -status, -source_name, -title) %>%
  mutate(across(c(age, pkyrs), as.numeric))

pDat20.3 %>% write_csv("data/pdata/gse19027.csv")

# paired data - complementary analysis?
# baseline (at least 2d no smok)
# acute smoke exposure (ASE)  = 24h post 3 cigarettes
gse18 <- getGEO("GSE97010")
pDat18 <- pData(gse18[[1]]) # BU, Adam Gower, GPL17244
pDat18.1 <- pDat18 %>% 
  select(geo_accession, title,source_name_ch1, contains("characteristics_ch1"), contains(":ch1")) %>%
  group_by(geo_accession) %>%
  mutate(across(contains("characteristics_ch1"), clean_col)) %>%
  ungroup() %>%
  rename(smok_exposure=characteristics_ch1, 
         rin=characteristics_ch1.1,
         sex=characteristics_ch1.2) 

pDat18.2 <- pDat18.1 %>%
  mutate(id=str_extract(title, "[0-9]+")) %>%
  rename(source_name=source_name_ch1) %>%
  select(-title)

pDat18.2 %>% filter(rin!=`rin:ch1`, sex!=`Sex:ch1`, smok_exposure!=`smoke exposure:ch1`) # none

pDat18.3 <- pDat18.2 %>% select(-contains(":ch1")) 

pDat18.3 %>% write_csv("data/pdata/gse97010.csv")


# --- alveolar macrophages, bronchial alveolar lavage
gse21 <- getGEO("GSE13896") # Cornell,  S-B
pDat21 <- pData(gse21[[1]])
pDat21 %>% select(contains("relation")) #  GSE46903, GSE119087
pDat21.1 <- pDat21 %>% select(geo_accession, title, source_name_ch1, contains("characteristics_ch1")) %>%
  group_by(geo_accession) %>%
  mutate(across(contains("characteristics_ch1"), clean_col)) %>%
  rename(source_name=source_name_ch1,
         age=characteristics_ch1,
         sex=characteristics_ch1.1,
         ancestry=characteristics_ch1.2,
         smok=characteristics_ch1.3
  ) %>%
  mutate(pkyrs=str_extract(smok, "[0-9]+")) %>%
  grab_copd_info() %>%
  mutate(smok=case_when(
    str_detect(smok, "non-smoker") ~ "NS",
    str_detect(smok, "smoker") ~ "S",
    TRUE ~ smok)) %>%
  mutate(id=str_extract(title, "[0-9]+")) %>%
  mutate(tissue=tolower(str_split(title, ",")[[1]][[1]])) %>%
  select(-title) %>%
  ungroup() %>%
  mutate(ethnic=case_when(
    ancestry=="African" ~ "black",
    ancestry=="European" ~ "white",
    ancestry=="Hispanic" ~ "hispanic",
    TRUE ~ ancestry))
pDat21.1 %>% write_csv("data/pdata/gse13896.csv")

gse22 <- getGEO("GSE13931") # Cornell,  S-B

pDat22 <- pData(gse22[[1]]) 
pDat22 %>% select(contains("relation")) # GSE60486, GSE46903, GSE119087

pDat22.1 <- pDat22 %>% select(geo_accession, title, source_name_ch1, contains("characteristics_ch1")) %>%
  group_by(geo_accession) %>%
  mutate(across(contains("characteristics_ch1"), clean_col)) %>%
  rename(source_name=source_name_ch1,
         age=characteristics_ch1,
         sex=characteristics_ch1.1,
         ethnic=characteristics_ch1.2,
         smok=characteristics_ch1.3
  ) %>%
  mutate(rma_mas=ifelse(str_detect(title, "RMA and MAS"), "y", "n")) %>%
  mutate(tissue=tolower(str_split(title, ",")[[1]][[1]])) %>%
  mutate(pkyrs=str_extract(smok, "[0-9]+")) %>%
  mutate(smok=case_when(
    str_detect(smok, "non-smoker") ~ "NS",
    str_detect(smok, "smoker") ~ "S",
    TRUE ~ smok)) %>%
  mutate(id=str_extract(title, "[0-9]+")) %>%
  select(-title) %>%
  mutate(ethnic=case_when(
    ethnic=="African" ~ "black",
    ethnic=="European" ~ "white",
    ethnic=="Hispanic" ~ "hispanic",
    TRUE ~ ethnic))

pDat22.1 %>% write_csv("data/pdata/gse13931.csv")



# ---- GSE37147 ---- #
gse23 <- getGEO("GSE37147") #  Steiling, Spira - BU, contains sex, age, resp metrics
pDat23 <- pData(gse23[[1]])
colnames(pDat23)
#pDat23 %>% select(contains("relation")) # NONE
pDat23.1 <- pDat23 %>% 
  select(geo_accession, title, source_name_ch1, 
         contains("characteristics_ch1")) %>% # did not check ":ch1" is the same
  group_by(geo_accession) %>%
  mutate(across(contains("characteristics_ch1"), clean_col)) %>%
  rename(id=characteristics_ch1,
         used_in_analysis=characteristics_ch1.1,
         fev1=characteristics_ch1.2,
         fev1_fvc=characteristics_ch1.3,
         copd=characteristics_ch1.4,
         age=characteristics_ch1.5,
         smok=characteristics_ch1.6,
         sex=characteristics_ch1.7,
         pkyrs=characteristics_ch1.8,
         hx_asthma=characteristics_ch1.9,
         inhaled_meds=characteristics_ch1.10,
         tissue=characteristics_ch1.11,
         source_name=source_name_ch1) %>%
  ungroup()

head(pDat23.1)
pDat23.2 <- pDat23.1 %>% 
  mutate(title=as.character(title),
         source_name=as.character(source_name),
         across(c(fev1, fev1_fvc, age, pkyrs), as.numeric)) %>%
  mutate(copd0=case_when(
    str_detect(source_name, "no COPD") ~ "n",
    str_detect(source_name, "COPD") ~ "y",
    TRUE ~ ""
  ),
  smok0=case_when(
    str_detect(source_name, "ex-smoker") ~ "FS",
    str_detect(source_name, "current smoker") ~ "S",
    TRUE ~ ""
  )) %>%
  rename(title_id=title)
pDat23.2 %>% group_by(copd, copd0) %>% count()
pDat23.2 %>% group_by(smok, smok0) %>% count()
unique(pDat23.2$tissue)
# THESE MATCH!
# todo - can we assume that  NA are NS?
pDat23.3 <- pDat23.2 %>% 
  select(-smok, -copd, -source_name) %>%
  rename(smok=smok0, copd=copd0)

pDat23.3 %>% write_csv("data/pdata/gse37147.csv")


# ------ GSE10038 ------ #
# ... all COPD ... #
gse24 <- getGEO("GSE10038") # Cornell: S-B
pDat24 <- pData(gse24[[1]])
#pDat24 %>% select(contains("relation")) # NONE
pDat24.1 <- pDat24 %>% 
  select(geo_accession, title, source_name_ch1,contains(":ch1")) %>%
  rename(age=`Age:ch1`, ethnic=`Ethnic group:ch1`, sex=`Sex:ch1`,
         source_name=source_name_ch1, smok=`Smoking status:ch1`) %>%
  mutate(pkyrs=str_extract(smok, "[0-9]+")) 

# ... all COPD, skip for now ... #

# ------ GSE8823 ------ #

gse25 <- getGEO("GSE8823")  # Cornell: S-B
pDat25 <- pData(gse25[[1]])
pDat25 %>% select(contains("relation"))
# GSE46903 , GSE119087

pDat25.1 <- pDat25 %>% select(geo_accession, title, source_name_ch1,
                              contains("characteristics_ch1")) %>%
  group_by(geo_accession) %>%
  mutate(across(contains("characteristics_ch1"), clean_col)) %>%
  mutate(tissue=str_split(title, ",")[[1]][[1]]) %>%
  rename(source_name=source_name_ch1,
         age=characteristics_ch1,
         sex=characteristics_ch1.1,
         ancestry=characteristics_ch1.2,
         smok=characteristics_ch1.3) %>%
  mutate(pkyrs=str_extract(smok, "[0-9]+"),
         id=str_extract(title, "[0-9]+"),
         smok=case_when(
           str_detect(smok,"non-smoker")~"NS",
           str_detect(smok,"smoker") ~ "S")
  ) %>%
  select(-source_name, -title) %>%
  mutate(ethnic=case_when(
    ancestry=="African" ~ "black",
    ancestry=="European" ~ "white",
    ancestry=="Hispanic" ~ "hispanic",
    TRUE ~ ancestry))

pDat25.1 %>% write_csv("data/pdata/gse8823.csv")