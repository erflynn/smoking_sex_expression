
# for each of the studies
#  - counts
#  - platform
#  - tissue type
#  - demographic info available


library(tidyverse)
source("code/00_utils.R")
load("data/study_sample_info.RData") # --> all_stats
incl_studies <- read_csv("data/incl_studies.csv")
my_studies <- all_stats[incl_studies$study_acc]


# ---- blood ---- #

blood_kept2 <- read_csv("data/blood_sample_labels.csv")

# in blood, 7 studies contained at least 5 smokers and 5 non-smokers (healthy)
# remove: SRP045500 (too many blood components), GSE103174 (LUNG), GSE54837 (all copd)
(blood_smok_studies <-blood_kept2 %>%
    group_by(study_acc, smok) %>%
    filter(smok %in% c("NS", "S")) %>%
    count() %>%
    pivot_wider(names_from="smok", values_from="n", values_fill=0) %>%
    filter(NS >= 5 & S >= 5) %>% 
    filter(! study_acc %in% c("SRP045500", "GSE103174", "GSE54837")))

# GSE18723 - menopause
gse18723 <- blood_kept2 %>% inner_join(
  my_studies[["GSE18723"]]$df %>%
  select(sample_acc, age, `cell type`, `menopausal`, `smoking`) %>%
  mutate(race_ethnicity="white")) %>%
  rename(tissue=`cell type`)

gse18723 %>% write_csv("data/pdata_filt/gse18723.csv")

# GSE21862
gse21862 <- blood_kept2 %>% inner_join(my_studies[["GSE21862"]]$df %>% 
  select(sample_acc, age, chip_id, gender, benzene, smoking, description, `subject id`, source_name_ch1) %>%
  rename(tissue=source_name_ch1))  
gse21862 %>% filter(gender!=sex_lab)   # none
gse21862 <- gse21862 %>% 
  select(-description, -gender, -benzene, -smoking) 
gse21862 %>% write_csv("data/pdata_filt/gse21862.csv")


# GSE87072
gse87072 <- blood_kept2 %>% 
           inner_join(my_studies[["GSE87072"]]$df) %>%
           select(study_acc, sample_acc, smok, gender, age) %>%
           mutate(tissue="pbmc") %>%
           rename(sex_lab=gender)
gse87072 %>% write_csv("data/pdata_filt/gse87072.csv")


# GSE20189 looks good, no addentional demo
blood_kept2 %>% 
  filter(study_acc =="GSE20189") %>% 
  left_join(my_studies[["GSE20189"]]$df ) %>%
  fct_summ()

gse20189 <- blood_kept2 %>% 
  filter(study_acc =="GSE20189") %>%
  mutate(tissue="whole blood")

gse20189 %>% write_csv("data/pdata_filt/gse20189.csv")

# "GSE42057" - alot 
gse42057 <- blood_kept2 %>% 
  filter(study_acc =="GSE42057")  %>% 
  left_join(my_studies[["GSE42057"]]$df ) %>%
  dplyr::select(-finalgold, -smokcignow, -source_name_ch1, -title, )
gse42057 %>% fct_summ()
gse42057 %>% 
  mutate(ref_sex=ifelse(gender==0, "female", "male")) %>%
  filter(sex_lab != ref_sex) # only 1 we could not label
gse42057 %>% write_csv("data/pdata_filt/gse42057.csv")

# "GSE55962" - has age
gse55962 <- blood_kept2 %>% 
  filter(study_acc =="GSE55962") %>% 
  left_join(my_studies[["GSE55962"]]$df) %>% 
  dplyr::select(-description, -source_name_ch1, -title, -`disease status`, -`smoking status`)
gse55962 <- gse55962 %>% filter(sex_lab == gender) %>% rename(tissue=`cell type`) 
gse55962 %>% write_csv("data/pdata_filt/gse55962.csv")

# "GSE56768" - remove stimulated
gse56768 <- blood_kept2 %>% filter(study_acc =="GSE56768") %>% 
  left_join(my_studies[["GSE56768"]]$df %>%  
              dplyr::select(sample_acc, source_name_ch1, tissue, stimulus, `cell marker`, `cell type`)) %>%
  filter(source_name_ch1=="whole blood", is.na(stimulus) | stimulus=="unstimulated") %>%
  select(-stimulus, -source_name_ch1, -`cell marker`, -`cell type`)

gse56768 %>% write_csv("data/pdata_filt/gse56768.csv")

## update counts
blood_kept_adj <- do.call(rbind, lapply(list(gse20189, gse42057, gse56768, gse18723, gse21862, gse87072, gse55962),
       function(x) x %>% select(study_acc, sample_acc, smok, sex_lab, tissue)))

blood_kept_adj %>% write_csv("data/blood_kept_phe.csv")

(blood_smok_counts <- blood_kept_adj %>% 
  filter(smok %in% c("S", "NS"), sex_lab %in% c("male", "female"), !is.na(sex_lab)) %>%
  group_by(study_acc, tissue, smok) %>% 
  count() %>%
  pivot_wider(names_from="smok", values_from="n") )

blood_sex_counts <- blood_kept_adj %>% 
  filter(!is.na(sex_lab),sex_lab != "unlabeled", !is.na(smok), smok != "FS" ) %>%
  group_by(study_acc, tissue, sex_lab, smok) %>% 
  count() %>%
  unite(grp, c(sex_lab, smok)) %>%
  pivot_wider(names_from="grp", values_from="n", values_fill=0)

# go thru the ones with sufficient samples
(suff_samples <- blood_sex_counts %>% filter(across(contains("S"), ~.>=3)))

# write out counts

# ---- lung ---- #


lung_kept <- read_csv("data/lung_sample_labels.csv")
(lung_smok <- lung_kept %>%
  group_by(study_acc, smok) %>%
  filter(smok %in% c("NS", "S")) %>%
  count() %>%
  pivot_wider(names_from="smok", values_from="n", values_fill=0) %>%
  filter(NS >= 5 & S >= 5) %>%
  filter(study_acc != "GSE40364") )

(lung_counts <- lung_kept %>% 
    semi_join(lung_smok) %>% 
  filter(!is.na(sex_lab),sex_lab != "unlabeled", smok %in% c("S", "NS") ) %>%
  group_by(study_acc, sex_lab, smok) %>% 
  count() %>%
  unite(grp, c(sex_lab, smok)) %>%
  pivot_wider(names_from="grp", values_from="n", values_fill=0))
lung_counts

(suff_samples_lung <- lung_counts %>% filter(across(contains("S"), ~.>=3)))

# GSE103174  - filter out COPD, actually LUNG
gse103174 <- blood_kept2 %>% 
  inner_join(
    my_studies[["GSE103174"]]$df %>% 
      filter(copd!="yes") %>%
      select(-copd, -source_name_ch1,-title,-treatment_protocol_ch1,-group,-smoking)) 

gse103174 <- gse103174 %>% filter(sex == sex_lab) %>% select(-sex) # 1
gse103174 %>% write_csv("data/pdata_filt/gse103174.csv")


# add age and bi(?)
gse31210 <- lung_kept %>% filter(study_acc=="GSE31210") %>% 
  left_join(my_studies[["GSE31210"]]$df %>%
              dplyr::select(sample_acc, `age (years)`, bi, gender))  %>%
  mutate(tissue="lung")
gse31210 %>% write_csv("data/pdata_filt/gse31210.csv")

# GSE32539
gse32539 <- lung_kept %>% filter(study_acc=="GSE32539") %>% 
  left_join(my_studies[["GSE32539"]]$df %>%
              dplyr::select(sample_acc, age, rin, `quit how many years ago`, `pack years`, `tissue source`, gender)) %>%
  mutate(tissue="lung")

gse32539 %>% filter(sex_lab != gender)
gse32539 <- gse32539 %>% filter(sex_lab == gender) 
gse32539 %>% write_csv("data/pdata_filt/gse32539.csv")

save(gse31210, gse32539, file="data/lung_kept.RData")


# GSE31908 
gse31908 <- lung_kept %>% 
           inner_join(my_studies[["GSE31908"]]$df %>% 
                        select(sample_acc, race, age, sex) %>%
                        mutate(tissue="lung"))
gse31908 %>% 
  mutate(sex=ifelse(sex=="f", "female", "male")) %>%
  filter(sex_lab != sex)

gse31908 %>% write_csv("data/pdata_filt/gse31908.csv")

# GSE37768 
fct_summ(lung_kept %>% inner_join(my_studies[["GSE37768"]]$df))
gse37768 <- lung_kept %>% filter(study_acc=="GSE37768") %>% mutate(tissue="peripheral lung")
gse37768 %>% write_csv("data/pdata_filt/gse37768.csv")

## ----- other tissues ! ------ #

# "GSE7895" -- bronchial
fct_summ(my_studies[["GSE7895"]]$df)
gse7895 <- my_studies[["GSE7895"]]$df %>%
  select(-description, -title) %>%
  rename(tissue=source_name_ch1) %>%
  mutate(smok=case_when(
    str_detect(characteristics_ch1, "current") ~ "S",
    str_detect(characteristics_ch1, "former") ~ "FS",
    str_detect(characteristics_ch1, "never") ~ "NS")) %>%
  separate(age, into=c("age", "pack_years", "mos_since_quit"), sep=":", fill="right") %>%
  mutate(across(c(age, pack_years, mos_since_quit), ~str_extract(., "[0-9\\.]+"))) %>%
  select(-characteristics_ch1) %>%
  add_sl() 
gse7895 %>% write_csv("data/pdata_filt/gse7895.csv")

# "GSE37147" - bronchial, all current or former smokers
# gse37147 <- my_studies[["GSE37147"]]$df %>%
#   filter(copd !="yes", `history of asthma`=="no", `inhaled medications`=="no") %>%
#   select(-description, -copd, -source_name_ch1, -title, -`used in analysis`, -`history of asthma`, -`inhaled medications`) %>%
#   add_sl() %>%
#   mutate(smok=case_when(
#     `smoking status`=="na" ~ "NS",
#     str_detect(`smoking status`, "cs") ~ "S",
#     str_detect(`smoking status`, "ex") ~ "FS"
#   )) %>%
#   mutate(ref_sex=case_when(
#     sex=="m" ~ "male",
#     sex=="f" ~ "female",
#   )) %>%
#   rename(pack_years=`pack years`) %>%
#   select(-`smoking status`,  -`sex`) %>%
#   mutate(across(everything(), ~ifelse(.=="na", NA, .)))
# 
# gse37147 %>% filter(sex_lab != ref_sex)
# gse37147 <- gse37147 %>% 
#   filter(sex_lab==ref_sex | is.na(ref_sex) | sex_lab=="unlabeled") 
# 
# gse37147 %>%
#   write_csv("data/pdata_filt/gse37147.csv")

# "GSE19027" - bronchial 
fct_summ(my_studies[["GSE19027"]]$df)
gse19027 <- my_studies[["GSE19027"]]$df %>%
  filter(status %in% c("ns", "snc")) %>%
  select(study_acc, sample_acc, age, patient_id, pkyrs, race, sex, status) %>%
  mutate(smok=ifelse(status=="ns", "NS", "S")) %>%
  select(-status) %>%
  mutate(race_ethnicity=case_when(
    race=="afa" ~ "black",
    race=="asi" ~ "asian",
    race=="his" ~ "hispanic",
    race=="cau" ~ "white",
    race=="oth" ~ "other"
  )) %>%
  select(-race) %>%
  add_sl() %>%
  group_by(study_acc, sample_acc) %>%
  summarize(across(everything(), ~paste(unique(.), collapse=";"))) %>% # condense by ID
  mutate(across(everything(),~ifelse(.=="NA", NA, .))) %>%
  ungroup() %>% mutate(tissue="bronchial epithelium")
gse19027 <- gse19027 %>% filter(sex_lab == sex | is.na(sex_lab))
gse19027 %>% distinct(patient_id, sex_lab, smok) %>% group_by(sex_lab, smok) %>% count()

gse19027 %>% write_csv("data/pdata_filt/gse19027.csv")

# "GSE8987" - buccal mucosa and nasal epithelium
gse8987 <- my_studies[["GSE8987" ]]$df %>% select(study_acc, sample_acc, source_name_ch1) %>%
  mutate(smok=ifelse(str_detect(source_name_ch1, "never"), "NS", "S"),
         tissue=ifelse(str_detect(source_name_ch1, "nasal"), "nasal epithelium", "buccal mucosa")) %>%
  select(-source_name_ch1) %>%
  add_sl()

gse8987 %>% group_by(tissue, smok) %>% count()
gse8987_1 <- gse8987 %>% filter(tissue=="buccal mucosa")
gse8987_2 <- gse8987 %>% filter(tissue=="nasal epithelium")
gse8987_1 %>% write_csv("data/pdata_filt/gse8987_1.csv")
gse8987_2 %>% write_csv("data/pdata_filt/gse8987_2.csv")

# "GSE40013" - stress, whole saliva
fct_summ(my_studies[["GSE40013" ]]$df)
gse40013 <- my_studies[["GSE40013" ]]$df %>%
  select(study_acc, sample_acc, age, `chronic stress level`, `race`, `smoking status`, `gender`, tissue) %>%
  rename(stress=`chronic stress level`, smok=`smoking status`, race_ethnicity=race) %>%
  mutate(smok=case_when(
    smok=="current" ~ "S",
    smok=="never" ~ "NS",
    smok=="former" ~ "FS"
  )) %>%
  add_sl()

gse40013 %>% filter(gender!=sex_lab & sex_lab != "unlabeled")

gse40013 <- gse40013 %>% filter(gender==sex_lab | sex_lab =="unlabeled") 
gse40013 %>% write_csv("data/pdata_filt/gse40013.csv")

# "GSE17913" - buccal mucosa, keep
gse17913 <- my_studies[["GSE17913" ]]$df %>% 
  select(study_acc, sample_acc, gender, `smoking status`, tissue) %>%
  add_sl() %>%
  mutate(smok=ifelse(`smoking status`=="never smoker", "NS", "S")) %>%
  select(-`smoking status`)

gse17913 <- gse17913 %>% filter(sex_lab == gender | sex_lab=="unlabeled")

gse17913 %>% write_csv("data/pdata_filt/gse17913.csv")

# M"GSE16149" - all female?, s + ns buccal mucosa
# keep: smoking study
fct_summ(my_studies[["GSE16149" ]]$df)
gse16149 <- my_studies[["GSE16149" ]]$df %>%
  select(study_acc, sample_acc, `smoking status`) %>%
  add_sl() %>%
  mutate(tissue="buccal mucosa") %>%
  mutate(smok=ifelse(`smoking status`=="non-smoker", "NS", "S")) %>%
  select(-`smoking status`) 
gse16149 %>% write_csv("data/pdata_filt/gse16149.csv")

#M "GSE42743"  (few S/NS female) - oral cancer
incl_studies %>% filter(study_acc=="GSE42743") %>% pull(description)
fct_summ(my_studies[["GSE42743" ]]$df)
gse42743 <- my_studies[["GSE42743" ]]$df %>%
  filter(`tu-vs-nl`=="normal") %>%
  select(study_acc, sample_acc, "age@dx", dxdate, gender, `smoking status`) %>%
  mutate(smok=case_when(
    str_detect(`smoking status`, "never") ~ "NS",
    str_detect(`smoking status`, "former") ~ "FS",
    str_detect(`smoking status`, "current") ~ "S")) %>%
  select(-`smoking status`) %>%
  add_sl() %>%
  mutate(tissue="oral")
gse42743 %>% 
  group_by(smok, sex_lab) %>% count()

gse42743 %>% write_csv("data/pdata_filt/gse42743.csv")


#? "GSE26549" - oral cancer
#  ALL come from lesions?? need to be careful about which treatment arm
# fct_summ(my_studies[["GSE26549"]]$df %>% dplyr::select(-description))
# gse26549 <- my_studies[["GSE26549"]]$df %>% 
#   filter(#outcome=="no oral cancer development" &
#     `time of biopsy`=="biopsy at baseline") %>% 
#   select(-`oral cancer-free survival time (years)`, -title, -description, 
#          -`treatment arm`, -`time of biopsy`, -treatment_protocol_ch1) %>%
#   add_sl() %>%
#   mutate(across(everything(), ~ifelse(.=="na", NA, .))) %>%
#   mutate(smok=case_when(
#     str_detect(`smoking habits`, "never") ~ "NS",
#     str_detect(`smoking habits`, "former") ~ "FS",
#     str_detect(`smoking habits`, "current") ~ "S")) %>%
#   rename(tissue="source_name_ch1") 
# 
# gse26549 <- gse26549 %>% filter(sex==sex_lab | sex_lab=="unlabeled")
# gse26549 %>% write_csv("data/pdata_filt/gse26549.csv")


# all alveolar macrophages after removing COPD
gse13896 <- my_studies[["GSE13896"]]$df %>% 
  filter(!str_detect(title, "copd")) %>%
  mutate(tissue=ifelse(str_detect(source_name_ch1, "airway"), "airway epithelium", "alveolar macrophages"),
         ref_sex=ifelse(sex=="f", "female", "male"),
         race_ethnicity=case_when(
           ancestry=="african" ~ "black",
           ancestry=="european" ~ "white",
           ancestry=="hispanic" ~ "hispanic"),
         smok=ifelse(`smoking status`=="non-smoker", "NS", "S"),
         pack_years=str_extract(`smoking status`, "[0-9\\.]+")) %>%
  add_sl() %>%
  select(-ancestry, -description, -`smoking status`, -sex, -source_name_ch1, -title)
  
gse13896 <- gse13896 %>% filter(ref_sex == sex_lab | sex_lab=="unlabeled")
gse13896 %>% group_by(smok, sex_lab) %>% count()
gse13896 %>% write_csv("data/pdata_filt/gse13896.csv")

# all AM
fct_summ(my_studies[["GSE2125"]]$df)
gse2125 <- my_studies[["GSE2125"]]$df %>%
  select(study_acc, sample_acc, status) %>%
  mutate(tissue="alveolar macrophages") %>%
  mutate(smok=ifelse(status=="nonsmoker", "NS", "S")) %>%
  add_sl() %>%
  select(-status)
gse2125 %>% group_by(smok, sex_lab) %>% count()
gse2125 %>% write_csv("data/pdata_filt/gse2125.csv")

# add: "GSE8823"  - actually alveolar macrophages from BL - matches GSE13896
#fct_summ(my_studies[["GSE8823"]]$df)


# hippocampus
fct_summ(my_studies[["GSE44456"]]$df %>% filter(phenotype=="control"))
gse44456 <- my_studies[["GSE44456"]]$df %>% 
  filter(phenotype=="control") %>%
  select(-phenotype, -source_name_ch1, -title, -case_id, -description, -cirrhosis) %>%
  add_sl() %>%
  mutate(smok=case_when(
    smoker=="no" ~ "NS",
    smoker=="ex" ~ "FS",
    smoker=="yes" ~ "S"
  )) %>%
  filter(!is.na(smok)) %>%
  mutate(tissue="hippocampus") %>%
  select( -smoker)

gse44456 %>% filter(gender!=sex_lab)
gse44456 %>% group_by(smok) %>% count()
gse44456 %>% write_csv("data/pdata_filt/gse44456.csv")

# LCL - sex labeling is a mess
# fct_summ(my_studies[["GSE36868"]]$df %>% filter(str_detect(source_name_ch1, "control")))
# gse36868 <- my_studies[["GSE36868"]]$df %>% 
#   filter(str_detect(source_name_ch1, "control")) %>%
#   select(-description, -donor, -title, -treatment_protocol_ch1, -source_name_ch1) %>%
#   add_sl() %>%
#   mutate(smok=ifelse(`smoking status`==0, "NS", "S")) %>%
#   select(-`smoking status`)
# fct_summ(gse36868)
# gse36868  %>% filter(gender != sex_lab & sex_lab !="unlabeled") 
# gse36868 %>% write_csv("data/pdata_filt/gse36868.csv")

# "kidney"
fct_summ(my_studies[["GSE46699"]]$df %>% filter(tissue =="normal"))
gse46699 <- my_studies[["GSE46699"]]$df %>%
  filter(tissue =="normal") %>% 
  select(-title, -source_name_ch1, -tissue) %>%
  add_sl() %>%
  mutate(smok=ifelse(smoking=="no", "NS", "S")) %>%
  select(-smoking) %>%
  mutate(tissue="kidney")

gse46699 %>%
  group_by(sex_lab, smok) %>%
  count()
gse46699  %>% write_csv("data/pdata_filt/gse46699.csv")

# muscle
fct_summ(my_studies[["GSE18732"]]$df %>% filter(str_detect(title, "normal")))


gse18732 <- my_studies[["GSE18732"]]$df %>% 
  select(-asa, -description, -dmverified, -dyslipidemia, -euglycemia, -fibrat,
         -glycemiagroup, -igt, -metabolic_syndrome, -sulfonyl, -title) %>%
  mutate(ref_sex=ifelse(gender==0, "female", "male"),
         smok=case_when(
           smokercurrent == 1 ~ "S",
           smokercurrent == 0 & smokerex==1 ~ "FS",
           smokercurrent == 0 & smokerex==0 ~ "NS"
         )) %>%
  add_sl()
gse18732 %>% filter(ref_sex!=sex_lab) %>% select(ref_sex, sex_lab)
gse18732 %>% group_by(sex_lab, smok) %>% count()
gse18732 %>% write_csv("data/pdata_filt/gse18732.csv")



