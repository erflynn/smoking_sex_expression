library('tidyverse')

### --- STEP TWO -- ###
# 47 yes, 69 maybe, 40 no
post_annot <- read_csv("data/manual_annot_extra_050221_v2.csv")
post_annot %>% filter(str_detect(n, "tobacco heating system"))  # 43 of the maybe
maybe_ds <- post_annot %>% filter(!str_detect(n, "tobacco heating system") & keep=="maybe")  # 20

kept_ds <- post_annot %>% filter(keep=="yes") 
kept_ds %>% filter(str_detect(study_acc, "GSE")) # 41
maybe_ds %>% filter(str_detect(study_acc, "GSE")) # 18

ae_ds <- kept_ds %>% filter(!str_detect(study_acc, "GSE")) %>%
  bind_rows(maybe_ds %>% filter(!str_detect(study_acc, "GSE")))

# there are 8 ArrayExpress datasets

# are any of these in Yang et al.?

# look at what these are
filt_gsm2 <- filt_gsm_long %>%
  mutate(value=case_when(
    !is.na(value1) ~ value1,
    TRUE ~ key
  ), key=case_when(
    is.na(value1) ~ name,
    TRUE ~ key
  )
  ) %>%
  select(gse, gsm, key, value)


filt_gsm3 <- filt_gsm2 %>% filter(gse %in% c(kept_ds$study_acc, maybe_ds$study_acc))
all_ae <- read_csv("data/ae_full_gsm.csv") 
filt_gsm3 %>% filter(gsm %in% all_ae$gsm) %>% distinct(gse) # GSE20250
filt_gsm4 <- filt_gsm3 %>% filter(!gsm %in% pDat5.1$geo_accession )
split_by_gse <- filt_gsm4 %>% group_split(gse)
names(split_by_gse) <- unique(filt_gsm4$gse)
study_dfs <- lapply(split_by_gse, function(x)
  x %>% pivot_wider(names_from=key, values_from=value))
names(study_dfs) <- names(split_by_gse)
save(study_dfs, file="data/tmp_study_dfs.RData")


# GSE101353 - looks like AE studies
# large AE + basal cells
p1 <- study_dfs[[1]] %>%
  filter(`time of serial bronchoscopy (month)`=="m0" | is.na(`time of serial bronchoscopy (month)`)) %>%
  filter(`smoking status`!="copd smoker") %>%
  mutate(smok=ifelse(`smoking status`=="nonsmoker", "NS", "S")) %>%
  mutate(tissue=ifelse(str_detect(`cell type`, "basal"), "basal cells",
                       "large airway epithelium")) %>%
  select(-`smoking status`, -`cell type`, -title, -source_name_ch1, 
         -`time of serial bronchoscopy (month)`, -`cell type`) %>%
  rename(dgm_id=description)
p1 %>% write_csv("data/pdata_filt2/gse101353.csv")
#intersect(all_ae$`department of genetic medicine id`, study_dfs[[1]]$description) # only 1


# GSE102556 - brain, KEEP
p2 <- study_dfs[[2]] %>% 
  filter(phenotype=="ctrl", smoking!="na") %>%
  filter(alcool!="yes", drugs!="yes") %>%
  select(gse, gsm, source_name_ch1, tissue, gender, smoking)  %>%
  mutate(smok=ifelse(smoking=="no", "NS", "S")) %>%
  select(-source_name_ch1, -smoking) 
p2 %>% write_csv("data/pdata_filt2/gse102556.csv")
p2 %>% group_by(tissue, smok) %>% count() # 7 smokers, 8 NS
# 45 no, 41 yes(heavy) 
#  22 controls (13 male; 9 females)
# 6 different brain regions

# GSE108712 - possibly keep tbd on drinkers
study_dfs[[3]] %>% 
  filter(source_name_ch1=="normal", smoker %in% c("current", "never")) %>% 
  distinct(`patient id`, smoker, drinker) %>% 
  group_by(smoker, drinker) %>%
  count() # ... hmm drinkers

# GSE109419 - ACL tears, 5 smoker + 23 non-smokers
study_dfs[[4]] %>% fct_summ()


# GSE110907 - lung adenocarcinoma
# 41 NS, 7 S, all female
# but is it microrna only?
p5 <- study_dfs[[5]] %>% filter(str_detect(source_name_ch1, "normal")) %>%
  select(gse, gsm, gender, age, `tumor stage`, smoking_status) %>%
  mutate(smok=ifelse(smoking_status=="never smoker", "NS", "S")) %>%
  select(-smoking_status)
p5 %>% write_csv("data/pdata_filt2/gse110907.csv")

# GSE111819 - diverticulitis?
# 12 NS, 6 S
p6 <- study_dfs[[6]] %>% 
  filter(`smoking status` %in% c("current smoker", "never smoked")) %>%
  select(-title, -source_name_ch1, -description) %>%
  mutate(smok=ifelse(`smoking status`=="never smoked", "NS", "S")) %>%
  select(-`smoking status`)
p6 %>% group_by(smok, sex)%>% count()
p6 %>% write_csv("data/pdata_filt2/gse111819.csv")

#  GSE112026 - no, 4 S, 19 NS -- mostly HPV
#study_dfs[[7]] %>% filter(`sample type`=="normal uppp") %>%
#  fct_summ()

# GSE11223 - removed UC + inflamed
# --> no smoking data :/ 
#unique(study_dfs[[8]]$source_name_ch1)
#study_dfs[[8]] %>% 
#  filter(str_detect(source_name_ch1, "normal uninflamed")) %>%
#  group_by(anatomic_location)

# GSE117973 - all head and neck
#study_dfs[[9]] 

# GSE123352 - maybe -- normal lung, ever vs never (no current info)
#study_dfs[[10]] %>% fct_summ()

# GSE12428 -- all current or ex
#study_dfs[[11]] %>% head()


# GSE131391 - keep! bronchial brushings, 6 ns, 6 s
p12 <- study_dfs[[12]] %>% select(gse, gsm, source_name_ch1, `smoking status`, sex, age, `pack years`) %>%
  rename(tissue=source_name_ch1) %>%
  mutate(smok=ifelse(`smoking status`=="never smoker", "NS", "S"))
p12 %>% write_csv("data/pdata_filt2/gse131391.csv")

#  GSE132607 - maybe -- all ipf, but blood
# baseline time pt --> 25 ns, 49 S
#study_dfs[[13]] %>% 
#  filter(`follow-up time`=="baseline", !is.na(smoker)) %>% 
#  group_by(smoker, gender) %>%
#  count()

# GSE134174 - make sure not duplicated!
# trachea - 6 NS, 8 S
p14 <- filt_gsm2 %>%
  filter(gse=="GSE134174", key!="description") %>%
  pivot_wider(names_from="key", values_from="value") %>%
  filter(`smoking status` %in% c("heavy smoker", "light smoker", "nonsmoker")) %>%
  mutate(smok=ifelse(`smoking status`=="nonsmoker", "NS", "S")) %>%
  rename(tissue=source_name_ch1) %>%
  select(-`cell type`, -title) 
p14 %>% write_csv("data/pdata_filt2/gse134174.csv")

#  GSE134692 - lung
#    5 smoker, 18 NS
# removed ipf, ali
p15 <- study_dfs[[15]] %>% 
  filter(diseasestatus=="normal",
         tobacco!="former smoker") %>%
  select(gse, gsm, tissue, batch, sample.property, age, race, gender, tobacco, height, weight) %>%
  mutate(smok=ifelse(tobacco=="active smoker", "S", "NS")) %>%
  select(-tobacco)
p15  %>% write_csv("data/pdata_filt2/gse134692.csv")

# GSE135055 - no smoking info for healthy donors
#study_dfs[[16]]$gse[[1]]
#filt_gsm2 %>%
#  filter(gse=="GSE135055", key!="description") %>%
#  pivot_wider(names_from="key", values_from="value") %>%
#  filter(source_name_ch1=="healthy donor") %>%
#  fct_summ()

# GSE135304 - whole blood from pt with benign, metastatic, and no pulmonary nodules
# possibly... 
# study_dfs[[17]] %>% 
#   filter(smoking %in% c("never", "smoker")) %>% 
#   distinct(patient, smoking, `nodule class`) %>%
#   group_by(smoking, `nodule class`) %>%
#   count()

# GSE136262 - oral mucosa
#  17 S, 21 NS
p18 <- study_dfs[[18]] %>%
  mutate(smok=ifelse(`smoking status`=="never smoker", "NS", "S")) %>%
  select(gse, gsm, smok, source_name_ch1) %>%
  rename(tissue=source_name_ch1)
p18 %>% write_csv("data/pdata_filt2/gse136262.csv")

# "GSE14633" - bronchial epithelium
# 11 current, 11 never
study_dfs[[19]]$gse[[1]]                           
p19 <- filt_gsm2 %>%
  filter(gse=="GSE14633", key!="description") %>%
  group_by(gse, gsm, key) %>%
  summarize(value=paste(value, collapse=";")) %>%
  pivot_wider(names_from="key", values_from="value") %>%
  mutate(smok=ifelse(str_detect(characteristics_ch1, "current"), "S", "NS"),
         pkyrs=ifelse(is.na(pkyrs), "5", pkyrs)) %>%
  select(gse, gsm, smok, pkyrs, race, sex) %>%
  mutate(tissue="bronchial epithelium")
p19 %>% write_csv("data/pdata_filt2/gse14633.csv")


# GSE150910 -lung, ever vs never
# what is an explant?
#study_dfs[[20]] %>% 
#  filter(diagnosis=="control", `ever_smoked` %in% c("yes", "no")) %>%
#  fct_summ()

# GSE152073 - blood, all female
# unclear what to filter out? (diabetes, stroke, etc)
# if filter only 3 smokers, 24 never
study_dfs[[21]]$gse[[1]]
p21 <- filt_gsm2 %>%
  filter(gse=="GSE152073", key!="description") %>%
  group_by(gse, gsm, key) %>%
  summarize(value=paste(value, collapse=";")) %>%
  pivot_wider(names_from="key", values_from="value") %>% 
  filter(smoking %in% c("never", "current smoker"))
p21 %>% filter(`stroke (personal background)`==0, 
               `diabetes mellitus (personal background)`==0,
               `metabolic syndrome`==0, 
               `acute myocardial infarction`==0)  %>%
  select(smoking) %>% fct_summ()


#  GSE15289 - blood, 211 NS, 74 S
#  peri vs postmenopausal, HMT
p22 <- study_dfs[[22]] %>%
  filter(!is.na(`smoking status`)) %>%
  mutate(smok=ifelse(`smoking status`=="no", "NS", "S")) %>%
  select(-title, -description, -replicate, -`smoking status`)
p22 %>% write_csv("data/pdata_filt2/gse15289.csv")

# GSE155588 - AE
# 11 S, 8 NS
p23 <- study_dfs[[23]] %>%
  filter(phenotype %in% c("non smoker", "smoker")) %>%
  mutate(smok=ifelse(phenotype=="smoker", "S", "NS")) %>%
  select(gse, gsm, source_name_ch1, smok) %>%
  rename(tissue=source_name_ch1)
p23 %>% write_csv("data/pdata_filt2/gse155588.csv")

# GSE16008 - BE + Nasal 13 NS, 21-26 S
# TODO check overlap
p24 <- study_dfs[[24]] %>% 
  mutate(smok=ifelse(smoker=="never smoker", "NS", "S"),
         tissue=ifelse(str_detect(source_name_ch1, "nasal"), "nasal epithelium",
                       "bronchial epithelium")) %>%
  select(gse, gsm, sex, race, age, pkyrs, smok, tissue)
p24 %>% write_csv("data/pdata_filt2/gse16008.csv")

# GSE19738 - whole blood, 14 S, 11 NS
# not looking at LPS stimulated
study_dfs[[25]]$gse[[1]]
p25 <- filt_gsm2 %>%
  filter(gse=="GSE19738", key!="description") %>%
  group_by(gse, gsm, key) %>%
  summarize(value=paste(value, collapse=";")) %>%
  pivot_wider(names_from="key", values_from="value") %>%
  filter(disease=="healthy") %>%
  filter(`smoking status`!="quit smoking", treatment=="none") %>%
  select(gse, gsm, age, gender, `smoking status`, `treatment`, `tissue`) %>%
  mutate(smok=ifelse(`smoking status`=="smoking", "S", "NS")) %>%
  select(-treatment, -`smoking status`)
p25 %>% write_csv("data/pdata_filt2/gse19738.csv")


# GSE20250 - alveolar macrophages
# TODO: check not duplicated
p26 <- study_dfs[[26]] %>% 
  select(-description) %>%
  separate(`smoking status`, into=c("smok", "pack_years"), sep=",") %>%
  mutate(smok=ifelse(smok=="smoker", "S", "NS")) %>%
  filter(str_detect(source_name_ch1, "alveolar macrophages")) %>%
  mutate(pack_years=str_extract(pack_years, "[0-9]+")) %>%
  mutate(race_ethnicity=case_when(
    is.na(ancestry) ~ `ethnic group`,
    ancestry == "african" ~ "black",
    ancestry=="european" ~ "white",
    TRUE ~ ancestry
  )) %>%
  mutate(tissue="alveolar macrophages") %>%
  select(gse, gsm, tissue, age, sex, race_ethnicity, smok, pack_years)
p26 %>% write_csv("data/pdata_filt2/gse20250.csv")


# GSE20681 - whole blood
#  48 NS, 14 S (only 3 female S)
p27 <- study_dfs[[27]] %>% 
  filter(str_detect(`disease state`, "control"),
         `smoking status` %in% c("never", "current"))  %>%
  mutate(smok=ifelse(`smoking status`=="never", "NS", "S")) %>%
  select(gse, gsm, `sex (gender)`, `date of birth`, smok) %>%
  rename(sex=`sex (gender)`)
p27%>% write_csv("data/pdata_filt2/gse20681.csv")

# GSE23323 - whole blood (22 S, 22 NS)
p28 <- study_dfs[[28]] %>% 
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  select(gse, gsm, tissue, smok)

p28%>% write_csv("data/pdata_filt2/gse23323.csv")

# "GSE23515" - blood, 12 NS, 12 S
# removed radiation data
study_dfs[[29]]$gse[[1]]
p29 <- filt_gsm2 %>%
  filter(gse=="GSE23515", key!="description") %>%
  group_by(gse, gsm, key) %>%
  summarize(value=paste(value, collapse=";")) %>%
  pivot_wider(names_from="key", values_from="value") %>%
  filter(`radiation dose`=="0gy") %>%
  select(gse, gsm, age, donor, gender, `smoking status`, tissue) %>%
  mutate(smok=ifelse(`smoking status`=="non-smoker", "NS", "S")) %>%
  select(-donor, -`smoking status`) %>%
  mutate(age=str_extract(age, "[0-9]+"))
p29%>% write_csv("data/pdata_filt2/gse23515.csv")

# GSE25326 - all tumor
#study_dfs[[30]] %>% fct_summ()


# GSE27002 - alveolar macrophages 10NS, 13S
p31 <- study_dfs[[31]] %>% 
  mutate(smok=ifelse(`sample type`=="non-smoker", "NS", "S")) %>%
  select(gse, gsm, `cell type`, smok) %>%
  rename(tissue=`cell type`)
p31 %>% write_csv("data/pdata_filt2/gse27002.csv")

# GSE30272, 166 NS, 56 S
# brain
# removed positive alcohol
p32 <- study_dfs[[32]] %>% 
  filter(`blood ethanol` !="positive" | is.na(`blood ethanol`)) %>%
  mutate(smok=case_when(
    `smoke at death`=="yes" ~ "S",
    `smoke at death` == "no" & `smoking history`=="no" ~ "NS"
  )) %>%
  filter(!is.na(smok)) %>%
  select(-treatment_protocol_ch1, -title, 
         -source_name_ch1, -contains("surrogate")) %>%
  mutate(tisuse="dorsolateral prefrontal cortex")
p32 %>% write_csv("data/pdata_filt2/gse30272.csv")


# GSE32504 - liver, 116 NS, 29 S
# should I remove alcohol or medicine? if I do --> 3 NS
p33 <- study_dfs[[33]] %>% select(-title, -description) %>%
  filter(smoking!="na") %>%
  mutate(smok=ifelse(smoking=="smoker", "S", "NS")) %>%
  select(-smoking) 
#p33 %>% 
#  filter(`alcohol intake`=="no", `presurgical medication`=="no") %>% fct_summ()
p33 %>% write_csv("data/pdata_filt2/gse32504.csv")


# GSE39366 - all H/N tumor
#study_dfs[[34]] %>% fct_summ()

# GSE39368 - all H/N tumor
#study_dfs[[35]] %>% fct_summ()

# GSE40419 - 34 NS, 40 S
# current vs smoker?
# lung adjacent normal
p36 <- study_dfs[[36]] %>% 
  filter(str_detect(source_name_ch1, "normal"),
         `smoking_status`!="na") %>%
  select(-title, -source_name_ch1, -stage) %>%
  mutate(smok=ifelse(`smoking_status`=="never smoker", "NS", "S")) 
p36 %>% write_csv("data/pdata_filt/gse40419.csv")

# GSE44133 - sperm
#study_dfs[[37]]$gse[[1]] 
#filt_gsm2 %>%
#  filter(gse=="GSE44133", key!="description") %>%
#  group_by(gse, gsm, key) %>%
#  summarize(value=paste(value, collapse=";")) %>%
#  pivot_wider(names_from="key", values_from="value")  %>%
#  fct_summ()


# GSE45847 - PBMC
# all asthma, 18 subjects with aspirin-tolerant asthma (ATA) + 24 aerd
# not enough after filter to ata
#study_dfs[[38]] %>% filter(condition=="ata") %>% fct_summ()


# GSE47415 - whole blood, 24NS 24 S (split m/f!)
p39 <- study_dfs[[39]] %>%
  mutate(smok=ifelse(other=="smoker", "S", "NS")) %>%
  select(gse, gsm, tissue, gender, age, smok)
p39 %>% write_csv("data/pdata_filt2/gse47415.csv")

# GSE47460
# only 2 current
#study_dfs[[40]]  %>%
#  filter(`disease state`=="control") %>%
#  fct_summ()

# GSE47718 - airway basal cells, 7 S, 10NS
# check there is some dgm info
study_dfs[[41]]$gse[[1]]
p41 <- filt_gsm2 %>%
  filter(gse=="GSE47718", key!="description") %>%
  group_by(gse, gsm, key) %>%
  summarize(value=paste(value, collapse=";")) %>%
  pivot_wider(names_from="key", values_from="value")  %>%
  mutate(smok=ifelse(`smoking status`=="smoker", "S", "NS")) %>%
  select(-source_name_ch1, -`smoking status`) %>%
  rename(tissue=`cell type`)
p41 %>% write_csv("data/pdata_filt2/gse47718.csv")

# GSE5056 - large airways, looks duplicated?
# 9 NS, 13 S (one brushing per lung)
p42 <- study_dfs[[42]] %>% 
  select(-source_name_ch1, -description) %>%
  mutate(tissue="airway epithelium") %>%
  separate(`smoking status`, into=c("smok", "pack_years"), sep=",") %>%
  mutate(pack_years=str_extract(pack_years, "[0-9]+")) %>%
  mutate(smok=ifelse(smok=="smoker", "S", "NS")) %>%
  mutate(id=str_extract(title, "[0-9]+")) %>%
  select(-title) %>%
  mutate(sex=ifelse(sex=="f", "female", "male"))
p42 %>% write_csv("data/pdata_filt2/gse5056.csv")

# GSE60424 - whole blood
# no smok info for ctl
#study_dfs[[43]] %>%filter(diseasestatus=="healthy control") %>%
#  fct_summ()

# GSE62182 - 27 NS, 43 S
# lung adj normal
# microrna?
study_dfs[[44]]$gse[[1]]
p44 <- filt_gsm2 %>%
  filter(gse=="GSE62182", key!="description") %>%
  group_by(gse, gsm, key) %>%
  summarize(value=paste(value, collapse=";")) %>%
  pivot_wider(names_from="key", values_from="value") %>%
  filter(str_detect(source_name_ch1, "non-")) %>%
  filter(`smoking_status` %in% c("ns", "cs")) %>%
  mutate(smok=ifelse(`smoking_status`=="ns", "NS", "S")) %>%
  select(-source_name_ch1, -`tissue subtype`, -title, -disease, -`molecule subtype`,
         -tumor_stage, -`patient id`) 
p44 %>% write_csv("data/pdata_filt/gse62182.csv")

#  GSE65213 - PBMC 34 NS, 55 S
# alcohol? psi = perceived social isolation
p45 <- study_dfs[[45]] %>% 
  select(-description, -`raw data file column`, -title) %>%
  filter(smoke != "missing") %>%
  mutate(smok=ifelse(smoke==0, "NS", "S")) %>%
  rename(tissue=source_name_ch1)  %>%
  select(-smoke)
p45 %>% write_csv("data/pdata_filt/gse65213.csv")


# another psi study
# GSE65317 - PBMC
p46 <- study_dfs[[46]] %>% 
  select(-description, -`raw data file`, -title) %>%
  filter(smoke != "missing") %>%
  mutate(smok=ifelse(smoke==0, "NS", "S")) %>%
  rename(tissue=source_name_ch1)  %>%
  select(-smoke)
p46 %>% write_csv("data/pdata_filt/gse65317.csv")

# GSE66863 - all tumor
#study_dfs[[47]] %>% 
#  fct_summ()

# GSE68549 - blood, 16 S, 24 NS
p48 <- study_dfs[[48]] %>%
  filter(smokingstatus %in% c("currently","never")) %>%
  mutate(smok=ifelse(smokingstatus=="never", "NS", "S")) %>%
  select(-title, -eversmoked, -neversmoked, -treatment_protocol_ch1, -smokingstatus) %>%
  rename(tissue=source_name_ch1)

p48 %>% write_csv("data/pdata_filt/gse68549.csv")

# GSE68559
# brain - multiple parts, all male (5 NS, 5 S)
p49 <- study_dfs[[49]] %>% 
  rename(tissue=source_name_ch1) %>%
  filter(tissue != "ba46", tissue!="raphae nuclei") %>% # remove tissues without enough samples
  mutate(smok=ifelse(smoker=="smoker", "S", "NS")) %>%
  select(-title, -smoker)
p49 %>% write_csv("data/pdata_filt/gse68559.csv")

# GSE68571 - normal does not have smok status
#study_dfs[[50]] %>% filter(disease_state=="normal", `smoking (pack years)` != "na") 


# GSE68787 - all tumor except 2
# study_dfs[[51]]$gse[[1]]
# filt_gsm2 %>%
#   filter(gse=="GSE68787", key!="description") %>%
#   group_by(gse, gsm, key) %>%
#   summarize(value=paste(value, collapse=";")) %>%
#   pivot_wider(names_from="key", values_from="value")  %>%
#   fct_summ()

# GSE72094 - all tumor
#study_dfs[[52]] %>% fct_summ()

# GSE76326 - airway basal cells, 10 NS, 7 S
# duplicated
# study_dfs[[53]]$gse[[1]]
# filt_gsm2 %>%
#   filter(gse=="GSE76326", key!="description") %>%
#   group_by(gse, gsm, key) %>%
#   summarize(value=paste(value, collapse=";")) %>%
#   pivot_wider(names_from="key", values_from="value") %>%
#   fct_summ()

# GSE81089
# all current (1) or ex(2)
#p54 <- study_dfs[[54]] %>%  
#  mutate(id=str_extract(title, "[0-9]+"))
#p54 %>% 
#  filter(id %in% c(p54 %>% filter(str_detect(title, "n")) %>% pull(id))) %>%
#  group_by(id) %>%
#  filter(!is.na(smoking)) %>%
#  fct_summ()


# GSE83527 - lung, only 3 ns non-tumor
# study_dfs[[55]] %>%
#   filter(str_detect(source_name_ch1, "non-malignant"),
#          smoking_status  %in% c("cs", "ns")) %>%
#   select(-title, -source_name_ch1, -disease, -description) %>%
#   fct_summ()


# GSE85121 - e-cigarette SAE + AM 
study_dfs[[56]] %>% fct_summ()

# GSE8569 - only 6 normal
#study_dfs[[57]] %>% filter(str_detect(title, "normal"))

# GSE86064 - possible duplicate?
# basal cells + AE cells  9 smok, 55 NS
study_dfs[[58]] %>% fct_summ()
p58 <- study_dfs[[58]]  %>% 
  filter(`smoking status` != "copd smoker") %>%
  filter(str_detect(`cell type`, "airway basal")) %>%
  select(-source_name_ch1, -`cell type`) %>%
  mutate(tissue="basal cells") %>%
  mutate(smok=ifelse(`smoking status`=="nonsmoker", "NS", "S")) %>%
  select(-title, -`smoking status`) %>%
  rename(dgm_id=description) 
p58 %>% write_csv("data/pdata_filt/gse86064.csv")


# "GSE99573" - stool, 37 S + 74 NS
study_dfs[[59]]$gse[[1]]
p59 <- filt_gsm2 %>%
  filter(gse=="GSE99573", key!="description") %>%
  group_by(gse, gsm, key) %>%
  summarize(value=paste(value, collapse=";")) %>%
  pivot_wider(names_from="key", values_from="value") %>%
  filter(`disease status`=="normal") %>%
  select(-source_name_ch1, -treatment_protocol_ch1, 
         -`disease status`, -`family history`) %>%
  mutate(smok=ifelse(`smoking status`=="yes", "S", "NS")) %>%
  select(-title, -`smoking status`)
p59 %>% write_csv("data/pdata_filt/gse99573.csv")


### now arrayexpress ones?


# curl ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-6667/E-MTAB-6667.sdrf.txt > E-MTAB-6667.sdrf.txt
# curl ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-5278/E-MTAB-5278.sdrf.txt > E-MTAB-5278.sdrf.txt
# curl ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-3604/E-MTAB-3604.sdrf.txt > E-MTAB-3604.sdrf.txt
# curl ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-5279/E-MTAB-5279.sdrf.txt > E-MTAB-5279.sdrf.txt
# curl ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-6043/E-MTAB-6043.sdrf.txt > E-MTAB-6043.sdrf.txt
# curl ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/TABM/E-TABM-305/E-TABM-305.sdrf.txt > E-TABM-305.sdrf.txt
# curl ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MEXP/E-MEXP-1277/E-MEXP-1277.sdrf.txt > E-MEXP-1277.sdrf.txt
# curl ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MEXP/E-MEXP-2813/E-MEXP-2813.sdrf.txt > E-MEXP-2813.sdrf.txt


#### ALLL THE SAMPLES!!! CEDAR COHORT
# 6 immune cell types and 3 colonic biopsis from 323 individuals i
ae1 <- read_tsv("data/array_exp/E-MTAB-6667.sdrf.txt")
ae1.1 <- ae1 %>% filter(`Characteristics[clinical history]` %in% c("smoker", "non-smoker"),
                        `Characteristics[disease]` =="normal")  %>%
  select(`Source Name`, contains("Characteristics"))
colnames(ae1.1) <- str_replace_all(str_replace_all(colnames(ae1.1), "Characteristics\\[", ""), "\\]", "")
ae1.2 <- ae1.1 %>% 
  mutate(smok=ifelse(`clinical history`=="smoker", "S", "NS")) %>%
  select(-organism, -disease, -`clinical information`) %>%
  rename(tissue=`organism part`)

ae1.2 %>%
  write_csv("data/pdata_filt2/e-mtab-6667.csv")
ae1.2 %>% 
  group_by(tissue,  `cell type`, `smok`, sex) %>% 
  count()

# "E-TABM-305" - already looked at
ae2 <- read_tsv("data/array_exp/E-TABM-305.sdrf.txt")
ae2.1 <- ae2 %>% select(`Source Name`, contains("Characteristics")) 

colnames(ae2.1) <- 
  str_replace_all(str_replace_all(colnames(ae2.1), "Characteristics \\[", ""), "\\]", "")
ae2.1 %>% select(-TimeUnit, -BioSourceType, -Organism, -DiseaseState, -OrganismPart)

# "E-MTAB-5278" - blood, 60 NS + 59 S
ae3 <- read_tsv("data/array_exp/E-MTAB-5278.sdrf.txt")
ae3.1 <-  ae3 %>% select(`Source Name`, contains("Characteristics")) 
colnames(ae3.1) <- 
  str_replace_all(str_replace_all(colnames(ae3.1), "Characteristics\\[", ""), "\\]", "")
ae3.2 <- ae3.1 %>% 
  filter(disease=="normal", `clinical history` != "former smoker") %>%
  mutate(smok=ifelse(`clinical history`=="never smoker", "NS", "S")) %>%
  select(-organism, -disease, -`clinical history`) %>%
  rename(tissue=`organism part`) 

ae3.2 %>% write_csv("data/pdata_filt2/e-mtab-5278.csv")

# "E-MTAB-3604" - sputum 40 S, 45 NS
# check duplicated!
ae4 <- read_tsv("data/array_exp/E-MTAB-3604.sdrf.txt")
ae4.1 <-  ae4 %>% select(`Source Name`, contains("Characteristics")) 
colnames(ae4.1) <- 
  str_replace_all(str_replace_all(colnames(ae4.1), "Characteristics\\[", ""), "\\]", "") 
ae4.2 <- ae4.1 %>% filter(`Study Group` %in% c("CS", "NS")) %>%
  select(-organism) %>%
  mutate(smok=ifelse(`Study Group`=="CS", "S", `Study Group`)) %>%
  rename(tissue=OrganismPart)
ae4.2 %>% write_csv("data/pdata_filt2/e-mtab-3604.csv") 


# "E-MTAB-5279" - blood 29 NS, 30 S
ae5 <- read_tsv("data/array_exp/E-MTAB-5279.sdrf.txt")
ae5.1 <-  ae5 %>% select(`Source Name`, contains("Characteristics")) 
colnames(ae5.1) <- 
  str_replace_all(str_replace_all(colnames(ae5.1), "Characteristics\\[", ""), "\\]", "") 

ae5.2 <- ae5.1 %>%  
  filter(disease=="normal", `clinical history` != "former smoker") %>%
  mutate(smok=ifelse(`clinical history`=="never smoker", "NS", "S")) %>%
  select(-organism, -disease, -`clinical history`) %>%
  rename(tissue=`organism part`) 
ae5.2 %>% write_csv("data/pdata_filt2/e-mtab-5279.csv")


# "E-MEXP-1277" - PBMC (36 no, 20 yes)
ae6 <- read_tsv("data/array_exp/E-MEXP-1277.sdrf.txt")
ae6.1 <-  ae6 %>% select(`Source Name`, contains("Characteristics")) 
colnames(ae6.1) <- 
  str_replace_all(str_replace_all(colnames(ae6.1), "Characteristics \\[", ""), "\\]", "") 
ae6.2 <- ae6.1 %>% 
  mutate(smok=ifelse(ClinicalHistory=="smok_no", "NS", "S")) %>%
  mutate(age=str_extract(Age, "[0-9]+")) %>%
  mutate(tissue="PBMC") %>%
  rename(sex=Sex) %>%
  select(-BioSourceType, -DevelopmentalStage, -Individual, -Age, -ClinicalHistory,
         -Organism, -OrganismPart, -TargetedCellType)
ae6.2 %>% write_csv("data/pdata_filt2/e-mexp-1277.csv")

# "E-MTAB-6043" - lung samples, 37 NS, 12 S
# multiple GEO datasets, check on this
# not enough per any dataset... hmmm
# ae7 <- read_tsv("data/array_exp/E-MTAB-6043.sdrf.txt")
# ae7.1 <-  ae7 %>% select(`Source Name`, contains("Characteristics")) 
# colnames(ae7.1) <- 
#   str_replace_all(str_replace_all(colnames(ae7.1), "Characteristics\\[", ""), "\\]", "") 
# ae7.2 <- ae7.1 %>%
#   filter(disease=="normal") %>%
#   mutate(smok=case_when(
#     str_detect(`clinical history`, "Never") ~ "NS",
#     str_detect(`clinical history`, "Non-smoking") ~ "NS",
#     str_detect(`clinical history`, "current") ~ "S"
#   )) %>%
#   filter(!is.na(smok)) %>%
#   select(-`clinical history`,-disease, -organism, -`ethnic group ` ) 
# ae7.2 %>% group_by(`GEO dataset`, smok) %>% count()


# "E-MEXP-2813" - wood smoke
# ae8 <- read_tsv("data/array_exp/E-MEXP-2813.sdrf.txt")
# ae8.1 <-  ae8 %>% select(`Source Name`, contains("Characteristics")) 
# colnames(ae8.1) <- 
#   str_replace_all(str_replace_all(colnames(ae8.1), "Characteristics\\[", ""), "\\]", "") 
# ae8.1 %>% fct_summ()

