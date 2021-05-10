library(tidyverse)


annot <- read_csv("data/manual_annot_1108_v2.csv")
annot_h <- annot %>% filter(organism=="Homo sapiens", ! str_detect(study_acc, "SRP|DRP|ERP"),
                            num_samples >=10) # 681 --> 624 --> 534
sum(annot_h$num_samples)
annot_h %>% filter(str_detect(study_acc, "GSE")) %>% nrow() # 551 --> 462
annot_h %>% filter(!str_detect(study_acc, "GSE")) # 73 --> 62
sum((annot %>% filter(str_detect(study_acc, "GSE"), organism=="Homo sapiens"))$num_samples)
sum(annot_h %>% filter(!str_detect(study_acc, "GSE")) %>% pull(num_samples))
annot2 <- annot %>% filter(organism=="Homo sapiens", num_samples >=10) # 572
annot3 <- annot2 %>% filter(num_samples >=20) # 495
head(annot3)

con <- dbConnect(SQLite(), "../GEOmetadb.sqlite") # 11/8/2020

non_gse <- annot2 %>% filter(!str_detect(study_acc, "GSE"))
gse2 <- annot2 %>% filter(str_detect(study_acc, "GSE")) # 462
gse_smok_list <- paste(gse2$study_acc, collapse=c("\',\'"))

library('GEOmetadb')
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite") # 11/8/2020

gsm_dat <- dbGetQuery(con, sprintf("SELECT gsm, molecule_ch1, title, source_name_ch1,
                      characteristics_ch1, treatment_protocol_ch1, description, 
                      channel_count, organism_ch1 FROM gsm;"))

gse_gsm <- dbGetQuery(con, sprintf("SELECT gse, gsm FROM gse_gsm WHERE gse IN ('%s');", gse_smok_list))

dbDisconnect(con)

gsm_dat2 <- gsm_dat %>% filter(organism_ch1== "Homo sapiens")
filt_gsm <- gsm_dat2 %>% right_join(gse_gsm, by="gsm")  

filt_gsm2 <- filt_gsm %>% select(-channel_count, -molecule_ch1, -organism_ch1)
head(filt_gsm2)
filt_gsm_long <- filt_gsm2 %>% 
  pivot_longer(c(-gse, -gsm)) %>%
  filter(!is.na(value)) %>%
  separate_rows(value, sep=";\t") %>%
  mutate(value=tolower(value)) %>%
  separate(value, into=c("key", "value1"), sep=": ", extra="merge")
filt_gsm_smok <- filt_gsm_long %>%
  filter(name!="description") %>%
  filter(str_detect(key, "smok") | str_detect(value1, "smok")) 

f2 <- filt_gsm_smok %>%
  distinct(gse, gsm) %>%
  group_by(gse) %>%
  count() # 338
f2 %>% filter(n < 10)   # get rid of?

# how many have we looked at before?
manual_annot <- read_csv("data/manual_annot_1231_v2.csv")
head(manual_annot)

non_gse2 <- non_gse %>% filter(!study_acc %in% manual_annot$study_acc) # 67

# new studies
f2_new <- f2 %>% anti_join(manual_annot, by=c("gse"="study_acc")) # 135

filt_gsm_smok_new <- filt_gsm_smok %>% filter(gse %in% f2_new$gse)

smok_key <- filt_gsm_smok_new %>% filter( str_detect(key, "smok")) 
smok_ch1 <- smok_key %>% filter(name=="characteristics_ch1" & !is.na(value1)) # 109
smok_nch1 <- smok_key %>% anti_join(smok_ch1, by="gse") # 25

nch1_keep <- smok_nch1 %>% filter(str_detect(key, "current") | str_detect(key, "never")) %>%
  distinct(gse)

smok_value <- filt_gsm_smok_new %>% filter(!str_detect(key, "smok"), 
                                           str_detect(value1, "smok")) # 15
smok_value %>% group_by(value1) %>% count() %>% arrange(desc(n)) %>% View()

keep_value <- smok_value %>% group_by(value1) %>% count() %>% arrange(desc(n)) %>%
  filter((str_detect(value1, "non")  | str_detect(value1, "^never")) & 
           ! str_detect(value1, "alcoholic") )

val_keep <- smok_value %>% semi_join(keep_value) %>%  distinct(gse) # 9

smok_ch1 %>% group_by(key) %>% count() %>% arrange(desc(n))
smok_ch1 %>% group_by(value1) %>% count() %>% arrange(desc(n))  %>% View()

map_smok <- smok_ch1 %>%
  distinct(value1) %>%
  mutate(smok=case_when(
    str_detect(value1, "non") ~ "NS",
    str_detect(value1, "never") ~ "NS",
    value1=="0" ~ "NS",
    str_detect(value1, "current") ~ "S",
    str_detect(value1, "former") ~ "FS",
    str_detect(value1, "ex") ~ "FS",
    str_detect(value1, "past") ~ "FS",
    str_detect(value1, "smoker") ~ "S",
    str_detect(value1, "quit") ~ "FS",
    value1 == "ns" ~ "NS",
    value1 == "s" ~ "S"
  ))


smok_ch1.1 <- smok_ch1 %>% 
  filter(!str_detect(key, "maternal")) %>%
  left_join(map_smok) %>%
  mutate(smok=case_when(
    !is.na(smok) ~ smok,
    key %in% c("smoker", "smoking status", "current smoking", 
               "current smoker", "tobacco smoking", "smoking") & value1 =="yes" ~ "S",
    key %in% c("smoker", "smoking status", "current smoking", 
               "current smoker", "tobacco smoking", "smoking") & value1 =="no" ~ "LNS",
    key %in% c("never smoker", "neversmoked") & value1=="yes" ~ "NS",
    key %in% c("never smoker", "neversmoked") & value1=="no" ~ "ES",
    key %in% c("eversmoked", "ever_smoked", "smoking history", "hx smoking") & 
      value1=="yes" ~ "ES",
    key %in% c("eversmoked", "ever_smoked", "smoking history", "hx smoking") & 
      value1=="no" ~ "NS",
    key=="smoking_status" & value1=="cs" ~ "S",
    key %in% c("smoker", "smoke", "smoking", "smoking_status(0=nonsmoker,1=smoker)") &
      value1=="1" ~ "S",
    str_detect(value1, "yes") ~ "S",
    value1=="ever" | value1=="es" ~ "ES",
    value1=="fs" | value1=="q" ~ "FS",
    value1=="cs" | value1=="active" ~ "S",
    value1=="smoking" | value1=="sm" ~ "S",
    key %in% c("smoker", "currently_smoking", "smoker (y/n)") & value1=="y" ~ "S",
    key %in% c("smoker", "smoking history", "smoker (y/n)", "smoking") & value1=="n" ~ "NS",
    key %in% c("smoked", "smoking history") & value1=="y" ~ "ES"
  ))

kv_counts <- smok_ch1.1 %>% 
  filter(is.na(smok) & value1 != "") %>% 
  group_by(key, value1) %>% 
  count() %>% 
  arrange(desc(n))

counts_by_study <- smok_ch1.1 %>% group_by(gse, smok) %>% count() %>% 
  pivot_wider(names_from="smok", values_from="n",
              values_fill=0)


counts_by_study %>% filter(`NA` >= 10)  

ch1_keep <- counts_by_study %>% filter(((NS >= 5 | LNS >= 5) & S >= 5 )  | `NA` >= 5)  # 76

nch1_keep
ch1_keep
val_keep
non_gse2
gses <- c(nch1_keep$gse, ch1_keep$gse, val_keep$gse, non_gse2$study_acc)

to_examine <- annot2 %>% filter(study_acc %in% gses)
to_examine %>% arrange(desc(num_samples)) %>%
  write_csv("data/manual_annot_extra_050221.csv")

### 47 yes, 69 maybe, 40 no
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



######## NOW WHAT DOES THIS LOOK LIKE? ########
# 27 studies!!!
new_kept <- toupper(str_replace_all(list.files("data/pdata_filt2/"), ".csv", ""))

kept_r2 <- post_annot %>% filter(study_acc %in% new_kept) %>%
  select(study_acc, title, description, num_samples, tissue2, n, X13)
gse_kept <- kept_r2 %>% filter(str_detect(study_acc, "GSE")) %>% pull(study_acc)
paste(gse_kept , collapse='","') # 23 are GSEs
# GSE20250 failed
kept_r2 %>% filter(!str_detect(study_acc, "GSE")) %>% pull(study_acc) # 5


# steps:
# 0) get more metadata
#   - platform
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite") # 11/8/2020

gse_gpl <- dbGetQuery(con, sprintf("SELECT gse, gpl FROM gse_gpl WHERE gse IN ('%s');",
                        paste()))#, collapse="','")))
gse_gpl2 <- gse_gpl %>%
  filter(gpl!="GPL6804") # skip SNP array
list_gpl <- gse_gpl2 %>% distinct(gpl) %>% pull(gpl)
gpl_info <- dbGetQuery(con, sprintf("SELECT * FROM gpl WHERE gpl IN ('%s');",
                        paste(list_gpl, collapse="','")))

# GPL21290, GPL6804
array_info <- gpl_info %>% 
  select(gpl, technology, manufacturer, data_row_count,title)

dl_info <- gpl_info %>% select(gpl, title, manufacturer, bioc_package, supplementary_file, web_link) 

gse_gpl2 %>% 
  group_by(gpl) %>%
  count() %>% arrange(desc(n)) %>%
  left_join(array_info) %>%
  select(-title)

# 1) try to download
# -- FAILED b/c RNA-seq -- #
# GSE101353 - SRP111742
# [GSE111819 - SRP135671  [possibly cut, diverticulitis], on RB]
# [GSE131391 -  SRP198757 [possibly cut, single cell]]
# GSE134174 - SRP214324, available on ds
# GSE134692 -  available on the ds

# GSE102556 -  	SRP115956, skip mouse (GPL13112)  available on ds
# GSE110907 -   	SRP133217, [is this miRNA only?] (looks like RAW avail?)     
# GSE47718 -  	SRP024274 (looks like RAW avail?) , on RB, on recount2
# GSE136262 -   	SRP219430, available on ds
# GSE155588 -  	SRP275616, available on ds

list_dl <- str_replace_all(list.files("data/small_ds2/"), ".RData", "")


# --- GOOD ---- #
# 12 

to_read <- c("GSE32504", "GSE15289", "GSE20681", "GSE30272", 
             "GSE14633","GSE16008","GSE27002", "GSE23323",
             "GSE23515", "GSE47415", "GSE19738", "GSE5056" )
gse_gpl %>% filter(gse %in% to_read) %>% left_join(array_info)
my_gses <- lapply(to_read, function(x){
  load(sprintf("data/small_ds2/%s.RData", x)); return(gse)})
lapply(my_gses, function(x) length(x$originalData)) # all 1
lapply(my_gses, function(x) x$originalData[[1]]$exp_comment)
lapply(my_gses, function(x) x$originalData[[1]]$key_comment)


# --- ArrayExpress ones --- #
unique(ae1$`Array Design REF`) # A-MEXP-2072 - Illumina HumanHT-12_V4_0_R2_15002873_B
unique(ae2$`Array Design REF`) # A-MEXP-691 - Illumina Human-6 v1 Expression BeadChip
unique(ae3$`Array Design REF`) # "A-AFFY-44", GPL570
unique(ae4$`Array Design REF`) # "A-AFFY-44", GPL570
unique(ae5$`Array Design REF`) # "A-AFFY-44", GPL570
unique(ae6$`Array Design REF`) # A-AGIL-9 - Agilent Human 1A Microarray (V2)

# get list of files to download
kept_r2 %>% filter(!str_detect(study_acc, "GSE")) %>% pull(study_acc) # 5

#"E-MTAB-6667" - need to download from raw, Illumina HumanHT
#"E-MTAB-5278" - processed available, GPL570
#"E-MTAB-3604" - processed available, GPL570
#"E-MTAB-5279" - processed available, GPL570
#"E-MEXP-1277" - processed available - agilent



######### 2) SEX LABEL ############
# get XY chromosome genes
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
chr_map <- biomaRt::getBM(attributes=c("hgnc_symbol", "chromosome_name"), mart=mart)
xy_chr_map <- chr_map %>% filter(chromosome_name %in% c("X", "Y"))
save(xy_chr_map, file="ref/xy_chr_map.RData")
# try labeling w both

library('MetaIntegrator')

# ---- GSE32504, liver
my_gse <- getGEOData("GSE32504")
my_gse <- my_gse$originalData$GSE32504
save(my_gse, file="data/small_gses2/GSE32504.RData")
pdat <- read_csv("data/pdata_filt2/gse32504.csv")

map_to_title <- my_gse$pheno %>% select(title, geo_accession)
meta_ds <- read_excel("ref/gse32504_meta.xls", skip=4)
meta_ds %>% filter(!is.na(`Sample ID`)) %>% 
  rename(title=`Sample ID`) %>%
  inner_join(map_to_title) %>%
  inner_join(pdat, by=c("geo_accession"="gsm")) %>%
  select(geo_accession, title, Sex, gender, Smoking, smok) %>%
  group_by(Sex, gender, Smoking, smok) %>%
  count()
# looks like they swapped M and F

#my_gse <- my_gses[[1]][[1]][[1]]
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
# appears to be probe --> one gene
probe_gene %>% filter(str_detect(gene, ","))
probe_gene %>% filter(str_detect(gene, ";"))

xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe), plot=T)

# woot - good separation
sl_df <- tibble("gsm"=names(sl),
       "expr_sex"=unlist(sl))

# uhhh everything is wrong except 2??? wtf?
# ... but the paper says 79 females, 71 males --> keep ours
pdat2 <- pdat %>% 
  mutate(gender_fixed=ifelse(gender=="male", "female", "male")) %>%
  left_join(sl_df) %>% 
  filter(gender_fixed==expr_sex)  %>%
  rename(tissue=source_name_ch1,
         sex=expr_sex,
         metadata_sex=gender_fixed) %>%
  mutate(race_ethnicity="white") %>%
  select(-`ethnic background`, -gender) 


table(pdat$gender)
table(sl_df$expr_sex)

plot(t(my_gse$expr[probe_gene %>% 
  filter(gene %in% c("XIST", "RPS4Y1")) %>%
  pull(probe),] ), col=ifelse(sl=="male", "blue", "red"))
# do we trust the smoking data tho?
# ... check -- it checks out in other ways?

pdat2 %>% write_csv("data/pdata_filt2_v2/gse32504.csv")

#  "GSE15289" -- possibly skip
# TODO: check annot
my_gse <- getGEOData("GSE15289")
my_gse <- my_gse$originalData$GSE15289
my_gse$platform
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
# appears to be probe --> one gene
probe_gene %>% filter(str_detect(gene, ","))
probe_gene %>% filter(str_detect(gene, ";"))

xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T)

# not good separation, try toker
probe_gene %>% filter(!is.na(gene)) # only 7.6k genes... hmm


# "GSE20681"
gse <- getGEOData("GSE20681")
my_gse <- gse$originalData$GSE20681
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T)
# decent separation
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

pdat <- read_csv("data/pdata_filt2/gse20681.csv")
save(my_gse, file="data/small_gses2/GSE20681.RData")
# all match!
pdat %>% left_join(sl_df) %>% filter(sex==expr_sex)
pdat %>% 
  left_join(sl_df) %>% 
  rename(metadata_sex=sex, sex=expr_sex) %>%
  write_csv("data/pdata_filt2_v2/gse20681.csv")

# "GSE30272"
gse <- getGEOData("GSE30272")
my_gse <- gse$originalData$GSE30272
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T) # ok sep
pdat <- read_csv("data/pdata_filt2/gse30272.csv")
save(my_gse, file="data/small_gses2/GSE30272.RData")
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))
pdat2 <- pdat %>% 
  mutate(metadata_sex=case_when(
    sex=="m" ~ "male",
    sex=="f" ~ "female"
  )) %>%
  left_join(sl_df) %>%
  select(-sex) %>%
  rename(sex=expr_sex, tissue=tisuse)  %>%
  mutate(race_ethnicity=case_when(
    race=="aa" ~ "black",
    race=="as" ~ "asian",
    race=="hisp" ~ "hispanic",
    race=="cauc" ~ "white")) %>%
  select(-race) %>%
  select(gse, gsm, smok, sex, metadata_sex, age, race_ethnicity, everything())

pdat2 %>%  
  select(gsm, metadata_sex, expr_sex) %>%
  filter(metadata_sex!=expr_sex) # all match
pdat2 %>% write_csv("data/pdata_filt2_v2/gse30272.csv")

# "GSE14633"
gse <- getGEOData("GSE14633")
my_gse <- gse$originalData$GSE14633
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T)  # good sep!
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))
pdat <- read_csv("data/pdata_filt2/gse14633.csv")
save(my_gse, file="data/small_gses2/GSE14633.RData")

pdat %>% left_join(sl_df) %>% filter(sex!=expr_sex) # all match
pdat2 <- pdat %>% 
  left_join(sl_df) %>%
  rename(metadata_sex=sex,
         sex=expr_sex) 
pdat2 %>% write_csv("data/pdata_filt2_v2/gse14633.csv")


# "GSE16008" - bronchial and nasal epithelium
# GPL5175, Exon Array data
gse <- getGEOData("GSE16008")
my_gse <- gse$originalData$GSE16008
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 


pdat <- read_csv("data/pdata_filt2/gse16008.csv")
pdat %>% group_by(tissue) %>% count()
nasal_samples <- pdat %>% filter(tissue=="nasal epithelium") %>% pull(gsm)
sl_n <- massiRAcc(my_gse$expr[,nasal_samples], xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T, threshold=4)  # looks like there are batch effects?

bronchial_samples <- pdat %>% filter(tissue=="bronchial epithelium") %>% pull(gsm)
sl_b <- massiRAcc(my_gse$expr[,bronchial_samples], xy_genes %>% 
                    filter(chromosome_name=="Y") %>% 
                    pull(probe),  plot=T, threshold=4)  # good sep

sl_df <- tibble("gsm"=names(sl_n),
                "expr_sex"=unlist(sl_n)) %>%
  bind_rows(tibble("gsm"=names(sl_b),
                  "expr_sex"=unlist(sl_b)))

# counts match for BE (minus 1 white f NS)
pdat %>% left_join(sl_df) %>% 
  filter(tissue=="bronchial epithelium")  %>%
  filter(sex!=expr_sex) # 8/26
lapply(pdat %>% left_join(sl_df) %>% 
  filter(tissue=="bronchial epithelium") %>%
  group_split(smok), function(x) fct_summ(x))

pdat %>% left_join(sl_df) %>% 
    filter(tissue=="nasal epithelium")  %>%
    filter(sex!=expr_sex) # 16/37

lapply(pdat %>% left_join(sl_df) %>% 
         filter(tissue=="nasal epithelium") %>%
         group_split(smok), function(x) fct_summ(x))


# some of these are from previous studies??

save(my_gse, file="data/small_gses2/GSE16008.RData")



# "GSE27002" 
gse <- getGEOData("GSE27002")
my_gse <- gse$originalData$GSE27002
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T)  # good sep
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))
pdat <- read_csv("data/pdata_filt2/gse27002.csv")
save(my_gse, file="data/small_gses2/GSE27002.RData")
 # no sex information tho
pdat %>% left_join(sl_df) %>% rename(sex=expr_sex) %>%
  write_csv("data/pdata_filt2_v2/gse27002.csv")


# "GSE23323",
gse <- getGEOData("GSE23323")
my_gse <- gse$originalData$GSE23323
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

y_probes <- xy_genes %>% 
  filter(chromosome_name=="Y") %>% 
  pull(probe)

count_nas <- apply(my_gse$expr[y_probes,],1, function(x) sum(is.na(x)))
y_probes_present <- names(count_nas[count_nas==0])

sl <- massiRAcc(my_gse$expr, y_probes_present,  plot=T)  # good sep
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

pdat <- read_csv("data/pdata_filt2/gse23323.csv")
save(my_gse, file="data/small_gses2/GSE23323.RData")
# no sex information tho
pdat %>% left_join(sl_df) %>% rename(sex=expr_sex) %>%
  write_csv("data/pdata_filt2_v2/gse23323.csv")

# "GSE23515"
# Paul and Amundson
gse <- getGEOData("GSE23515")
my_gse <- gse$originalData$GSE23515
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 


y_probes <- xy_genes %>% 
  filter(chromosome_name=="Y") %>% 
  pull(probe)

count_nas <- apply(my_gse$expr[y_probes,],1, function(x) sum(is.na(x)))
y_probes_present <- names(count_nas[count_nas==0])


sl <- massiRAcc(my_gse$expr, y_probes_present,  plot=T)  # good sep
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

pdat <- read_csv("data/pdata_filt2/gse23515.csv")
save(my_gse, file="data/small_gses2/GSE23515.RData")
pdat %>% 
  left_join(sl_df) %>%
  filter(gender!=expr_sex) # all match
pdat %>% left_join(sl_df) %>%
  rename(metadata_sex=gender,
         sex=expr_sex) %>%
  write_csv("data/pdata_filt2_v2/gse23515.csv")

# "GSE47415" - I think this might be a duplicate
# Paul and Amundson
# gse <- getGEOData("GSE47415")
# my_gse <- gse$originalData$GSE47415
# probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
# xy_genes <- probe_gene %>% 
#   inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 
# 
# 
# y_probes <- xy_genes %>% 
#   filter(chromosome_name=="Y") %>% 
#   pull(probe)
# 
# count_nas <- apply(my_gse$expr[y_probes,],1, function(x) sum(is.na(x)))
# y_probes_present <- names(count_nas[count_nas==0])
# 
# 
# sl <- massiRAcc(my_gse$expr, y_probes_present,  plot=T)  # good sep
# sl_df <- tibble("gsm"=names(sl),
#                 "expr_sex"=unlist(sl))
# 
# pdat <- read_csv("data/pdata_filt2/gse47415.csv")
# 
# save(my_gse, file="data/small_gses2/GSE47415.RData")
# pdat %>% left_join(sl_df) %>%
#   rename(metadata_sex=gender,
#          sex=expr_sex) %>%
#   fct_summ()



# "GSE19738"
gse <- getGEOData("GSE19738")
my_gse <- gse$originalData$GSE19738
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

pdat <- read_csv("data/pdata_filt2/gse19738.csv")
pdat %>% fct_summ()
y_probes <- xy_genes %>% 
  filter(chromosome_name=="Y") %>% 
  pull(probe)

expr2 <- my_gse$expr[,pdat$gsm]

count_nas <- apply(expr2[y_probes,],1, function(x) sum(is.na(x)))
y_probes_present <- names(count_nas[count_nas==0])
y_probes_present

sl <- massiRAcc(expr2, y_probes_present,  plot=T, threshold=4)  
# not great sep


save(my_gse, file="data/small_gses2/GSE19738.RData")

sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

unlist(my_gse$pheno %>%
  filter(`treatment:ch1`=="none",
                       `sample type:ch1`=="control") %>%
  filter(!str_detect( characteristics_ch1.3, "quit")) %>%
  mutate(title=as.character(title)) %>% pull(title) )
    
  #fct_summ() # 20f, 14 m
pdat %>% fct_summ() # study says 23f, 12 m

# we have 16 and 9
pdat %>% 
  left_join(sl_df) %>%
  filter(gender!=expr_sex) # 11/25 do not match

#  write_csv("data/pdata_filt2_v2/gse19738.csv")

# "GSE5056" - note 2 per each
gse <- getGEOData("GSE5056")
my_gse <- gse$originalData$GSE5056
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 


y_probes <- xy_genes %>% 
  filter(chromosome_name=="Y") %>% 
  pull(probe)

count_nas <- apply(my_gse$expr[y_probes,],1, function(x) sum(is.na(x)))
y_probes_present <- names(count_nas[count_nas==0])

sl <- massiRAcc(my_gse$expr, y_probes_present,  plot=T, threshold=2)  
# ok sep at 2
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))
pdat <- read_csv("data/pdata_filt2/gse5056.csv")
save(my_gse, file="data/small_gses2/GSE5056.RData")
pdat %>% 
  left_join(sl_df) %>%
  filter(sex!=expr_sex) # all match
pdat2 <- pdat %>% left_join(sl_df) %>%
  rename(metadata_sex=sex,
         sex=expr_sex,
         race_ethnicity=`ethnic group`)  
pdat2 %>% fct_summ()
pdat2 %>% write_csv("data/pdata_filt2_v2/gse5056.csv")

# download the array express ones
#"E-MTAB-6667" - need to download from raw, Illumina HumanHT
res <- sapply(1:16, function(i)
       sprintf("wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6667/E-MTAB-6667.raw.%s.zip", i))

res %>% write_tsv("list_cmd.txt")
#"E-MTAB-3604" - processed available, GPL570
# added a comma at the start
exp <- read_csv("data/array_exp/E-MTAB-3604/COPD_SPUTUM_GCRMANormalizedMatrix.txt")
# check parsing fails?
pdat <- read_csv("data/pdata_filt2/e-mtab-3604.csv")
full_p <- read_tsv("data/array_exp/E-MTAB-3604.sdrf.txt") 
pdat2 <- pdat %>% left_join(full_p %>% select(`Source Name`, `Array Data File`))

exp_sm <- exp[,c("X1", pdat2$`Array Data File`)] %>% rename(probes=X1)
head(exp_sm[,1:5])
exp_mat <- as.matrix(exp_sm %>% select(-probes))
rownames(exp_mat) <- exp_sm$probes

.GEM_log_check_v1 <- function(GEM){
  
  #pick the first column that is not all NAs
  ref_p <- which(sapply(1:ncol(GEM$expr),
                        function(i)
                          !all(is.na(GEM$expr[,i]))))[1]
  
  #this range is obserbed only in log-normal case [in this part context]
  if(range(GEM$expr[,ref_p],na.rm=T)[1]<0){
    return(TRUE)
  }
  #if it's all positives all the bets are off
  else{
    #The data should be normally distributed
    #once in log scale
    #Get qq-norm objects
    obj_real <- qqnorm(GEM$expr[,ref_p],plot.it=FALSE);
    obj_exp  <- qqnorm(exp(GEM$expr[,ref_p]),plot.it=FALSE);
    obj_log  <- qqnorm(log(abs(GEM$expr[,ref_p])+0.00001),plot.it=FALSE);
    
    #remove infinite
    obj_exp$x[which(is.infinite(obj_exp$x))] <- NA;
    obj_exp$y[which(is.infinite(obj_exp$y))] <- NA;  
    
    #look at correlations
    cor_real <- cor(obj_real$x,obj_real$y,use='pairwise.complete');
    cor_exp  <- cor(obj_exp$x, obj_exp$y, use='pairwise.complete');
    cor_log  <- cor(obj_log$x ,obj_log$y ,use='pairwise.complete');
    
    #compute R^2 difference and take ratio
    log_check_score <- log2(abs(cor_real**2 - cor_log**2)/abs(cor_real**2 - cor_exp**2))
    
    #dataset is in log scale
    if(log_check_score<=0){
      return(TRUE);
    }
    #dataset is not in log scale
    else{
      return(FALSE);
    }
  }
}

.GEM_log_check_v1(my_gse) # log scale
my_gse <- list("expr"=exp_mat, "platform"="GPL570")
save(my_gse, file="data/small_gses2/E-MTAB-3604.RData")



load("ref/gpl570_probe_gene.RData") # probe_gene
xy_genes <- probe_gene %>% inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 
sl <- massiRAcc(exp_mat, xy_genes %>% filter(chromosome_name=="Y") %>%pull(probes),  plot=T)  # good sep
sl_df <- tibble("sample_acc"=names(sl),
                "expr_sex"=unlist(sl))

pdat3 <- sl_df %>% inner_join(pdat2, by=c("sample_acc"="Array Data File"))  %>%
  rename(metadata_sex=Sex,
         sex=expr_sex,
         age=Age,
         race_ethnicity=`ethnic group`,
         BMI=`body mass index`,
         pack_years=`Pack Year`) %>%
  mutate(race_ethnicity=ifelse(str_detect(race_ethnicity, "Black"), "black", "white")) %>%
  mutate(tissue=tolower(tissue))

pdat3 %>% 
  select(sample_acc, sex, metadata_sex, smok, tissue, age, race_ethnicity, pack_years, everything()) %>%
  write_csv("data/pdata_filt2_v2/e-mtab-3604.csv")
  
#  filter(expr_sex != Sex) # all match!

#"E-MTAB-5278" - processed available, GPL570
exp <- read_tsv("data/array_exp/E-MTAB-5278/H1_QASMC_Blood_Hs_MRNA_normMatrix.txt")
pdat <- read_csv("data/pdata_filt2/e-mtab-5278.csv")
full_p <- read_tsv("data/array_exp/E-MTAB-5278.sdrf.txt") 
pdat2 <- pdat %>% inner_join(full_p %>% select(`Source Name`, `Array Data File`)) %>% distinct()
length(unique(pdat2$`Source Name`)) # 118
length(unique(pdat2$`Array Data File`)) # 119
length(intersect(pdat2$`Array Data File`, colnames(exp))) # 116
setdiff(pdat2$`Array Data File`, colnames(exp)) # 3 missing
pdat2.2 <- pdat2 %>% filter(`Array Data File` %in% colnames(exp))
length(unique(pdat2.2$`Source Name`)) # 115 - one person has two samples
pdat2.2[duplicated(pdat2.2$`Source Name`),]
pdat2.2 %>% filter(`Source Name`=="DD3056")


exp_sm <- exp[,c("Hybridization Name", pdat2.2$`Array Data File`)] %>% 
  rename(gene=`Hybridization Name`) %>% filter(gene!="Reporter REF")
exp_mat <- as.matrix(apply(exp_sm %>% select(-gene), c(1,2), as.numeric))
rownames(exp_mat) <- exp_sm$gene
colnames(exp_mat) <- str_replace_all(colnames(exp_mat), "-", "\\.")
my_gse <- list("expr"=exp_mat)

my_gse <- list("expr"=exp_mat, platform="GPL570")
save(my_gse, file="data/small_gses2/E-MTAB-5278.RData")



sl <- massiRAcc(exp_mat, xy_genes %>% 
                  filter(gene %in% rownames(exp_mat)) %>%
                  filter(chromosome_name=="Y") %>% 
                  pull(gene),  plot=F)  # good sep

sl_df <- tibble("sample_acc"=names(sl),
                "expr_sex"=unlist(sl))

pdat3 <- sl_df %>% inner_join(pdat2.2 %>% 
                       mutate(`Array Data File`=str_replace_all(`Array Data File`, "-", "\\.")), 
                     by=c("sample_acc"="Array Data File"))

pdat4 <- pdat3 %>% filter(sex == expr_sex) %>%  # 4 mismatch
  rename(metadata_sex=sex,
         sex=expr_sex,
         id=`Source Name`,
         race_ethnicity=`ethnic group`) %>%
  mutate(race_ethnicity=ifelse(race_ethnicity=="African American", "black", "white")) 
pdat4 %>% write_csv("data/pdata_filt2_v2/e-mtab-5278.csv")

#"E-MTAB-5279" - processed available, GPL570
exp <- read_tsv("data/array_exp/E-MTAB-5279/H2_BLDSMK01_Blood_Hs_MRNA_normMatrix.txt")
pdat <- read_csv("data/pdata_filt2/e-mtab-5279.csv")
full_p <- read_tsv("data/array_exp/E-MTAB-5279.sdrf.txt") 

pdat2 <- pdat %>% inner_join(full_p %>% select(`Source Name`, `Array Data File`)) %>% distinct()
colnames(exp)
setdiff(pdat2$`Array Data File`, colnames(exp)) # 4 missing
pdat2.2 <- pdat2 %>% filter(`Array Data File` %in% colnames(exp))

exp_sm <- exp[,c("Hybridization Name", pdat2.2$`Array Data File`)] %>% 
  rename(gene=`Hybridization Name`) %>% filter(gene!="Reporter REF")
exp_mat <- as.matrix(apply(exp_sm %>% select(-gene), c(1,2), as.numeric))
rownames(exp_mat) <- exp_sm$gene

my_gse <- list("expr"=exp_mat, platform="GPL570")
save(my_gse, file="data/small_gses2/E-MTAB-5279.RData")



sl <- massiRAcc(exp_mat, xy_genes %>% 
                 filter(gene %in% rownames(exp_mat)) %>%
                 filter(chromosome_name=="Y") %>% 
                 pull(gene),  plot=T)  # good sep

sl_df <- tibble("sample_acc"=names(sl),
                "expr_sex"=unlist(sl))

pdat3 <- sl_df %>% inner_join(pdat2.2, 
                              by=c("sample_acc"="Array Data File"))

pdat3 %>% 
  rename(
    metadata_sex=sex,
    sex=expr_sex,
    id=`Source Name`,
    race_ethnicity=`ethnic group`
  ) %>%
  mutate(race_ethnicity="white") %>%
  write_csv("data/pdata_filt2_v2/e-mtab-5279.csv")


#"E-MEXP-1277" - processed available - agilent
load.Rdata("data/array_exp/E-MEXP-1277.eSet.r", "my_eset")
colnames(fData(my_eset)) # hgnc_symbol: "Composite.Element.Database.Entry.locus."
exp_mat <- eData(my_eset)
head(exp_mat[,1:5])
#my_eset@assayData$
# G, Gb, R, Rb

# 3) check duplicated???



# 4) update study data frame
#   - make sure separated by tissue
# 5) update count lists 

all_pdat <- lapply(list.files("data/pdata_filt2_v2/"), 
       function(x) read_csv(sprintf("data/pdata_filt2_v2/%s", x)))
studies <- toupper(str_replace_all(list.files("data/pdata_filt2_v2/"), ".csv", ""))
names(all_pdat) <- studies

lapply(all_pdat, colnames)
new_study_df <- do.call(rbind, lapply(studies , function(study) {
  x <- all_pdat[[study]]
  if (!"metadata_sex" %in% colnames(x)){
    x <- x %>% mutate(metadata_sex="")
  }
  if (!"tissue" %in% colnames(x)){
    x <- x %>% mutate(tissue="")
  }
  if ("gsm" %in% colnames(x)){
    x <- x %>% rename(sample_acc=gsm)
  }
  x %>% 
    mutate(study_acc=study) %>%
    select(study_acc, sample_acc, smok, metadata_sex, sex, tissue)
  })
)

df <- data.frame(new_study_df)
df2 <- df %>% 
  mutate(tissue=case_when(
    study_acc=="GSE20681" | tissue=="whole blood" ~ "blood",
    tissue=="bronchial epithelium" ~ "airway epithelium",
    str_detect(tissue, "cortex") ~ "brain",
    TRUE ~ tissue)) %>%
  group_by(tissue, study_acc, smok, sex)  %>%
  count() %>%
  unite(grp, c("smok", "sex")) %>%
  pivot_wider(names_from=grp, values_from=n) %>%
  mutate(NS=NS_female + NS_male,
         S=S_female+S_male,
         total=NS+S) 
df2 %>% write_csv("table1_part2.csv")

# get title, platform
post_annot <- read_csv("data/manual_annot_extra_050221_v2.csv")

df3 <- df2 %>% inner_join(post_annot %>% filter(study_acc %in% studies) %>%
  select(study_acc, title) %>%
  left_join(gse_gpl, by=c("study_acc"="gse")) %>%
  mutate(gpl=ifelse(is.na(gpl), "GPL570", gpl))) %>%
  select(tissue, study_acc, title, gpl, total, S, NS, everything())


df3 %>% write_csv("table_part2_w_gpl.csv")
df3 <- read_csv("table_part2_w_gpl.csv")
annot_part2 <- read_csv("table1_part2_v2.csv")


# load the platform data and save it
# gse_to_gpl <- df3 %>% filter(gpl !="GPL570") %>% ungroup() %>%
#   select(study_acc, gpl) %>%
#   group_by(gpl) %>%
#   sample_n(1)
# lapply(gse_to_gpl$study_acc, function(gse){
#   load(sprintf("data/small_gses2/%s.RData", gse))
#   probe_gene <- tibble("probes"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
#   save(probe_gene, file=sprintf("ref/%s_probe_gene.RData", tolower(my_gse$platform)))
# })


library('readxl')
study_annot <- read_excel("table1_study_info.xlsx") %>%
  select(-remove, -`duplicate?`, -notes)
form_list <- read_csv("data/supp_tables/s4_list_smok_studies_formatted.csv")
head(form_list)
form_list %>% distinct(tissue)
form_list %>%
  filter(tissue %in% c("whole blood", "pbmc")) 
form_list %>%
  filter(tissue %in% c("alveolar macrophages"))

ae_be <- form_list %>%
  filter(tissue %in% c("bronchial epithelium", "airway epithelium")) %>%
  View()

part1_info <- form_list %>% select(study, title, tissue, platform, mismatch_sex_lab, `additional phenoty[es`) %>%
  rename(covars=`additional phenoty[es`, num_mismatch=mismatch_sex_lab)

part1_counts <- form_list %>% select(study, samples, smokers, `non-smokers`,
                                     contains("frac"))
annot2 <- annot_part2 %>% select(tissue, study_acc, num_mismatch, covars)
part2_info <- annot2 %>% left_join(df3 %>% select(study_acc, title, gpl)) %>% rename(study=study_acc, platform=gpl) 
part2_counts <- df3 %>% 
  mutate(frac_female_s=S_female/S,
         frac_female_ns=NS_female/NS) %>%
  rename(study=study_acc, samples=total, smokers=S, `non-smokers`=NS) %>%
  select(colnames(part1_counts))

part1.2 <- part1_info %>% mutate(
  tissue=case_when(
    tissue=="bronchial epithelium" ~ "airway epithelium",
    tissue=="whole blood" ~ "blood - whole",
    tissue=="pbmc" ~ "blood - pbmcs",
    tissue =="b cell" ~ "blood - b cells",
    tissue =="leukocyte" ~ "blood - leukocytes",
    tissue=="trachea" ~ "trachea epithelium",
    tissue=="oral" ~ "oral cavity",
    tissue=="hippocampus" ~ "brain - hippocampus",
    TRUE ~ tissue
  )
)

comb_info <- part1.2 %>%
  bind_rows(part2_info %>% select(colnames(part1.2))) %>%
  mutate(tissue=ifelse(tissue=="brain", "brain - prefrontal cortex", tissue)) %>%
  select(tissue, everything()) %>% 
  arrange(tissue, study)
counts <- part1_counts %>% bind_rows(part2_counts) %>% 
  filter(study!="GSE44456",study!="GSE55962") %>%
  mutate(female_NS=round(`non-smokers`*frac_female_ns), 
         female_S=round(smokers*frac_female_s),
         male_NS=`non-smokers`-female_NS,
         male_S=smokers-female_S) %>%
  dplyr::select(-contains("frac"))
counts %>% write_csv("data/supp_tables/counts_by_study.csv")
median(counts$smokers)
summary(counts)
counts2 <- counts  %>% group_by(study) %>%
  mutate(s_ns=chisq.test(c(smokers, `non-smokers`))$p.value,
         cat=chisq.test(matrix(c(female_NS, female_S, male_NS, male_S), nrow=2, byrow=T))$p.value)
counts2 %>% dplyr::select(study, smokers, `non-smokers`, s_ns, cat) %>%
  filter(s_ns < 0.05)

counts2 %>% dplyr::select(study, female_NS, female_S, male_NS, male_S, s_ns, cat) %>%
  filter(cat < 0.05)

comb_info2 <- comb_info %>% 
  inner_join(part1_counts %>% bind_rows(part2_counts), by="study") %>%
  select(-num_mismatch, -covars, everything()) %>%
  filter(study != "GSE55962") %>%
  mutate(smokers=sprintf("%s (%s)", smokers, round(frac_female_s*smokers), 0),
         `non-smokers`=sprintf("%s (%s)", `non-smokers`, round(frac_female_ns*`non-smokers`), 0)) %>%
  select(-contains("frac")) 

comb_info2 %>% write_csv("data/table1b_updated.csv")

comb_info3 <- comb_info2 %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                               "nasal epithelium", "oral cavity",
                               "buccal mucosa", "sputum", "alveolar macrophages", 
                               "lung", 
                        "blood - b cells", "blood - pbmcs", "blood - whole",
                        "liver", "kidney", 
                        "brain - prefrontal cortex", "brain - hippocampus")))


