# code for creating the updated list for filtering

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








