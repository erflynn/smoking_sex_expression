# 01d_extract_small_ds.R
#
# Code for extracting the covariate data for smaller tissues.

library('tidyverse')
library('GEOquery')
library('GEOmetadb')

# buccal mucosa
# GSE17913
# GSE40013
# GSE8987

# nasal epithelium
# GSE8987
# alveolar macrophages
# GSE46903 --> GSE2125

# implant adherent cells - GSE42288
# urothelial - GSE21142
# kidney - GSE46699

my_studies <- c("GSE17913", "GSE40013", "GSE8987", "GSE2125",
                "GSE42288","GSE21142", "GSE46699")
my_studies_str <- paste(my_studies, collapse="','")

con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
gse_info <- dbGetQuery(con, sprintf("SELECT gse.gse, title, 
submission_date, gpl, pubmed_id, contributor, summary, overall_design 
FROM gse JOIN gse_gpl ON gse.gse=gse_gpl.gse
WHERE gse.gse IN ('%s');", my_studies_str))


gsm_info <- dbGetQuery(con,  sprintf("SELECT gsm.gsm, gse_gsm.gse, title, source_name_ch1,
description, gpl, channel_count, organism_ch1, molecule_ch1, characteristics_ch1,
supplementary_file
FROM gsm JOIN gse_gsm ON
                       gse_gsm.gsm=gsm.gsm WHERE gse_gsm.gse IN ('%s');", 
                                     my_studies_str))

dbDisconnect(con)


# all human, total RNA, 1 channel
unique(gsm_info$channel_count)
unique(gsm_info$molecule_ch1)
unique(gsm_info$organism_ch1)

gsm_info2 <- gsm_info %>% 
  select(-channel_count, -organism_ch1, -molecule_ch1,
         -supplementary_file)

# look at the most common platforms
gsm_info2 %>% group_by(gpl) %>% count() %>% arrange(desc(n))

#gpl          n
#<chr>    <int>
#  1 GPL570     254  Affymetrix Human Genome U133 Plus 2.0 Array
#2 GPL11532    84 - Affymetrix Human Gene 1.1 ST Array
#3 GPL6244     47 - Affymetrix Human Gene 1.0 ST Array
#4 GPL10274    24 - Affymetrix GeneChip Human Genome U133 Plus 2.0 Array
#5 GPL571      15  - Affymetrix Human Genome U133A 2.0 Array
#6 GPL96       10 - Affymetrix Human Genome U133A Array

download_info <- gsm_info %>% select(gsm, gse, gpl, supplementary_file)
download_info %>% write_csv("data/download_info_small_n.csv")
download_info %>% select(supplementary_file) %>% 
  separate_rows(supplementary_file, sep=";\t") %>%
  filter(str_detect(supplementary_file, "CEL")) %>%
  write_tsv("data/small_n_list_files.csv", col_names=FALSE)

no_smok <- gsm_info2 %>%
  filter(!str_detect(characteristics_ch1, "S|NS|smok|Smok"),
         !str_detect(source_name_ch1, "S|NS|smok|Smok"),
         !str_detect(title, "S|NS|smok|Smok"),
         !str_detect(description, "S|NS|smok|Smok"))  # 14
no_smok_gses <- no_smok %>% distinct(gse) %>% pull(gse) # 1
smok <- gsm_info2 %>%
  filter(str_detect(characteristics_ch1, "S|NS|smok|Smok") |
           str_detect(source_name_ch1, "S|NS|smok|Smok") |
           str_detect(title, "S|NS|smok|Smok") |
           str_detect(description, "S|NS|smok|Smok")) # 420

smok_gses <- smok %>% distinct(gse) %>% pull(gse) # 7

# TODO add a tissue label??

# split
gse_w_info <- gsm_info2 %>% group_split(gse) # --> 7 studies
phe_tabs <- lapply(gse_w_info, get_pheno_table) 
names(phe_tabs) <- unique(smok$gse)

# go through each

# ------ GSE17913 ------ #
# buccal mucosa
phe_tab1 <- phe_tabs[["GSE17913"]]
phe_tab1.2 <- phe_tab1 %>%
  select(-description) %>%
  rename(smok=`smoking status`,
         metadata_sex=gender) %>%
  # extract number from title
  mutate(id=str_extract(title, "[0-9]+")) %>%
  # extract smok-status, sex from smok, title, source_name_ch1
  mutate(title_sex=ifelse(str_detect(title, "female"), "female", "male"),
         source_sex=ifelse(str_detect(source_name_ch1, "female"), "female", "male")) %>%
  mutate(across(c(smok, title, source_name_ch1),
                ~ifelse(str_detect(., "Never"), "NS", "S"))) 

# make sure it all matches - it does
table(phe_tab1.2$metadata_sex==phe_tab1.2$title_sex)
table(phe_tab1.2$metadata_sex==phe_tab1.2$source_sex)
table(phe_tab1.2$smok==phe_tab1.2$title)
table(phe_tab1.2$smok==phe_tab1.2$source_name_ch1)
phe_tab1.3 <- phe_tab1.2 %>% select(-title, -source_name_ch1, -title_sex, -source_sex)

phe_tab1.3 %>% group_by(smok, metadata_sex) %>% count()
phe_tab1.3 %>% write_csv("data/pdata/gse17913.csv")


# ----- "GSE40013" ----- #
# saliva
phe_tab2 <- phe_tabs[["GSE40013"]]
table(phe_tab2$source_name_ch1)
phe_tab2.2 <- phe_tab2 %>% 
  rename(metadata_sex=gender,
         smok=`smoking status`,
         race_ethnicity=race,
         stress=`chronic stress level`) %>%
  mutate(smok=case_when(
    smok=="former" ~ "FS",
    smok=="current" ~ "S",
    smok=="never" ~ "NS"
  )) %>%
  mutate(age=as.numeric(age)) %>%
  # unclear if this is an ID or what?
  mutate(id=str_extract(title, "[0-9]+")) %>% # I think for pairing the data?
  # extract stress level from title, etc
  mutate(across(c(title, source_name_ch1), 
                ~ifelse(str_detect(.,"high"), "high", "low")))

# compare - does it match?
table(phe_tab2.2$title==phe_tab2.2$stress)
table(phe_tab2.2$source_name_ch1==phe_tab2.2$stress)
phe_tab2.3 <- phe_tab2.2 %>% 
  select(-title, -source_name_ch1) %>%
  rename(id1=description,
         id2=id)

phe_tab2.3 %>% 
  select( age, metadata_sex,race_ethnicity, smok, stress) %>%
  mutate(across(c(metadata_sex,race_ethnicity, smok, stress), as.factor)) %>%
  summary()

# small counts - this is complicated by stress level... and mostly non-smokers 
phe_tab2.3 %>% filter(smok %in% c("NS", "S")) %>%
  group_by(metadata_sex, smok, stress) %>%
  count()

phe_tab2.3 %>% filter(smok %in% c("NS", "S")) %>%
  group_by(metadata_sex, smok) %>%
  count()
# TODO: figure out more abt ID, experiment setup
phe_tab2.3 %>% write_csv("data/pdata/gse40013.csv")


# --------- GSE46903 ----- #
# this is peripheral BLOOD
# GSE13896 (COPD patients, smokers and non-smokers)  - already downloaded
# GSE2125 (asthmatic patients, smokers and non-smokers) 
#phe_tab3 <- phe_tabs[["GSE46903"]]
# .... we aren't interested in this, but try GSE2125

# ------- "GSE2125" ------- #
phe_tab3 <- phe_tabs[["GSE2125"]]
# alveolar macrophages
phe_tab3.1 <- phe_tab3 %>%
  select(-source_name_ch1, -description) %>%
  mutate(tissue="alveolar macrophages") %>%
  rename(smok=status) %>%
  mutate(smok=case_when(
    smok=="Smoker" ~ "S",
    smok=="Nonsmoker" ~ "NS",
    smok=="Asthmatic" ~ "asthma")) %>%
  mutate(has_asthma=ifelse(smok=="asthma", "y", "n")) %>%
  mutate(batch=ifelse(str_detect(title, "Koth"), "Koth", "Woodruff"))
# //TODO: does this overlap with the other Woodruff study??
phe_tab3.1 %>% group_by(smok) %>% count()
phe_tab3.1 %>% write_csv("data/pdata/gse2125.csv")

# ------- "GSE8987" ----- #
phe_tab4 <- phe_tabs[["GSE8987"]]
phe_tab4.2 <- phe_tab4 %>%
  # remove parsed badly
  select(-`Current Smoker`, -`Never Smoker`) %>%
  mutate(id=str_extract(title, "[0-9]+")) %>%
  mutate(across(c(source_name_ch1, description, title), tolower)) %>%
  mutate(source_tissue=source_name_ch1,
         descript_tissue=description,
         title_tissue=title) %>%
  rename(source_smok=source_name_ch1,
         descript_smok=description,
         title_smok=title) %>%
  mutate(across(contains("tissue"), ~ifelse(
    str_detect(., "buccal mucosa") | str_detect(., "mouth"), "buccal mucosa", 
    "nasal epithelium"))) %>%
  mutate(across(contains("smok"), ~ifelse
    (str_detect(., "never"), "NS", "S")))

# check matches  -- all matches
table(phe_tab4.2$title_smok==phe_tab4.2$source_smok)
table(phe_tab4.2$title_smok==phe_tab4.2$descript_smok)
table(phe_tab4.2$title_tissue==phe_tab4.2$source_tissue)
table(phe_tab4.2$title_tissue==phe_tab4.2$descript_tissue)

phe_tab4.3 <- phe_tab4.2 %>% 
  select(-contains("source"), -contains("descript")) %>%
  rename(smok=title_smok,
         tissue=title_tissue)

# small study...
phe_tab4.3 %>% group_by(smok, tissue) %>% count()

phe_tab4.3 %>% write_csv("data/pdata/gse8987.csv")



# ---- GSE42288 ---- #
# implant adherent cells --> discard this?
phe_tab5 <- phe_tabs[["GSE42288"]]
head(unique(phe_tab5$title))
head(unique(phe_tab5$source_name_ch1))
phe_tab5.1 <- phe_tab5 %>% 
  select(-description) %>%
  unite(tissue, c(tissue, `cell type`), sep=": ") %>%
  rename(smok=subject) %>%
  mutate(title_day=tolower(str_extract(title, "DAY[0-9]+"))) %>%
  mutate(source_day=str_extract(source_name_ch1, "day[0-9]+")) %>%
  mutate(surface=case_when(
    str_detect(source_name_ch1, "micro") ~ "microroughsurface",
    str_detect(source_name_ch1, "nano") ~ "nanosurface")) %>%
  mutate(type=case_when(
    str_detect(title, "Osseospeed") ~ "osseospeed",
    str_detect(title, "TiOblast") ~ "tioblast"
  )) %>%
  mutate(id=str_trim(str_extract(title, "[0-9]+ "))) %>%
  mutate(title_smok=ifelse(str_detect(title, "Smoking"), "S", "NS")) %>%
  rename(source_smok=source_name_ch1) %>%
  mutate(across(c(smok, source_smok), ~ifelse(str_detect(., "Non-"), "NS", "S")))

table(phe_tab5.1$title_day==phe_tab5.1$source_day)
table(phe_tab5.1$smok==phe_tab5.1$source_smok) 
table(phe_tab5.1$smok==phe_tab5.1$title_smok) 

phe_tab5.2 <- phe_tab5.1 %>% 
  rename(day=title_day) %>% 
  select(-contains("source"), -contains("title"))

phe_tab5.2 %>% write_csv("data/pdata/gse42288.csv")

# ---- GSE21142 ---- #
# urothelial 
smok_status <- function(x) {
  case_when(
  str_detect(x, "former smoker") ~ "FS",
  str_detect(x, "non-smoker") ~ "NS",
  str_detect(x, "smoker") ~ "S"
)}
phe_tab6 <- phe_tabs[["GSE21142"]]
table(phe_tab6$title==phe_tab6$source_name_ch1)
phe_tab6.1 <- phe_tab6 %>% 
  select(-source_name_ch1) %>%
  rename(tissue_type=type,
         smok=`smoking habit`, 
         carcinoma_type=disease) %>%
  mutate(id=str_extract(title, "Nr[0-9]+")) %>%
  mutate(title_smok=title,
         descript_smok=description) %>%
  mutate(across(c(smok, title_smok, descript_smok), ~smok_status(.))) %>%
  mutate(across(c(title, description), ~ifelse(str_detect(., "superficial"), 
                                               "superficial urothelial carcinoma",
                                               "invasive urothelial carcinoma")))

table(phe_tab6.1$title==phe_tab6.1$carcinoma_type)
table(phe_tab6.1$description==phe_tab6.1$carcinoma_type)
table(phe_tab6.1$smok==phe_tab6.1$title_smok)
table(phe_tab6.1$smok==phe_tab6.1$descript_smok)

phe_tab6.2 <- phe_tab6.1 %>% select(-contains("title"), -contains("descript"))
phe_tab6.2 %>%
  write_csv("data/pdata/gse21142.csv")

# ---- GSE46699 ---- #
# kidney
phe_tab7 <- phe_tabs[["GSE46699"]]
table(phe_tab7$description)


phe_tab7.1 <- phe_tab7 %>%
  select(-description) %>%
  rename(tissue_type=tissue,
         smok=smoking,
         id=patient,
         tissue=source_name_ch1) %>%
  mutate(smok=ifelse(smok=="no", "NS", "S")) %>%
  mutate(title_id=str_trim(str_extract(title, " [0-9]+ "))) %>%
  mutate(title_tissue_type=ifelse(str_detect(title, "normal"), "normal", "tumor")) %>%
  mutate(title_smok=ifelse(str_detect(title, "yes smoking"), "S", "NS"))

table(phe_tab7.1$tissue_type==phe_tab7.1$title_tissue_type)
table(phe_tab7.1$id==phe_tab7.1$title_id)
table(phe_tab7.1$smok==phe_tab7.1$title_smok)
phe_tab7.2 <- phe_tab7.1 %>%
  select(-contains("title")) 
phe_tab7.2 %>% write_csv("data/pdata/gse46699.csv")
