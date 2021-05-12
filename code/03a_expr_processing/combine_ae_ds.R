# Code for putting the AE ds together into one ds
#

# others: GSE43949. GSE7895, GSE37147

library('tidyverse')
library('readxl')


# yang et al version
supp_files1 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=1, skip=2, col_names=TRUE)
incl_files <- supp_files1 %>%
  filter(`Microarray Platform`=="U133 Plus 2.0 Array")
#supp_files1 %>%
#  filter(`Microarray Platform`!="U133 Plus 2.0 Array")
#  GSE7895 - GPL96
#  GSE37147 - gene chip
# TODO: do I add these somewhere?

list_geo_studies <- incl_files %>% pull(`GEO Accession`)
ae_phe <- read_csv("data/ae_phe_data.csv") # 365
fct_summ(ae_phe) # only complete
ae_phe2 <- read_csv("data/rep_full_gsm.csv") # with incomplete!

# rb version
load("data/ae_eset_pheno.RData") # --> eset_comb3, comb_metadata3 (after dedup)
fct_summ(comb_metadata3)
new_gsm <- comb_metadata3 %>% 
  filter(str_detect(tissue, "airway")) %>% # 484
  anti_join(ae_phe2, by=c("geo_accession"="gsm"))

# note -- missing expr_sex?

# ae GSEs
ae_gses <- comb_metadata3 %>% 
  filter(str_detect(tissue, "airway")) %>% 
  distinct(study) %>%
  separate_rows(study, sep=";") %>%
  distinct(study) %>%
  pull(study)
setdiff(list_geo_studies, ae_gses) # "GSE5058"  "GSE8545"  "GSE13933" "GSE4498"  "GSE7832" 
setdiff(ae_gses, list_geo_studies) # "GSE64614"  "GSE108134" "GSE43949"  "GSE19722" 



meta <- read_csv("data/comb_sample_epithelium_metadata.csv")
ae_gses_before_dedup <- meta %>%
  filter(str_detect(tissue, "airway")) %>% 
  distinct(study) %>%
  separate_rows(study, sep=";") %>%
  distinct(study) %>%
  pull(study)
setdiff(ae_gses_before_dedup, ae_gses) # GSE4302
setdiff( ae_gses, ae_gses_before_dedup) # none


# - which of these are GPL96?
new_gses <- setdiff(ae_gses_before_dedup, list_geo_studies)
all_gses <- union(ae_gses_before_dedup, list_geo_studies)
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
list_studies <- paste(all_gses, collapse="','")
plat <- dbGetQuery(con, sprintf("SELECT gse,gpl FROM gse_gpl WHERE gse IN ('%s');",
                                list_studies))
dbDisconnect(con)
plat # GSE43949 - not GPL570 TODO include
gpl570_gses <- plat %>% filter(gpl=="GPL570") %>% pull(gse)
setdiff(gpl570_gses, list_geo_studies) # 4 studies
# "GSE108134" "GSE19722"  "GSE4302"   "GSE64614" 

# are any of them entirely encompassed by the others? NOPE
meta %>% distinct(study) %>% filter(str_detect(study, "GSE108134"))
meta %>% distinct(study) %>% filter(str_detect(study, "GSE19722"))
meta %>% distinct(study) %>% filter(str_detect(study, "GSE4302"))
meta %>% distinct(study) %>% filter(str_detect(study, "GSE64614"))


# next step:
# 0. make a list of new ones to download
# get the lists of files associated with these data

library(GEOmetadb)
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
gse_gsm <- dbGetQuery(con, sprintf("SELECT gse, gsm FROM gse_gsm WHERE gse IN ('%s');", 
                                   paste(list_studies, collapse="','")))

gsm2 <- dbGetQuery(con, sprintf("SELECT gsm, title, source_name_ch1, description, 
characteristics_ch1, gpl, submission_date, supplementary_file FROM gsm 
                        WHERE gsm IN ('%s');", paste(unique(gse_gsm$gsm), collapse="','")))
dbDisconnect(con)

# additional info
download_info <- gsm2 %>% select(gsm, gpl, 
                                 submission_date, 
                                 supplementary_file)

# clean up AE phenotype data
gsm2.1 <- gsm2 %>%
  select(-gpl, -submission_date, -supplementary_file) %>%
  separate_rows(characteristics_ch1, sep=";\t") %>%
  mutate(characteristics_ch1=tolower(characteristics_ch1)) %>%
  separate(characteristics_ch1, into=c("key", "value"), sep=": ") %>%
  dplyr::select(-title, -source_name_ch1, -description) %>%
  pivot_wider(names_from=key, values_from=value)

gsm2_kv <- gsm2 %>%
  select(-gpl, -submission_date, -supplementary_file) %>%
  separate_rows(characteristics_ch1, sep=";\t") %>%
  mutate(characteristics_ch1=tolower(characteristics_ch1)) %>%
  filter(str_detect(characteristics_ch1, ": ")) %>%
  separate(characteristics_ch1, into=c("key", "value"), sep=": ") %>%
  dplyr::select(-title, -source_name_ch1, -description) 

gsm2_v <- gsm2 %>%
  select(-gpl, -submission_date, -supplementary_file) %>%
  separate_rows(characteristics_ch1, sep=";\t") %>%
  mutate(characteristics_ch1=tolower(characteristics_ch1)) %>%
  filter(!str_detect(characteristics_ch1, ": ")) %>%
  pivot_longer(-gsm, names_to="key", values_to="value")

gsm2.1 <- gsm2_kv %>% 
  bind_rows(gsm2_v)  %>%
  pivot_wider(names_from=key, values_from=value)
  

# TODO: check for pheno duplicates with other vals!
gsm2.2 <- gsm2.1 %>%
  mutate(race_ethnicity=case_when(
    `ethnic group` == "hispnaic" ~ "hispanic",
    `ethnic group`=="afr" ~ "black",
    `ethnic group`=="eur" ~ "white",
    `ancestry`=="african" ~ "black",
    `ancestry`=="european" ~ "white",
    `ancestry`=="hispanic" ~ "hispanic",
    `ethnicity`=="afr" ~ "black",
    `ethnicity`=="eur" ~ "white",
    TRUE ~ `ethnic group`
  )) %>%
  dplyr::select(-ethnicity, -ancestry, -`ethnic group`) %>%
  separate(`smoking status`, into=c("smoking", "pack_years"), 
           sep=", ", extra="merge") %>% 
  mutate(copd=case_when(
    smoking == "copd smoker" ~ "yes",
    `copd status`=="yes" ~ "yes",
    smoking == "copd" ~ "yes",
    smoking == "early-copd" ~ "early",
    TRUE  ~ "no"
  )) %>%
  mutate(smoking=case_when(
    smoking %in% c("non-smoker", "nonsmoker", "ns") ~ "NS",
    smoking %in% c("smoker", "s", "copd smoker") ~ "S"
  )) %>%
  dplyr::select(-`copd status`) %>%
  mutate(pack_years=as.numeric(str_replace_all(pack_years, " pack-years", ""))) %>%
  dplyr::select(gsm, age, sex, smoking, race_ethnicity, 
                copd, pack_years, everything())

gsm2.2 %>% select(description, title, source_name_ch1) %>% fct_summ() 
# nothing in description or source_name_ch1

#gsm2.2$title[!is.na(gsm2.2$title)]
#gsm2.2 %>% filter(str_detect(title, "non-smoker")) %>% pull(smoking)

# clean up based on title, source_name_ch1, description, characteristics_ch1
gsm2.3 <- gsm2.2 %>%
  mutate(smoking=ifelse(str_detect(characteristics_ch1,"smoker") & 
                          !is.na(characteristics_ch1), "S", 
                                   smoking)) %>%
  mutate(smoking=ifelse(str_detect(title, "non-smoker") & !is.na(characteristics_ch1), 
                        "NS", smoking)) %>%
  mutate(smoking=ifelse(characteristics_ch1=="healthy control" & 
                          !is.na(characteristics_ch1), "NS", smoking)) %>%
  mutate(pack_years=ifelse(str_detect(characteristics_ch1, "1 pack-years") & !is.na(characteristics_ch1),
                           "1", pack_years)) %>%
  mutate(asthma_treatment=ifelse(str_detect(characteristics_ch1, "asthma"), characteristics_ch1, NA))   %>%
  select(-description, -source_name_ch1, -characteristics_ch1, -title)


# remove copd, trachea, and repeated measures
gsm2.4 <- gsm2.3 %>%
  filter(is.na(`time of serial bronchoscopy (month)`) | 
           `time of serial bronchoscopy (month)`=="m0") %>%
  filter(is.na(strain), is.na(asthma_treatment)) %>%
  select(-strain, -`cell type`, -`chip antibody`, -`chip antibody manufacturer`,
         -`time of serial bronchoscopy (month)`, -`asthma_treatment`) %>%
  filter(copd=="no") %>%
  select(-copd) %>%
  mutate(age=ifelse(age=="unknown", NA, age)) %>%
  filter(!is.na(smoking)) %>%
  filter(is.na(tissue) | !str_detect(tissue, "trachea")) %>%
  select(-tissue)

trachea_gsm <- gsm2.3 %>% filter(str_detect(tissue, "trachea")) %>% pull(gsm)

# remove GSE4302 b/c other lab
gse4302_gsm <- meta %>% filter(study=="GSE4302")  %>% pull(geo_accession)
gsm2.5 <- gsm2.4 %>% filter(!gsm %in% gse4302_gsm)


gsm2.5 %>% write_csv("data/ae_full_gsm.csv") # 743

# remove previously downloaded
list_downloaded <- read_csv("data/list_ae_to_download_replication.csv")
rem_download <- download_info %>% 
  filter(gpl=="GPL570") %>%
  anti_join(list_downloaded, by="gsm")
rem_download %>% write_csv("data/list_ae_rem_download.csv")

# 1. download new ones
# ---> download server side
# while IFS=',' read -r gsm gpl submission_date fpath
# do
# wget "$fpath";
# done < "../list_ae_rem_download.csv"
# <--- read in

# 2. put together new ones
# copy over previously downloaded
#  cp cel_files_rep/* cel_files_rem_ae/

# filter files
# keep_f <- read_csv("data/ae_full_gsm.csv")
# setwd("data/cel_files_rem_ae/")
# my_f <- list.files()
# gsm_name <- lapply(my_f, function(x) str_split(str_replace_all(x, ".CEL.gz", ""), "_")[[1]][[1]])
# names(my_f) <- gsm_name
# to_remove <- setdiff(gsm_name, keep_f$gsm) 
# missing <- setdiff(keep_f$gsm, gsm_name) # 7
# 
# lapply(my_f[unlist(to_remove)], function(fname){
#   system(sprintf("rm %s", fname))
# })


# ---> load server side
# library(affy)
# 
# setwd("data/cel_files_rem_ae/")
# data1 <- ReadAffy(filenames=list.celfiles()) 
# save(data1, file="../data_rep.RData")
# 
# eset1 <- affy::rma(data1)
# save(eset1, file="../eset_rep.RData")


# 3. deduplicate


# ---- other GSEs ---- #
# pull table 1 data for these: only GSE4302 + GSE994
load("data/study_sample_info.RData") # --> all_stats
incl_studies <- read_csv("data/incl_studies.csv")
my_studies <- all_stats[incl_studies$study_acc]

# AE other platforms: GSE43949 - actually MOUSE

# AE other lab: GSE4302
gse4302 <- my_studies[["GSE4302"]]$df %>%
  filter(characteristics_ch1 %in% c("smoker", "healthy control")) %>%
  mutate(smok=ifelse(characteristics_ch1=="smoker", "S", "NS")) %>%
  select(study_acc, sample_acc, smok) %>%
  mutate(tissue="airway epithelium") %>%
  add_sl()
gse4302 %>% write_csv("data/pdata_filt/gse4302.csv")

#fct_summ(my_studies[["GSE7895"]]$df) # actually BE - already present
# fct_summ(my_studies[["GSE37147"]]$df)  # also BE - all current or former


# TODO - other GSEs, which are already included in small ds?
other_gses <-  comb_metadata3 %>% 
  filter(!str_detect(tissue, "airway")) %>% 
  group_by( tissue, study, smok) %>%
  count() %>%
  pivot_wider(names_from="smok", values_from="n", values_fill=0)
other_gses

# GSE19027 bronchial epithelium - already present

# GSE994   bronchial epithelium - removed overlap and added
fct_summ(all_stats[["GSE994"]]$df) # overlap with
intersect(all_stats[["GSE19027"]]$df$sample_acc, all_stats[["GSE994"]]$df$sample_acc) # only 10
gse994 <- all_stats[["GSE994"]]$df %>% 
  filter(!sample_acc %in% all_stats[["GSE19027"]]$df$sample_acc) %>%
  select(-description) %>%
  mutate(smok=case_when(
    str_detect(title, "current") ~ "S",
    str_detect(title, "never") ~ "NS",
    str_detect(title, "former") ~ "FS"
  )) %>%
  select(study_acc, sample_acc, smok) %>%
  add_sl()  %>% 
  mutate(tissue="trachea")

gse994 %>% group_by(sex_lab, smok) %>% count()
gse994 %>% write_csv("data/pdata_filt/gse994.csv")

# GSE11906 trachea epithelium  
#fct_summ(my_studies[["GSE11906"]]$df ) # mostly AE, only 8 trachea

# GSE64614 trachea epithelium  - all nonsmokers

# gse <- getGEO("GSE64614")
# pDat <- pData(gse$GSE64614_series_matrix.txt.gz)
# pDat2 <- pDat %>% mutate(source_name_ch1=tolower(source_name_ch1)) %>% filter(str_detect(source_name_ch1, "trachea"))
# pDat2 %>% select(geo_accession, contains("characteristics_ch1")) %>%
#   pivot_longer(-geo_accession) %>%
#   select(-name) %>%
#   filter(value!="") %>%
#   separate(value, into=c("key", "value"), sep=": ") %>%
#   pivot_wider(names_from="key", values_from="value") %>%
#   fct_summ()

# GSE13896 alveolar macrophages - already present
# GSE13931 alveolar macrophages = AE, AM already present in GSE13896

# GSE8823  alveolar macrophages - excluded b.c overlap


