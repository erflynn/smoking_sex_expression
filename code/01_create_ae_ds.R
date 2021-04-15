
library(tidyverse)
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
# save(data1, file="../data_ae.RData")
# 
# eset1 <- affy::rma(data1)
# save(eset1, file="../eset_ae.RData")

# -- DEDUPLICATE
load("data/eset_ae.RData") # --> eset1
expDat <- exprs(eset1) # 736

eset1_mat <- expDat[,!duplicated(t(expDat))] # 595

# figure out which are duplicated...
eset1_gene2 <- data.frame(t(expDat[1:2,order(expDat[1,])]))
eset1_gene2$sample <- rownames(eset1_gene2)
dup_split <- eset1_gene2 %>% 
  as_tibble() %>% 
  group_by(`X1007_s_at`) %>% 
  summarize(n=n(), sample=paste(str_extract(sample, "^[0-9A-Za-z]+"), 
                                collapse=";")) %>%
  filter(n>=2) %>%
  arrange(desc(n), sample) %>%
  dplyr::select(-`X1007_s_at`)



condense_cols <- function(x){
  paste(sort(unique(x[!is.na(x)])), collapse=";")
}
ae_only_meta <- read_csv("data/ae_only_meta.csv")

comb_metadata <- read_csv("data/ae_full_gsm.csv") %>% 
  rename(geo_accession=gsm) %>%
  mutate(sample_acc=geo_accession) %>%
  add_sl() %>%
  select(-sample_acc) %>%
  left_join(
    ae_only_meta %>%
      mutate(dgm_id=tolower(str_extract(id, "DGM-[0-9]+")) ) %>%
      select(geo_accession, dgm_id)
  ) %>% 
  unite(dgm_id,c(dgm_id, `department of genetic medicine id`), na.rm=T ) %>%
  group_by(geo_accession) %>%
  mutate(dgm_id=unique(str_split(dgm_id, "_")[[1]])) %>%
  ungroup() %>%
  mutate(across(everything(), ~ifelse(.=="", NA, .)))


dup_dat <- dup_split %>% 
  mutate(grp=1:n()) %>%
  select(-n) %>%
  rename(geo_accession=sample) %>%
  separate_rows(geo_accession) %>% 
  full_join(comb_metadata, by=c("geo_accession")) %>%
  mutate(grp=ifelse(is.na(grp), geo_accession, grp)) %>%
  group_by(grp) %>%
  summarize(across(everything(), ~condense_cols(.))) 

# for duplicated columns, put together the phenotype data
dup_dat %>% filter(str_detect(sex,";") |
                     str_detect(race_ethnicity,";") |
                     str_detect(smoking,";") | # one
                     str_detect(age,";")  | #one
                     str_detect(sex_lab, ";")) 
# age we keep
# do not keep the S/NS

dup_dat2 <- dup_dat %>%
  mutate(
    across(c(age, pack_years), ~ case_when(
      ! str_detect(., ";") | is.na(.) ~ .,
      TRUE ~ as.character(mean(as.numeric(str_split(., ";")[[1]])))
    )) # take the mean of the ages
  ) %>%
  filter(!str_detect(smoking, ";"))  # remove conflicting
# STOPPED

dgm_id_to_grp <- dup_dat2 %>% 
  group_by(dgm_id) %>%
  count() %>% filter(dgm_id!="" & n>1)

condensed_id <- dup_dat2 %>% semi_join(dgm_id_to_grp) %>%
  group_by(dgm_id) %>%
  mutate(across(everything(), ~ifelse(.=="", NA, .))) %>%
  summarize(across(everything(), ~condense_cols(.)))

dup_dat3 <- dup_dat2 %>% 
  anti_join(dgm_id_to_grp) %>%
  bind_rows(condensed_id) # 464

dup_dat3$keep_acc <- sapply(dup_dat3$geo_accession, function(x) str_split(x, ";")[[1]][[1]])
length(unique(dup_dat3$keep_acc )) # 466
nrow(dup_dat3) # 466
colnames(expDat) <- str_extract_all(colnames(expDat), "^[0-9A-Za-z]+")
dup_dat4 <- dup_dat3 %>% filter(keep_acc %in% colnames(expDat)) # 457


expDat2 <- as.matrix(expDat[,dup_dat4$keep_acc]) # 457
pDat2 <- dup_dat4 %>% 
  mutate(geo_accession=keep_acc) %>% 
  select(-grp, -keep_acc) %>%
  mutate(sex=case_when(
    sex=="f" ~ "female",
    sex=="m" ~ "male")) %>%
  mutate(across(everything(), ~ifelse(.=="", NA, .)))

fct_summ(pDat2)
pDat2 %>% filter(sex!=sex_lab & sex_lab!="unlabeled") # 11 are mismatched

# 3. redo sex lab to confirm

# download GPL570
gse4 <- getGEO("GSE2125")
f_df <- fData(gse4$GSE2125_series_matrix.txt.gz)

# toker
xist_gpl570 <- f_df %>% filter(str_detect(`Gene Symbol`, "XIST")) %>% pull(ID)
rps_gpl570 <- f_df %>% filter(str_detect(`Gene Symbol`, "RPS4Y1")) %>% pull(ID)
kdm_gpl570 <- f_df %>% filter(str_detect(`Gene Symbol`, "KDM5D")) %>% pull(ID)
ts <- tokerSexLab(expDat2, f.genes =xist_gpl570, m.genes=c(rps_gpl570, kdm_gpl570))
xmeans <- colMeans(expDat2[xist_gpl570,])
ymeans <- colMeans(expDat2[c(rps_gpl570, kdm_gpl570),])
dat_info <- tibble(x=xmeans, y=ymeans, sample_acc=colnames(expDat2))
dat_info2 <- dat_info %>% 
  left_join(pDat2 %>% select(geo_accession, sex, sex_lab), by=c("sample_acc"="geo_accession")) %>%
  rename(metadata_sex=sex, expr_sex=sex_lab)%>%
  mutate(across(contains("sex"), ~ifelse(is.na(.), "unlabeled", .))) %>%
  left_join(
    tibble(
      sample_acc=names(ts),
      toker_sex=unlist(ts)
    )
  )
ggplot(dat_info2, aes(x=x, y=y, col=metadata_sex, shape=expr_sex))+
  geom_point(alpha=0.7)+
  theme_bw()+
  xlab("xist probes")+
  ylab("kdm5d + rps4y1 probes")
ggsave("figures/paper_figs/ae_sex_lab_viz.png")
dat_info2 %>% filter(expr_sex!=toker_sex & expr_sex!="unlabeled")

# 4. remove conflicting
# remove 11 samples that have mismatched labels + 3 that are outside clusters
# 1 mismatched
to_remove_sl <- dat_info2 %>%
  filter((toker_sex != metadata_sex &
           metadata_sex != "unlabeled") |
           (5 < x & x < 7)) %>%
  pull(sample_acc)
dat_info3 <- dat_info2 %>%
  filter(toker_sex==metadata_sex | metadata_sex == "unlabeled") %>%
  filter(x < 5 | x > 7)

pDat3 <- pDat2 %>% inner_join(dat_info3 %>% 
                      select(sample_acc, toker_sex),
                     by=c("geo_accession"="sample_acc"))
# 444

# clean up, add submission date, save
pDat4 <- pDat3 %>% 
  select(-sex, -sex_lab) %>% 
  rename(sex=toker_sex) %>%
  left_join(download_info %>% select(gsm, submission_date), by=c("geo_accession"="gsm"))
expDat4 <-expDat2[,pDat4$geo_accession]
fct_summ(pDat4) # 105 missing RE, 109 missing age
save(expDat4, pDat4, file="data/ae_full_exp.RData")
pDat4 %>% write_csv("data/ae_pdat_full.csv")

# 4. PC plot to demo effects

# 5. table 1


# 6. infer missing labels?
# decide which variable selection method to use?
# - stepwise or glmnet, correlation, varImp?
# - 
# - how to iteratively chain