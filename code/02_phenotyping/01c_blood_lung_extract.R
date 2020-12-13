# 01c_blood_lung_extract.R
#
# Updated phenotyping code that works by extracting information directly
# from GEOmetadb instead of downloading all the datasets individually.
#
# Note: we are also including the "maybe" and "smoking history"
# datasets here.
#
# This only has one example study completed - lots to do!


library('tidyverse')
library('GEOquery')
library('GEOmetadb')

# clean up the phenotype data
get_pheno_table <- function(df){
  df %>% 
    separate_rows(characteristics_ch1, sep=";\t") %>%
    separate(characteristics_ch1, into=c("attr", "value"), sep=":") %>%
    mutate(value=str_squish(value)) %>%
    pivot_wider(names_from=attr, values_from=value)
}


lb_studies <- read_tsv("data/list_lung_blood.tsv", col_names=F)
lb_studies_str <- paste(lb_studies$X1, collapse="','")

# connect and grab all the GSE info
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
gse_info <- dbGetQuery(con, sprintf("SELECT gse.gse, title, 
submission_date, gpl, pubmed_id, contributor, summary, overall_design 
FROM gse JOIN gse_gpl ON gse.gse=gse_gpl.gse
WHERE gse.gse IN ('%s');", lb_studies_str))


# grab all the GSM info
# the characteristics_ch1 column is what has the covariate info
gsm_info <- dbGetQuery(con,  sprintf("SELECT gsm.gsm, gse_gsm.gse, title, source_name_ch1,
description, gpl, channel_count, organism_ch1, molecule_ch1, characteristics_ch1,
supplementary_file
FROM gsm JOIN gse_gsm ON
                       gse_gsm.gsm=gsm.gsm WHERE gse_gsm.gse IN ('%s');", 
                       lb_studies_str))

dbDisconnect(con)

# all human, total RNA, 1 channel
gsm_info %>% 
    select(gsm, gse, gpl, supplementary_file) %>% 
    write_csv("data/download_info_lung_blood.csv")
gsm_info2 <- gsm_info %>% 
  select(-channel_count, -organism_ch1, -molecule_ch1,
         -supplementary_file)

# look at the most common platforms
gsm_info2 %>% group_by(gpl) %>% count() %>% arrange(desc(n))
# # A tibble: 11 x 2
# # Groups:   gpl [11]
# gpl          n
# <chr>    <int>
#   1 GPL10904  1082 - Illumina HumanHT v4 expression beadchip
# 2 GPL570     440
# 3 GPL10399   253  - Illumina HumanRef-8 v3.0
# 4 GPL13667   159 - [HG-U219] Affymetrix Human Genome U219 Array
# 5 GPL6102     97 - Illumina human-6 v2.0 expression beadchip
# 6 GPL16384    36
# 7 GPL571      36
# 8 GPL6480     36
# 9 GPL6244     24
# 10 GPL96       23
# 11 GPL16686    10



# --- look at whether there is smoking information available --- #
no_smok <- gsm_info2 %>%
  filter(!str_detect(characteristics_ch1, "S|NS|smok|Smok"),
          !str_detect(source_name_ch1, "S|NS|smok|Smok"),
         !str_detect(title, "S|NS|smok|Smok"),
         !str_detect(description, "S|NS|smok|Smok")) # 84

no_smok_gses <- no_smok %>% distinct(gse) %>% pull(gse) # 5

smok <- gsm_info2 %>%
  filter(str_detect(characteristics_ch1, "S|NS|smok|Smok") |
         str_detect(source_name_ch1, "S|NS|smok|Smok") |
         str_detect(title, "S|NS|smok|Smok") |
         str_detect(description, "S|NS|smok|Smok")) # 2040

smok_gses <- smok %>% distinct(gse) %>% pull(gse)
intersect(no_smok_gses, smok_gses) # --> rescue these two

no_smok_keep <- no_smok %>% filter(gse %in% smok_gses )
smok2 <- smok %>% bind_rows(no_smok_keep)

# label tissue
smok3 <- smok2 %>% mutate(
  tissue2=case_when(
    str_detect(title,"Lung|lung") ~ "lung",
    str_detect(description,"Lung|lung") ~ "lung",
    str_detect(characteristics_ch1,"Lung|lung") ~ "lung",
    str_detect(source_name_ch1,"Lung|lung") ~ "lung",
    str_detect(title,"blood|Blood|PBMC|leukocytes") ~ "blood",
    str_detect(description,"blood|Blood|PBMC|leukocytes") ~ "blood",
    str_detect(characteristics_ch1,"blood|Blood|PBMC|leukocytes") ~ "blood",
    str_detect(source_name_ch1,"blood|Blood|PBMC|leukocytes") ~ "blood",
    TRUE ~ ""))
# > smok3 %>% group_by(tissue2) %>% count()
# blood: 1863 samples, 17 studies
# lung: 819 samples, 9 studies

lung_dat <- smok3 %>% filter(tissue2=="lung") %>% group_split(gse) 
blood_dat <- smok3 %>% filter(tissue2=="blood") %>% group_split(gse) 

gse_w_info <- smok2 %>% group_split(gse) # --> 26 studies

# count studies / samples
smok2 %>% group_by(gpl) %>% mutate(n_samples=n()) %>%
  group_by(gpl, n_samples) %>%
  distinct(gse) %>% summarize(n_studies=n()) %>% arrange(desc(n_studies)) 
# gpl      n_samples n_studies
# <chr>        <int>     <int>
# 1 GPL10904      1082         7 - Illumina HumanHT-12 V4.0 expression beadchip 
# 2 GPL570         529         6
# 3 GPL10399       253         3 - Illumina HumanRef-8 v3.0 
# 4 GPL96          160         3
# 5 GPL13667       159         2 - Affymetrix Human Genome U219 Array
# 6 GPL6244        327         2 - Affymetrix Human Gene 1.0 ST Array
# 7 GPL6102         97         1 - Illumina human-6 v2.0 expression beadchip	
# 8 GPL6480         36         1 - Agilent-014850 Whole Human Genome Microarray 4x44K G4112F 
# 9 GPL6884         39         1 - Illumina HumanWG-6 v3.0 expression beadchip


phe_tabs <- lapply(gse_w_info, get_pheno_table) 
names(phe_tabs) <- unique(smok2$gse)
lung_tabs <-lapply(lung_dat, get_pheno_table)
names(lung_tabs) <- unique(smok3 %>% filter(tissue2=="lung") %>% pull(gse))
save(lung_tabs, file="data/lung_pheno_data.RData")

blood_tabs <-lapply(blood_dat, get_pheno_table)
names(blood_tabs) <- unique(smok3 %>% filter(tissue2=="blood") %>% pull(gse))
save(blood_tabs, file="data/blood_pheno_data.RData")

# look at one particular study
gse10072 <- lung_tabs["GSE10072"][[1]] %>%
  mutate(tissue="lung") %>%
  dplyr::rename('smok'=`Cigarette Smoking Status`,
                'age'=`Age at Diagnosis`,
                'metadata_sex'=Gender,
                "cancer_stage"=Stage) %>%
  mutate(metadata_sex=tolower(metadata_sex),
         cancer_type="adenocarcinoma",
         cancer=ifelse(source_name_ch1=="Adenocarcinoma of the Lung", "y", "n"),
         smok=case_when(
           smok=="Never Smoked" ~ "NS",
           smok=="Current Smoker" ~ "S",
           smok=="Former Smoker" ~ "FS"
         )) %>%
  group_by(gsm) %>%
  mutate(id=str_split(title, "_")[[1]][[2]]) %>%
  select(-source_name_ch1, -title, -description)
gse10072 %>% mutate(across(c(smok, metadata_sex, cancer), as.factor),
                    age=as.numeric(age)) %>%
  summary() 
gse10072 %>% filter(cancer!="y") %>% group_by(smok, metadata_sex) %>% count()
#smok  metadata_sex     n
#<chr> <chr>        <int>
#  1 FS    male            18
#2 NS    female          11
#3 NS    male             4
#4 S     female           4
#5 S     male            12