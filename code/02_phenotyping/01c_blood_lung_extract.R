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


blood_studies <- read_tsv("data/list_lung_blood.tsv", col_names=F)
blood_studies_str <- paste(blood_studies$X1, collapse="','")

# connect and grab all the GSE info
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
gse_info <- dbGetQuery(con, sprintf("SELECT gse.gse, title, 
submission_date, gpl, pubmed_id, contributor, summary, overall_design 
FROM gse JOIN gse_gpl ON gse.gse=gse_gpl.gse
WHERE gse.gse IN ('%s');", blood_studies_str))


# grab all the GSM info
# the characteristics_ch1 column is what has the covariate info
gsm_info <- dbGetQuery(con,  sprintf("SELECT gsm.gsm, gse_gsm.gse, title, source_name_ch1,
description, gpl, channel_count, organism_ch1, molecule_ch1, characteristics_ch1,
supplementary_file
FROM gsm JOIN gse_gsm ON
                       gse_gsm.gsm=gsm.gsm WHERE gse_gsm.gse IN ('%s');", 
                       blood_studies_str))

dbDisconnect(con)

# all human, total RNA, 1 channel
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

    
gse_w_info <- smok2 %>% group_split(gse) # --> 26 studies

# clean up the phenotype data
get_pheno_table <- function(df){
  df %>% 
  separate_rows(characteristics_ch1, sep=";\t") %>%
  separate(characteristics_ch1, into=c("attr", "value"), sep=":") %>%
  mutate(value=str_squish(value)) %>%
  pivot_wider(names_from=attr, values_from=value)
}

phe_tabs <- lapply(gse_w_info, get_pheno_table) 
names(phe_tabs) <- unique(smok2$gse)


# look at one particular study
gse10072 <- phe_tabs["GSE10072"][[1]] %>%
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