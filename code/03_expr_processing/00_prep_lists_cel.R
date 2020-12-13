library('GEOmetadb')
library('tidyverse')
# put together the lists for lung, blood, and small datasets
small_ds <- read_csv("data/download_info_small_n.csv")
lb_ds <- read_csv("data/download_info_lung_blood.csv")
all_ds <- small_ds %>% bind_rows(lb_ds) # 3247
all_ds2 <- all_ds %>%
  separate_rows(supplementary_file, sep=";\t") %>%
  filter(str_detect(supplementary_file, "CEL") |
           str_detect(supplementary_file, "cel")) # hmm this does not work

# there are ton of missing supplementary files. shit
table(is.na(all_ds$supplementary_file))
#FALSE  TRUE 
#1775  1472 

# add in the ones that didn't have a CEL file to begin with
no_cel <- all_ds %>% 
  filter(!str_detect(supplementary_file, "CEL|cel"))

all_ds3 <- all_ds2 %>% bind_rows(no_cel)
dim(all_ds3) # 1775

# ------ get the data for the samples with missing supplementary files ----- #
missing_supp <- all_ds %>% filter(is.na(supplementary_file)) # 1472 samples, 12 studies

missing_supp %>% group_by(gse, gpl) %>% count() %>% arrange(desc(n)) 
missing_supp %>% distinct(gse, gpl) %>% group_by(gpl) %>% count() %>% arrange(desc(n)) 

# GPL10904
# GPL10399
# GPL6102
# GPL6884
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")

# doesn't help :/
missing_samples_str <- paste(missing_supp %>% pull(gsm), collapse="','")
missing_samples_supp <- dbGetQuery(con, sprintf("SELECT gsm, supplementary_file FROM gsm WHERE gsm IN ('%s');",
                                                missing_samples_str))
table(missing_samples_supp$supplementary_file=="")

# does help!
missing_studies <- missing_supp %>% distinct(gse) %>% pull(gse)
missing_studies_str <- paste(missing_studies, collapse="','")
missing_studies_supp <- dbGetQuery(con, sprintf("SELECT gse, supplementary_file FROM gse WHERE gse IN ('%s');",
                                                missing_studies_str))
table(missing_studies_supp$supplementary_file=="")


# ---- get GPL info ---- #
list_gpls <- all_ds %>% distinct(gpl) %>% pull(gpl) # 14
platform_str <- paste(list_gpls, collapse="','")
platform_info <- dbGetQuery(con, sprintf("SELECT gpl, title, manufacturer,data_row_count, bioc_package, supplementary_file FROM gpl WHERE gpl IN ('%s');",
                        platform_str)) %>%
  as_tibble()
dbDisconnect(con)

platform_info %>% write_csv("data/platform_info.csv")
# --------- #
# group into affy vs illumina
affy_platforms <- platform_info %>% filter(manufacturer=="Affymetrix") # 9
illumina_platforms <- platform_info %>% filter(manufacturer=="Illumina Inc.") # 4
agilent_platforms <- platform_info %>% filter(manufacturer=="Agilent Technologies") # 1

affy_data <- all_ds3 %>% semi_join(affy_platforms, by="gpl") # 1739
agilent_data <- all_ds3 %>% semi_join(agilent_platforms, by="gpl") # 36

# --- group affy further by platform --- #
affy_data %>% group_by(gpl) %>% count() %>% arrange(desc(n))


# "GPL16384" <-- remove this b/c miRNA expression
# for missing: 
#  GPL16686: "hugene20sttranscriptcluster.db"
#  brainarray one? GSE62699 -- unclear if can use plus2 or not

# all of the illumina data does not have listed CEL files
# -- maybe this is all in the RAW.tar?
illumina_data <- all_ds3 %>% semi_join(illumina_platforms, by="gpl") # 0
missing_supp %>% 
  distinct(gpl) %>% 
  left_join(platform_info %>% 
              select(gpl, manufacturer, title))
illumina_studies_raw <- missing_studies_supp %>% 
  left_join(missing_supp %>% distinct(gse, gpl)) %>%
  left_join(platform_info %>% 
              select(gpl, manufacturer)) %>%
  separate_rows(supplementary_file, sep=";\t") 
illumina_studies_raw2 <- illumina_studies_raw %>% 
  select(manufacturer, gpl, gse, supplementary_file) %>%
  arrange(gpl, gse)

# set up for download

# manufacturer, platform, study, samples
samp_files_raw <- all_ds3 %>% left_join(platform_info %>% 
                        select(gpl, manufacturer), by="gpl")  %>%
  select(manufacturer, gpl, gse, supplementary_file) %>%
  arrange(gpl, gse)

list_files <- samp_files_raw %>% bind_rows(illumina_studies_raw)
list_files %>% mutate(manufacturer=ifelse(manufacturer=="Illumina Inc.", "Illumina",
                                          manufacturer)) %>% write_csv("data/list_f_to_download.csv")
