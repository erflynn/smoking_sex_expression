# 01_get_list_studies.R
# 11/8/2020
# Updated code for getting the full counts and list of studies
# from GEO, SRA, and ArrayExpress.
#  - combines code from "02_new_studies.R" and "smoking_list.R"
# GOAL:
#  - get counts in the process
#  - done at the end of this!

# lists --> deduplications --> organism breakdown
# --> add in old
# ??? filter based on availability???

library('tidyverse')
library('GEOmetadb')

SMOK.STR <- "smok|nicotine|tobacco|cigarette"
RB.PATH <- "../drug_trt/data/"

# 0. Load the RB data, we are filtering for this
rb_meta_m <- read_csv(sprintf("%s/01_sample_lists/rb_metadata/human_microarray_experiment_metadata.csv", RB.PATH),
                      col_types="cccc") 
h_rb_metadata <- read_csv(sprintf("%s/01_sample_lists/rb_metadata/human_rnaseq_experiment_metadata.csv", RB.PATH) )
h_rb_metadata2 <- h_rb_metadata %>%
  bind_rows(rb_meta_m %>% anti_join(h_rb_metadata, by="study_acc") )
h_rb_metadata2 %>% filter(str_detect(study_acc, "GSE")) %>% nrow() # 14,752 studies
h_rb_metadata2 %>% filter(str_detect(study_acc, "DRP|SRP|ERP")) %>% nrow() # 8,357 studies


# 1. search GEO using GEOmetadb
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite") # 11/8/2020
gse_dat <- dbGetQuery(con, "SELECT gse, title, summary, overall_design FROM gse;")
gse_smok <- gse_dat %>% 
  as_tibble() %>%
  group_by(gse) %>% 
  mutate(str=paste(c(title, overall_design, summary), collapse=" ")) %>%
  filter(grepl(SMOK.STR, str, ignore.case=TRUE))
# --> 757. a lot

gse_smok_list <- paste(gse_smok$gse, collapse=c("\',\'"))


gsm_dat <- dbGetQuery(con, "SELECT gsm, molecule_ch1, title, source_name_ch1,
                      characteristics_ch1, treatment_protocol_ch1, description, 
                      channel_count, organism_ch1 FROM gsm 
                      WHERE organism_ch1 IN ('Homo sapiens');")
unique(gsm_dat$molecule_ch1)


gsm_dat2 <- gsm_dat %>% 
  filter(channel_count==1, molecule_ch1 %in% c("total RNA", "polyA RNA")) %>%
  as_tibble()


# extremely SLOW - abt 10min
gsm_dat_smok <- gsm_dat2 %>%
  group_by(gsm) %>% 
  mutate(str=paste(c(title, source_name_ch1, characteristics_ch1, 
                     treatment_protocol_ch1, description), collapse=" ")) %>%
  filter(grepl(SMOK.STR, str, ignore.case=TRUE))

length(unique(gsm_dat_smok$gsm)) # 36,267

gsm_smok_list <- paste(gsm_dat_smok$gsm, collapse=c("\',\'"))

# A) studies with mentions in the study
gse_gsm <- dbGetQuery(con, sprintf("SELECT gse,gsm FROM gse_gsm WHERE gse IN ('%s');", gse_smok_list))
gsm_dat_smok_gse <- gsm_dat2 %>% inner_join(gse_gsm, by="gsm")
# 373 studies, 24,319 samples
gsm_dat_smok_gse2 <- gsm_dat_smok_gse %>% semi_join(h_rb_metadata2, by=c("gse"="study_acc"))
# 185 studies, 13,256 samples

gsm_gse0 <- gsm_dat_smok_gse2 %>% distinct(gse, gsm)

# B) studies with mentions in the samples
gsm_gse <- dbGetQuery(con, sprintf("SELECT gse, gsm FROM gse_gsm WHERE gsm IN ('%s');", gsm_smok_list)) 
gsm_dat_smok2 <- gsm_dat_smok %>% inner_join(gsm_gse, by="gsm") %>%
  semi_join(h_rb_metadata2, by=c("gse"="study_acc"))
# 481 studies, 35,752 samples -->
# 252 studies, 19,547 samples

#  i) the study/sample pair does not have associated study metadata
non_overlapping <- gsm_dat_smok2 %>% anti_join(gsm_gse0, by=c("gse", "gsm"))  
# 223 studies, 22,418 samples -->
# 113 studies, 11,900 samples

#  ii) overlapping - included in both
gse_and_gsm_smok_dat <- gsm_dat_smok_gse2 %>% semi_join(gsm_gse, by="gsm")
# 258 studies, 14,494 samples -->
# 148 studies, 8,644

#  iii) only in the study but not in the sample?
gse_dat_smok_no_gsm_dat <- gsm_dat_smok_gse2 %>% anti_join(gsm_gse, by="gsm") 
# 146 studies, 9,825 samples -->
# 62 studies, 4,821 samples


# C) put it all together
gsm_gse2 <- gsm_gse0 %>% bind_rows(gsm_gse) %>% distinct(gse, gsm) %>%
  semi_join(h_rb_metadata2, by=c("gse"="study_acc"))
# 596 studies, 45,577 samples -->
# 298 studies, 24,368 samples

# write this out
gsm_gse2 %>% write_csv("data/gse_gsm_smok.csv")
save(gsm_dat_smok_gse2, gsm_dat_smok2, file="data/geo_smok_attr_pair.RData")
# add in additional information

# get GPL data
gpl_data <-dbGetQuery(con, sprintf("SELECT 
  gse.gse, gse.title, gse.summary, gse.overall_design, organism, 
  gse.type, technology, gse_gpl.gpl, gpl.title AS platform FROM gse 
  JOIN gse_gpl ON gse.gse=gse_gpl.gse 
  JOIN gpl ON gse_gpl.gpl=gpl.gpl
  WHERE gse.gse IN (\'%s\');", paste(unique(gsm_gse2$gse), collapse=c("\',\'"))))
# "organism", "type", "technology"

gpl_filt <- gpl_data %>% 
  select(gse, gpl, organism, type, technology, platform) %>%
  filter(str_detect(organism, "Homo sapiens") | str_detect(organism,"Mus musculus")) %>%
  as_tibble() %>%
  separate_rows(type, sep=";\\t") %>%
  separate_rows(technology, sep=";\\t") 
gpl_filt2 <- gpl_filt %>% 
  filter(type %in% c("Expression profiling by array", "Expression profiling by high throughput sequencing")) 
gpls_only <- gpl_filt2 %>% distinct(gpl, platform, technology)

# remove methylation + microRNA
gpls_part2 <- gpls_only %>% 
  filter(!grepl("methylation", platform, ignore.case=T)) %>%
  filter(!grepl("miR|MiR|microRNA", platform, ignore.case=F))

gpl_data2 <- gpl_data %>% select(gse, gpl) %>% semi_join(gpls_part2, by="gpl")
length(unique(gpl_data2$gpl))

gsm_gse3 <- gsm_gse2 %>% inner_join(gpl_data2) #doesn't remove any
length(unique(gsm_gse3$gse))
length(unique(gsm_gse3$gsm))

sample_counts <- gsm_gse3 %>% 
  distinct(gse, gsm) %>% 
  group_by(gse) %>% 
  summarise(num_samples=n())

geo_smoking_studies <- gse_dat %>% 
  semi_join(gsm_gse3) %>% 
  left_join(gpl_data2 %>% group_by(gse) %>% summarise(gpl=paste(gpl, collapse=";"))) %>%
  left_join(gpl_data %>% distinct(gse, organism) %>% group_by(gse) %>% summarise(organism=paste(organism, collapse=";"))) %>%
  left_join(sample_counts)  %>%
  as_tibble()
stopifnot(nrow(geo_smoking_studies)==length(unique(geo_smoking_studies$gse)))

# TODO: add back in type, add in a date column!
gse_info <- dbGetQuery(con, sprintf("SELECT gse, type, submission_date, pubmed_id, contributor FROM gse
           WHERE gse IN ('%s')", paste(unique(geo_smoking_studies$gse), collapse=c("\',\'"))))
geo_smoking_studies2 <- geo_smoking_studies %>% 
  left_join(gse_info %>% select(gse, type, submission_date, pubmed_id)) %>%
  dplyr::rename(study_acc=gse)

geo_smoking_studies2 %>% select(-type, -organism) %>% write_csv("data/geo_smok_studies_1231.csv")

# 2. search ArrayExpress using the browser and download
# search: "smoking OR nicotine OR smoker OR cigarette OR smoke"
ae_experiments <- read_tsv("data/ArrayExpress-Experiments-201110-192305.txt") # 544
filtered_ae <- ae_experiments %>% 
  filter(str_detect(Type, "transcription profiling by array") |
           str_detect(Type, "RNA-seq of coding RNA"))  %>% # 449
  filter(str_detect(Organism, "Homo sapiens")) %>%  # 345
  filter(!grepl("miR|MiR|microRNA|MicroRNA", Title, ignore.case=F)) # 312

filtered_ae2 <- filtered_ae %>% 
  mutate(Accession=str_replace_all(Accession, "E-GEOD-", "GSE")) %>%
  semi_join(h_rb_metadata2, by=c("Accession"="study_acc")) # --> 176

ae_gses <- filtered_ae2 %>% filter(str_detect(Accession, "GSE")) %>% pull(Accession) # 238 --> 171
length(intersect(ae_gses, geo_smoking_studies$gse)) # 170 intersect
length(setdiff(ae_gses, geo_smoking_studies$gse)) # 1 new

array_exp_studies <- filtered_ae2 %>% 
  anti_join(geo_smoking_studies, by=c("Accession"="gse")) %>% 
  dplyr::rename("study_acc"="Accession", "title"="Title", "organism"="Organism", 
         "num_samples"="Assays", "type"="Type") %>%
  select(-`Present in Atlas`, -`ArrayExpress URL`)

length(unique(array_exp_studies$study_acc)) # 105 studies --> 6 studies
sum(array_exp_studies$num_samples) # 57097 samples --> 189 samples
stopifnot(nrow(array_exp_studies)==length(unique(array_exp_studies$study_acc)))
array_exp_studies %>% write_csv("data/array_exp_studies_1231.csv")

# additional fields: description?

# ------ 3. search SRA ----- #


# A) search MetaSRA
# MetaSRA for human + mouse
# 821 samples from 24 studies annotated with "smoking behavior" | "nicotine dependence"
sra_human <- read_csv("data/metaSRA-human_samples.csv")

metasra_studies_h <- sra_human %>% 
  group_by(study_id, study_title) %>% 
  summarise(num_samples=n()) %>%
  mutate(organism="Homo sapiens") # 24

metasra_studies <- metasra_studies_h  %>%
  dplyr::rename(study_acc=study_id, title=study_title) %>% 
  semi_join(h_rb_metadata2, by="study_acc") %>%
  ungroup() # 17

# B) Search refine-bio metadata for study title/description
rnaseq_rb2 <- h_rb_metadata2 %>% filter(str_detect(study_acc, "DRP|SRP|ERP")) # 8,357 studies
rnaseq_rb_smok <- rnaseq_rb2 %>%
  filter(grepl(SMOK.STR, title, ignore.case = T) |
           grepl(SMOK.STR, description, ignore.case = T)) # 43



# C) Look at RNA-seq covariate data
load(sprintf("%s/data_old/01_sample_lists/all_sample_attrib_clean.RData", RB.PATH) ) # --> all_attrib_clean

rb_exp_sample_map <- read_csv(sprintf("%s/01_sample_lists/rb_metadata/human_microarray_exp_to_sample.csv", RB.PATH)) %>%
  bind_rows(read_csv(sprintf("%s/01_sample_lists/rb_metadata/human_rnaseq_exp_to_sample.csv", RB.PATH)))


rb_smok_attr <- all_attrib_clean %>% 
  filter(str_detect(key_clean, SMOK.STR) |  
           str_detect(value_clean, SMOK.STR)) %>%
  semi_join(rb_exp_sample_map, by="sample_acc") 
length(unique(rb_smok_attr$sample_acc)) # 12186 samples
rb_smok_attr_rnaseq <- rb_smok_attr %>%
  filter(str_detect(sample_acc, "ERR|SRR|DRR"))  %>%
  semi_join(rb_exp_sample_map, by="sample_acc") 
length(unique(rb_smok_attr_rnaseq$sample_acc)) # 1386 samples


ae_study_sample_acc <- all_attrib_clean %>%
  distinct(sample_acc) %>%
  filter(str_detect(sample_acc, "^E-")) %>%
  group_by(sample_acc) %>%
  mutate(study_acc=str_split(sample_acc, "_")[[1]][[1]]) %>%
  ungroup() %>%
  select(study_acc, sample_acc)

sample_map <-  rb_exp_sample_map %>% 
  bind_rows(ae_study_sample_acc)
rb_exp_sample <-sample_map %>%
  semi_join(rb_smok_attr, by="sample_acc")
rb_exp <- rb_exp_sample %>% 
  group_by(study_acc) %>% 
  summarise(num_samples=n()) # --> 214 studies

sra_rb_exp <- rb_exp %>% 
  filter(str_detect(study_acc, "ERP|SRP|DRP")) # 36
# --> adds one AE study

sra_study_acc <- metasra_studies %>% 
  filter(organism=="Homo sapiens") %>% select(study_acc) %>%
  bind_rows(sra_rb_exp %>% select(study_acc)) %>%
  bind_rows(rnaseq_rb_smok %>% select(study_acc)) %>%
  distinct(study_acc) # 70

# put together a table
rb_metadata <- h_rb_metadata2 %>%
  mutate(organism="Homo sapiens") %>%
  select(-date)
length(setdiff(sra_study_acc$study_acc, rb_metadata$study_acc) ) 
metasra_studies %>% anti_join(rb_metadata, by="study_acc") 
sra_metadata_filt <- rb_metadata %>% semi_join(sra_study_acc)
sra_sample_counts <- rb_exp_sample_map %>% 
  semi_join(sra_study_acc) %>% 
  group_by(study_acc) %>%
  summarise(num_samples=n())
sra_metadata <- sra_metadata_filt %>% 
  left_join(sra_sample_counts, by=c("study_acc"))
sra_metadata2 <- sra_metadata %>% bind_rows(
  metasra_studies %>% anti_join(rb_metadata, by="study_acc") %>% 
    mutate(description="") %>%
    select(colnames(sra_metadata))
)
sum(sra_metadata2$num_samples) # 6148 samples, 70 studies

# additional fields:
# - submission_date
# - type
# - platform

sra_metadata %>% write_csv("data/sra_studies_1231.csv")

# ---- 4. Put the data together ----- #
combined_data <- geo_smoking_studies2 %>% 
  group_by(study_acc) %>%
  mutate(description=paste(c(summary, overall_design), collapse=";")) %>%
  ungroup() %>%
  select(colnames(sra_metadata)) %>%
  mutate(source="GEO") %>%
  bind_rows(array_exp_studies %>%
              mutate(description="") %>%
              select(colnames(sra_metadata)) %>%
              mutate(source="ArrayExpress")) %>%
  bind_rows(sra_metadata %>%
              mutate(source="SRA")) %>%
  filter(str_detect(organism, "Homo sapiens")) %>%
  select(-organism)

nrow(combined_data) # 367
(study_counts <- combined_data %>%
  group_by(source) %>%
  summarise(num_studies=n()) %>%
  pivot_wider(names_from=source, values_from=num_studies, values_fill=0))
study_counts %>% write_csv("data/study_counts_before_filtering_1231.csv")

(sample_counts <- combined_data %>%
  group_by(source) %>%
  summarise(tot_samples=sum(num_samples)) %>%
  pivot_wider(names_from=source, values_from=tot_samples, values_fill=0))
# TODO this includes duplicates!!
sample_counts %>% write_csv("data/sample_counts_before_filtering_1231.csv")

# ---- 5. Add in previous notes ----- #

prev_manual_annot <- read_csv("data/smok_data_manual_annot_0408.csv")  %>% dplyr::rename(study_acc=gse) # 274
prev_manual_annot %>% anti_join(combined_data, by="study_acc") # lost 46... :/
combined_data %>% semi_join(prev_manual_annot, by="study_acc")  %>% nrow() # 228
combined_data %>% anti_join(prev_manual_annot, by="study_acc")  %>% nrow() # 139

combined_data_w_annot <- combined_data %>% 
  left_join(prev_manual_annot %>% 
              select(study_acc, keep, design, type, treatment, tissue2), by="study_acc")
combined_data_w_annot %>% write_csv("data/manual_annot_1231.csv")


# ---- 6. Annotate if there is available sample metadata ----- #

# ---- 7. What data is actually available? ---- #

####### SANITY CHECKING SECTION #######
# STEPS TO BE DONE WITH THIS:
#  [x] SRAdb search - prefer not? I don't think this will help + missingness
#  [x] rb metadata -- very little added --> NOPE
#  [x] ENA biosamples? GEO biosamples?  --> NOPE
#  - get additional data attributes from ENA
#      pro: RB is incomplete, also we could look at "source_name", "title", etc, would give us ALL the info
#      con: I have to download these...

# best way to do this -- use ENA and crawl all available info
# TODO - first double check if we *ARE* looking at all studies vs only present studies
# ... we're looking at all runs FROM refine-bio
# -- BioSample search? NOPE --
# public + in SRA
#biosample <- vroom("~/Downloads/biosample_result.txt", delim="\n", col_names="bios_list") # GROSS
# 
# # ---- SRAdb ---- #
# # we know is missing info, do this server-side
# 
# library('SRAdb')
# sra_con <- dbConnect(SQLite(), "SRAmetadb.sqlite")
# sra_studies <- dbGetQuery(sra_con, "SELECT study_accession, study_title, study_type, 
# study_abstract, study_attribute, submission_accession FROM study;")
# sra_studies2 <- sra_studies %>% filter(is.na(study_type) |
#                                          study_type %in% c("Transcriptome Seqencing", "Transcriptome Analysis", "Other"))
# 
# sra_smok_studies <- sra_studies2 %>%
#   filter(grepl(SMOK.STR, study_title, ignore.case = T) |
#            grepl(SMOK.STR, study_abstract, ignore.case = T) |
#            grepl(SMOK.STR, study_attribute, ignore.case = T)) # 370
# 
# sra_samples <- dbGetQuery(sra_con, "SELECT sample_accession, taxon_id, description, 
#                           sample_attribute, submission_accession FROM sample;")
# sra_experiments <- dbGetQuery(sra_con, "SELECT experiment_accession, title, 
#                               design_description, platform, experiment_attribute, 
#                               submission_accession FROM experiment;")
# sra_runs <- dbGetQuery(sra_con, "SELECT run_accession, run_attribute, 
#                        submission_accession FROM run;")
# 
# 
# 
# 
# # ----- check: does rb sample data add any? NOPE ----- #
# rb_smok_rnaseq_sample <- read_csv("../drug_trt/data/01_sample_lists/rb_metadata/human_rnaseq_sample_metadata.csv",
#                                   col_types="ccccccccccc") %>%
#   bind_rows(read_csv("../drug_trt/data/01_sample_lists/rb_metadata/mouse_rnaseq_sample_metadata.csv", col_types="ccccccccccc")) %>%
#   select(acc, title, platform, compound, trt) %>%
#   rename(sample_acc=acc) %>%
#   group_by(sample_acc) %>%
#   mutate(description=paste(c(title, compound, trt), collapse=";")) %>%
#   filter(grepl(SMOK.STR, description, ignore.case = T) ) # 126
# 
# rb_rnaseq_sample2 <- rb_exp_sample_map %>% 
#   semi_join(rb_smok_rnaseq_sample, by="sample_acc") %>% 
#   distinct(study_acc) # 9
# rb_rnaseq_sample2 %>% 
#   anti_join(metasra_studies)  %>% # 8
#   anti_join(rb_from_acc) %>% # 8
#   anti_join(h_rb_smok %>% bind_rows(m_rb_smok)) # no new studies
# 
# # ----- check: does refine-bio microarray add any? ------ #
# h_rb_m <- read_csv("../drug_trt/data/01_sample_lists/rb_metadata/human_microarray_experiment_metadata.csv", col_types = "cccc")
# m_rb_m <- read_csv("../drug_trt/data/01_sample_lists/rb_metadata/mouse_microarray_experiment_metadata.csv", col_types = "cccc")
# h_rb_smok_m <- h_rb_m %>%
#   filter(grepl(SMOK.STR, title, ignore.case = T) |
#            grepl(SMOK.STR, description, ignore.case = T)) # 192
# m_rb_smok_m <- m_rb_m %>%
#   filter(grepl(SMOK.STR, title, ignore.case = T) |
#            grepl(SMOK.STR, description, ignore.case = T)) # 54
# 
# # 1 mouse AE study
# missing_micro <- h_rb_smok_m %>% bind_rows(m_rb_smok_m) %>% # 246
#   anti_join(geo_smoking_studies2, by="study_acc") %>%  # 62
#   anti_join(array_exp_studies, by="study_acc") %>%  # 59
#   filter(!str_detect(study_acc, "ERP|SRP|DRP")) # 1 study
# 
# # 7 SRP studies missing from RNA-seq
# h_rb_smok_m %>% bind_rows(m_rb_smok_m) %>% 
#   filter(str_detect(study_acc, "ERP|SRP|DRP")) %>%
#   anti_join(m_rb_smok %>% bind_rows(h_rb_smok), by="study_acc")
# 
# 
# # there are a ton of studies that just have smoking as a covariate -- how do I ID them?
# # e.g. E-MTAB-1708 -- though this doesn't have smoking status labels so it is unclear how I'd deal w this
