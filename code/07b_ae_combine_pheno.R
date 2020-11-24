# load all of the airway epithelium data
require('tidyverse')
require('GEOmetadb')
prev_data <- read_csv("data/sae_sl_mapped.csv")  # 10
#prev_studies
#[1] "GSE11784"  "GSE11952"  "GSE17905"  "GSE18385"  "GSE19407"  "GSE20257" 
#[7] "GSE63127"  "GSE19667"  "GSE64614"  "GSE108134"

# -- the rest of the pdata files -- #
new_studies <- setdiff(toupper(str_replace_all(list.files("data/pdata/"), ".csv", "")), "GSE97010") # 14
# remove: GSE97010 - acute cigaerette exposure
intersect(prev_studies, new_studies) # "GSE17905" "GSE18385" - two overlap :/
ae_studies <- union(prev_studies, new_studies)


new_studies_read <- lapply(new_studies, function(x) 
  read_csv(sprintf("data/pdata/%s.csv", tolower(x))))

list_cols <- unique(unlist(sapply(new_studies_read, function(x) colnames(x))))
my_cols <- setdiff(list_cols, c("fev1", "fev1_fvc", "diff_formatted", "inhaled_meds",
                                 "status", "copd_pheno", "ancestry",
                                 "hx_asthma", "rma_mas", "used_in_analysis"))

names(new_studies_read) <- new_studies
# -- fill in missing columns -- #
add_missing_cols <- function(sheet, study_name){
  missing_cols <- setdiff(all_cols, colnames(sheet))
  for (new_col in missing_cols){
    sheet[,new_col] <- NA
  }
  sheet$study <- study_name
  return(sheet %>% select(study, my_cols))
}

new_reform_sheets <- lapply(1:length(new_studies_read), function(i) 
  add_missing_cols(new_studies_read[[i]], names(new_studies_read)[[i]]))

# -- put together -- #
new_studies_comb <- do.call(rbind, new_reform_sheets)

# -- normalize all the columns -- #
# normalize tissue
new_studies_comb2 <- new_studies_comb %>% 
  mutate(tissue=case_when(
    !is.na(tissue) ~ tolower(tissue),
    is.na(tissue) ~ tolower(source_name))) %>%
  mutate(tissue=str_replace_all(tissue, "brushing|cells", ""),
         tissue=str_replace_all(tissue, "epithelial", "epithelium"),
         tissue=str_replace_all(tissue, "airways", "airway epithelium"),
         tissue=str_squish(tissue))
new_studies_comb2 %>% group_by(tissue) %>% count()
# remove the obtained from mskcc


prev_data1 <- prev_data %>% 
  mutate(tissue=tolower(source_name),
         tissue=case_when(
           str_detect(tissue, "small airway") ~ "small airway epithelium",
           str_detect(tissue, "large airway") ~ "large airway epithelium",
           str_detect(tissue, "airway") ~ "airway epithelium",
           str_detect(tissue, "trachea epithelial") ~ "trachea epithelium",
           TRUE ~ tissue
         ))
prev_data1 %>% group_by(tissue) %>% count()
# // TODO is trachea all trachea epithelium?

# combine ID cols
new_studies_comb3 <- new_studies_comb2 %>%
  unite("id", contains("id"), sep=";", na.rm=TRUE)

# history/cancer cols?
new_studies_comb4 <- new_studies_comb3 %>% 
  mutate(cancer=case_when(
    is.na(cancer) & str_detect(history, "without cancer") ~ "no",
    TRUE ~ cancer
  )) %>%
  select(-history) %>%
  # fix others
  mutate(sex=case_when(
    sex=="female" ~ "F",
    sex=="male" ~ "M",
    TRUE ~ sex
  )) %>%
  mutate(ethnic=ifelse(str_detect(ethnic, "hispnaic"), "hispanic", ethnic)) %>%
  rename(metadata_sex=sex, race_ethnicity=ethnic,
         pack_years=pkyrs)%>%
  mutate(age=as.numeric(age))

prev_data2 <- prev_data1 %>% 
  mutate(race_ethnicity=ifelse(race_ethnicity=="black/hispanic;black",
                               "black/hispanic", race_ethnicity)) %>%
  rename(id=dgm_id) 
prev_data2$cancer <- rep(NA, nrow(prev_data2))

prev_data3 <- prev_data2 %>% select(colnames(new_studies_comb4))
overlap_dat <- new_studies_comb4 %>% 
  semi_join(prev_data3, by=c("geo_accession"))
# over 700 :/

combine_vals <- function(x) paste(sort(unique(x[!is.na(x)])), collapse=";")
comb_data <- prev_data3 %>% 
  bind_rows(new_studies_comb4) %>%
  mutate(smok=ifelse(smok=="unknown", NA, smok))  %>%
  group_by(geo_accession) %>% 
  summarize(across(everything(),combine_vals)) 


# deal w conflicting data

# easy to resolve
comb_data2 <- comb_data %>% 
  group_by(geo_accession) %>%
  mutate(tissue=case_when( 
    str_detect(tissue, ";") & str_detect(tissue, "small airway") ~ "small airway epithelium",
    str_detect(tissue, ";") & str_detect(tissue, "large airway") ~ "large airway epithelium",
    TRUE ~ tissue),
  race_ethnicity=ifelse(race_ethnicity=="black;black/hispanic",
                          "black/hispanic", race_ethnicity),
  pack_years=ifelse(str_detect(pack_years, ";"), str_split(pack_years, ";")[[1]][[2]], pack_years)
  ) %>%
  ungroup() %>%
  mutate(metadata_sex=case_when(
    is.na(metadata_sex) ~ "",
    metadata_sex=="M" ~ "male",
    metadata_sex=="F" ~"female",
    TRUE ~ "")) %>%
  # replace "" with NA
  mutate(across(everything(),~ifelse(.=="", NA, .)))


comb_data2 %>% filter(str_detect(smok, ";")) %>% distinct(smok) # this is concerning


# what to do about the data that does not have smoking information?!

# --  add in expr_sex -- #
require('googlesheets4')
require('googledrive')
sl_sample <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1lEUVsyXLDcyUQB7mXk0B-LtsWkBn66gKaIV74OFAGVs/edit#gid=463642465")

comb_data_sl <- comb_data2 %>% 
  left_join(sl_sample %>% distinct(sample_acc, pred, expr_sex), 
            by=c("geo_accession"="sample_acc"))


stopifnot(nrow(comb_data_sl)==length(unique(comb_data_sl$geo_accession)))

summary(comb_data_sl %>%
          select(-source_name) %>%
          mutate(across(c(study, tissue, metadata_sex, expr_sex, race_ethnicity, smok, copd, cancer), as.factor)))
# we are missing 124 expr_sex :/, 52 smok, 28 smok are conflicting! :/ BOO

# remove data that doesn't go
comb_data_sl_filt <- comb_data_sl %>%
  filter(!is.na(smok), !is.na(expr_sex), !smok %in% c("NS;S"),
         copd == "n" | is.na(copd), cancer=="no" | is.na(cancer),
         tissue != "obtained from mskcc")

comb_data_sl_filt2 <- comb_data_sl_filt %>%
          select(-source_name, -copd, -cancer) %>%
          mutate(across(c(study, tissue, metadata_sex, expr_sex, 
                          race_ethnicity, smok), as.factor)) %>%
  mutate(across(c(age, pack_years), as.numeric))

summary(comb_data_sl_filt2)
nrow(comb_data_sl_filt2) #1209 samples, 97 are former smokers 
comb_data_sl_filt2 %>% filter(expr_sex!=metadata_sex) %>% nrow() # 10 mismatch
comb_data_sl_filt2 %>% write_csv("data/comb_sample_epithelium_metadata.csv")

# -- summarize each of the studies -- #
#  name, year, location, source, tissue, counts
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")

ae_studies_str <- paste(ae_studies, collapse="','")
gse_info <- dbGetQuery(con, sprintf("SELECT gse.gse, title, submission_date, gpl, pubmed_id, contributor, summary, overall_design 
FROM gse JOIN gse_gpl ON gse.gse=gse_gpl.gse
WHERE gse.gse in ('%s');", ae_studies_str))
gse_info %>% select(gse,  submission_date, gpl, pubmed_id)
# 23, each with own
table(gse_info$gpl) # almost all GPL570
gse_info %>% filter(gpl %in% c("GPL570", "GPL96")) %>% nrow() # 20

# list of cel files
gsm_str <- paste(comb_data_sl_filt2 %>% filter(smok!="FS") %>% pull(geo_accession), collapse="','")
cel_list <- dbGetQuery(con, sprintf("SELECT gsm, supplementary_file FROM gsm WHERE gsm in ('%s');", gsm_str))
head(cel_list) # 1087
# remove extra listed file
cel_list2 <- cel_list %>% separate_rows(supplementary_file, sep=";\\t") %>%
  filter(str_detect(supplementary_file, "CEL"))
cel_list2$supplementary_file[[1]]
cel_list2 %>% select(-gsm) %>% write_tsv("data/epithelium_list_cel_files.txt", 
                                         col_names=FALSE)

# gse to gsm + gpl
# separate into disc/valid
head(comb_data_sl_filt2)
studies_w_overlap <- comb_data_sl_filt2 %>% 
  distinct(study) %>%
  mutate(gse=study) %>%
  separate_rows(gse, sep=";") %>%
  left_join(gse_info %>% select(gse, submission_date, gpl)) %>%
  group_by(study) %>%
  mutate(num_studies=n()) %>%
  ungroup()

gse_info %>% group_by(gpl) %>% count()
# GPL9185 - mouse
# GPL13243 - [HuGene10stv1_Hs_ENSG] Affymetrix GeneChip Human Gene 1.0 ST Array [Brainarray Version 11.0.1]
gse_info %>% filter(gpl=="GPL13243") %>% pull(gse) # "GSE37147" 

# before I used samples assigned to "GSE11784" as validation
comb_data_sl_filt3 <- comb_data_sl_filt2 %>% filter(smok!="FS") 
comb_data_sl_filt3 %>% filter(str_detect(study, "GSE11784")) %>% nrow() #123
comb_data_sl_filt3 %>% filter(!str_detect(study, "GSE11784")) %>% nrow() #964



study_to_platform <- studies_w_overlap %>% 
  distinct(study, num_studies, gpl) %>%
  filter(gpl!="GPL9185") %>% # remove mouse platform one
  group_by(study) %>%
  summarize(gpl=paste(unique(gpl), collapse=";"), 
            num_studies=unique(num_studies)) %>%
  ungroup()


# count the number of former smokers
comb_data_sl_filt2 %>% filter(smok=="FS") %>% group_by(study) %>% count()
# GSE37147  <-- has the largest number of former smokers
#   this study is also a different platform! good to add to validation
# GSE7895 has a number, 
# GSE19027 has 8

# count the number of copd patients
comb_data_sl %>% filter(copd=="y") %>% group_by(study) %>% count() %>%
  arrange(desc(n))
# GSE108134  has the most
# GSE37147 has a large number (also a good one for the alternate)
# GSE64614 has 22

# which studies are from what/a different lab?

# this is the only one that is a different lab... Woodruff lab UCSF
# rest are BU/Cornell
comb_data_sl %>% filter(study=="GSE4302")  %>% 
  group_by(expr_sex, smok) %>%
  count() # only two female smokers

# which was used for validation previously?
comb_data_sl %>% filter(study=="GSE7895") %>%
  group_by(expr_sex, smok) %>%
  count()#?
comb_data_sl %>% filter(study=="GSE37147") %>%
# no non-smokers...

comb_data_sl %>% filter(study %in% c("GSE37147", "GSE4302")) %>% nrow()
# 313


# do we have the same set of studies included total too?
# yang et al studies:
yang_disc_studies <- c( "GSE19667", "GSE5058", "GSE11784", "GSE11906", 
                        "GSE8545", "GSE11952", "GSE63127", "GSE13933", 
                        "GSE4498", "GSE18385","GSE17905", "GSE13931", 
                        "GSE19407") # 13
# 22 airway studies, 9 overlap
intersect(yang_disc_studies, ae_studies ) 
setdiff(yang_disc_studies, ae_studies )  # 4
# why don't we get these? "GSE5058"  "GSE8545"  "GSE13933" "GSE4498" 
# they're all included in the "manual annotation" file but not the final list
# they also have some missing expression sex labels
# were they removed b/c overlap? unclear to me...


setdiff(ae_studies, yang_disc_studies) # 13 --> 11 (if remove replication)
# they required demographic info:  sex, smoking status, age, COPD status, ethnicity and pack-years 
# were available for 211 subjects 
# after removing duplicate data and outliers


# now put together
sl_data_ds <- comb_data_sl_filt3 %>% 
  left_join(study_to_platform) %>%
  mutate(ds=case_when(
    study %in% c("GSE37147", "GSE11784" , "GSE4302", "GSE7895" ) ~ "validation",
    gpl == "GPL570" ~ "discovery_gpl570",
    gpl == "GPL96" ~ "discovery_gpl96"
    
  )) %>%
  filter(!is.na(gpl)) # remove GSE43949 which is mouse actually

sl_data_ds %>% group_by(ds) %>% count()


sl_data_ds %>% write_csv("data/ae_samples_studies_ds.csv")
