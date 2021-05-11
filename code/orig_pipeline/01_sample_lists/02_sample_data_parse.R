# 02_sample_data_parse.R 
# Clean up the data and parse the sample metadata --> get a "rough" phenotype data type
# Use this in a couple different ways
#

library(tidyverse)
RB.PATH <- "../drug_trt/data/"

load("data/geo_smok_attr_pair.RData") # --> gsm_dat_smok_gse2 (gsm data for GSEs with mentions), gsm_dat_smok2 (gsm data with mentions)
load(sprintf("%s/data_old/01_sample_lists/all_sample_attrib_clean.RData", RB.PATH) ) # --> all_attrib_clean

SMOK.STR <- "smok|nicotine|tobacco|cigarette"
SEX.STR <- "sex|gender|male"
RACE_ETH.STR <- "race|ethnic|ancestry"
AGE.STR <- "\\bage"
PACK_YEAR.STR <- "pack|pkyr|pyr"

# helpful function for condensing/deduplication
condense_dedup <- function(x) paste(unique(str_split(x,";")[[1]]), collapse=";")

find_col <- function(df1, df2, SEARCH.STR){
  # find a column in the long or wide version of a table
  name_col <- colnames(df2)[str_detect(colnames(df2), SEARCH.STR)]
  val_col <- df1 %>% filter(str_detect(value, SEARCH.STR)) %>% distinct(attr) %>% pull(attr)
  my_cols = unique(c(name_col, val_col))
  return(my_cols)
}

rename_col <- function(df2, my_cols, cols_str){
  # rename the found column to the col_str
  if (length(my_cols)==0){
    return(df2 %>% mutate({{ cols_str }} := "-999"))
  }
  if (length(my_cols)==1){
    if (cols_str %in% colnames(df2)){
      return(df2)
    } else {
      df4 <- df2 %>% 
        dplyr::rename({{ cols_str }} := my_cols)
    }
  }
  if (length(my_cols)>1){
    if (cols_str %in% colnames(df2)){
      print(" reduced the columns")
      return(df2)
    }
    non_generic <- setdiff(my_cols, c("description", "source_name_ch1", "title"))
    if (length(non_generic)==1){
      df4 <- df2 %>% 
        dplyr::rename({{ cols_str }} := non_generic[[1]])
    } else {
      if (length(non_generic) > 1){ 
        my_cols <- non_generic
      }
      # add suffixes
      df3 <- df2
      for(i in 1:length(my_cols)){
        new_str <- paste(cols_str, i, sep="_");
        df3 <- df3 %>% 
          dplyr::rename( {{ new_str }} :=  my_cols[[i]] ) ;
      }
      
      # put them together
      df4 <- df3 %>%
        unite(new_col, contains(sprintf("%s_", cols_str)), sep=";", na.rm=TRUE) %>%
        dplyr::rename({{cols_str}} := new_col)  
      
    }
  }
  return(cbind(df2, df4 %>% select({{cols_str}})))
}


find_col_stats <- function(df){
  # find the columns with the covariate info and rename them
  df1 <- df %>%
    mutate(attr=tolower(attr),
           value=tolower(value)) 
  df2 <- df1 %>%
    unique() %>%
    pivot_wider(names_from=attr, values_from=value, values_fn=unlist)#, values_fill="-999")
  smok_cols <- find_col(df1, df2, SMOK.STR)
  sex_cols <- find_col(df1, df2, SEX.STR)
  re_cols <- find_col(df1, df2, RACE_ETH.STR)
  age_cols <- find_col(df1, df2, AGE.STR)
  pkyrs_cols <- find_col(df1, df2, PACK_YEAR.STR)
  copd_cols <- find_col(df1, df2, "copd")
  
  #tissue_cols <- find_col(df1, df2, "tissue|organ")
  #cell_cols <- find_col(df1, df2, "\\bcell")
  #trt_cols <- find_col(df1, df2, "treatment|compound|drug")
  #time_cols <- find_col(df1, df2, "\\btime|hours|days")
  col_list <- list("smok"=smok_cols, 
                   "sex"=sex_cols, 
                   "race_ethnicity"=re_cols, 
                   "age"=age_cols, 
                   "pkyrs"=pkyrs_cols, 
                   "copd"=copd_cols)
  df2.1 <- rename_col(df2, smok_cols, "smok")
  df2.2 <- rename_col(df2.1, sex_cols, "sex")
  df2.3 <- rename_col(df2.2, re_cols, "race_ethnicity")
  df2.4 <- rename_col(df2.3, age_cols, "age")
  df2.5 <- rename_col(df2.4, pkyrs_cols, "pkyrs")
  df2.6 <- rename_col(df2.5, copd_cols, "copd")
  
  
  df4 <- df2.6 %>% group_by(sample_acc) %>% mutate(across(names(col_list), condense_dedup))
  
  # rearrange the data frame
  # --> these labels + any extra columns
  return(list("df"=df2, 
              "df_clean"= df4 %>% select(study_acc, sample_acc, names(col_list), everything()),
              "cols"=col_list))
}


# --- apply to the GEO data ---  #
geo_gsm_dat <- gsm_dat_smok2 %>% 
  select(-str) %>%
  bind_rows(gsm_dat_smok_gse2 %>% anti_join(gsm_dat_smok2, by=c("gse", "gsm"))) %>%
  select(-molecule_ch1, -organism_ch1, -channel_count) %>%
  select(gse, everything()) %>%
  ungroup() %>%
  unique()
length(unique(geo_gsm_dat$gsm)) # 24,368

split_by_gse <- geo_gsm_dat %>% 
  pivot_longer(title:description, names_to="column", values_to="contents") %>%
  separate_rows(contents, sep=";\t") %>%
  separate_rows(contents, sep="; ") %>%
  separate(contents, into=c("attr", "value"), sep=":", extra="merge") %>%
  mutate(value=str_squish(value)) %>%
  filter(!(is.na(attr) & is.na(value))) %>%
  mutate(value=ifelse(is.na(value), attr, value),
         attr=ifelse(attr==value, column, attr)) %>%
  select(-column) %>%
  dplyr::rename(study_acc=gse, sample_acc=gsm) %>%
  group_by(study_acc, sample_acc, attr) %>% 
  summarise(value=paste(unique(value), collapse=";")) %>%
  ungroup() %>%
  group_split(study_acc)

col_stats <- lapply( 1:length(split_by_gse), function(i) {
  print(i); find_col_stats(split_by_gse[[i]])})
names(col_stats) <- sapply(col_stats, function(i) i$df$study_acc[[1]])


# --- SRA + AE studies --- #
rb_exp_sample_map <- read_csv(sprintf("%s/01_sample_lists/rb_metadata/human_microarray_exp_to_sample.csv", RB.PATH)) %>%
  bind_rows(read_csv(sprintf("%s/01_sample_lists/rb_metadata/human_rnaseq_exp_to_sample.csv", RB.PATH)))

ae_studies <- read_csv("data/array_exp_studies_1231.csv")
sra_studies <- read_csv("data/sra_studies_1231.csv")
ae_sra <- ae_studies %>% select(study_acc) %>% bind_rows(sra_studies %>% select(study_acc))
ae_sra2 <- ae_sra %>% inner_join(rb_exp_sample_map) 
ae_sra2 %>% filter(str_detect(study_acc, "ER|SR|DR")) %>% nrow() # 5989
ae_sra2 %>% filter(!str_detect(study_acc, "ER|SR|DR")) %>% nrow() # 189
length(unique(ae_sra2$sample_acc)) # 3461

ae_sra_sample_dat <- ae_sra2 %>% left_join(all_attrib_clean, by=c("sample_acc"))
missing <- ae_sra_sample_dat %>% filter(is.na(key) & is.na(value))
(missing%>% nrow()) # 2240 are missing
present_sample_dat <- ae_sra_sample_dat %>% filter(!(is.na(key) & is.na(value)))

study_split_aes <- present_sample_dat %>%
  select(-key_clean, -value_clean) %>%
  dplyr::rename(attr=key) %>%
  unique()  %>%
  group_split(study_acc) 

aes_stats <- lapply( 1:length(study_split_aes), function(i) {
  print(i); find_col_stats(study_split_aes[[i]])})
names(aes_stats) <- sapply(aes_stats, function(i) i$df$study_acc[[1]])

# ---- put data together ---- #
all_stats <- c(col_stats, aes_stats)
save(all_stats, file="data/study_sample_info.RData")
cleaned_dat <- do.call(rbind, lapply(all_stats, function(x) x$df_clean %>% 
                                       select(study_acc, sample_acc, "smok", "pkyrs",
                                       "sex", "race_ethnicity", "age", "copd")))

cleaned_dat_comb <- cleaned_dat %>% 
  mutate(across(smok:copd, ~ifelse(. %in% c("-999", "", "NA"), NA, .))) 

summary(cleaned_dat_comb %>%
          mutate(across(c(smok, sex, race_ethnicity, age, copd), as.factor), 
                 across(c(age, pkyrs), as.numeric)))


# do the columns vary
no_vary_comb <- cleaned_dat_comb %>%
  distinct(study_acc, smok) %>% group_by(study_acc) %>% count() %>%
  filter(n==1)

# get all the distinct values for a type
get_distinct_vals <- function(df, col){
  df %>%
    filter(!is.na({{col}})) %>%
    group_by({{col}}) %>%
    count() %>%
    arrange(desc(n)) 
}


cleaned_dat_comb %>% separate_rows(smok, sep=",|;")


extract_from_str <- function(my_str, SEARCH.STR){
  if (!str_detect(my_str, ",|;")){
    return(my_str)
  }
  my_str2 <- str_replace_all(my_str, "[0-9]+", "")
  terms <- unique(sapply(str_split(my_str2, ",|;")[[1]], str_squish))
  my_term <- terms[str_detect(terms, SEARCH.STR)]
  if(length(my_term)==0){
    return(my_str)
  } else {
    return(paste(my_term, collapse=";"))
  }
}

smok_vals <- cleaned_dat_comb %>% get_distinct_vals(smok)
smok_vals2 <- smok_vals %>%
  group_by(smok) %>%
  mutate(smok_clean=extract_from_str(smok, "smok|current|former|never|quit")) %>%
  ungroup() %>% 
  group_by(smok_clean) %>%
  summarise(smok=paste(smok, collapse=";;"),
            n=sum(n)) %>%
  arrange(desc(n))


sex_vals <- cleaned_dat_comb %>% get_distinct_vals(sex)
race_ethnicity_vals <- cleaned_dat_comb %>% get_distinct_vals(race_ethnicity)
age_non_na <- cleaned_dat_comb %>% filter(is.na(as.numeric(age))) %>% get_distinct_vals(age)

smok_vals2 %>% write_csv("data/smok_distinct_vals.csv")
sex_vals %>% write_csv("data/sex_distinct_vals.csv")
race_ethnicity_vals %>% write_csv("data/race_ethnicity_distinct_vals.csv")
age_non_na %>% write_csv("data/age_non_num_distinct_vals.csv")
# --> THEN map

# <-- read in the annotated data

smok_annot <- read_csv("data/smok_distinct_vals_annot1.csv")
sex_annot <- read_csv("data/sex_distinct_vals_annot.csv")
smok_map <- smok_annot %>% select(-smok) %>%
  left_join(smok_vals2 %>% select(smok_clean, smok)) %>%
  separate_rows(smok, sep=";;") %>%
  distinct(annot, smok) 
cleaned_dat_comb2 <- cleaned_dat_comb %>% 
  left_join(smok_map, by=c("smok")) %>% 
  dplyr::rename(smok_orig=smok, smok=annot) %>%
  group_by(study_acc, sample_acc) %>%
  summarise(smok=paste(smok, collapse=";;"),
            smok_orig=paste(smok_orig, collapse=";;")) %>%
  mutate(smok=ifelse(smok=="NA", NA, smok)) %>%
  ungroup()

cleaned_dat_comb2 %>% group_by(study_acc, smok) %>% count() 
diff_smok_stats <- cleaned_dat_comb2 %>% 
  group_by(study_acc,smok) %>%  
  count() %>%
  ungroup() %>%
  group_by(study_acc) %>%
  summarise(smok_grps=paste(unique(smok), collapse=";"),
            num_diff=n())

# --- add this as a column --> will help ID studies that have smoking history vs not! --- #
manual_annot <- read_csv("data/manual_annot_1231.csv")
manual_annot2 <- manual_annot %>% 
  left_join(diff_smok_stats %>% select(study_acc, smok_grps)) 
manual_annot2 %>% write_csv("data/manual_annot_1231_v2.csv")
