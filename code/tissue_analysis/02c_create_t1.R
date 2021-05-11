# --- Create table 1 --- #
library("tidyverse")

# ---- 0. resolve messed up SL? grab missing? ---- #
# (see table1_check_sl.R)
filt_pdat <- str_replace_all(list.files("data/pdata_filt/"), ".csv", "")
my_f <- lapply(filt_pdat, function(x) read_csv(sprintf("data/pdata_filt/%s.csv", x)))
names(my_f) <- filt_pdat
names(my_f)
# deal with multi-tissue
my_f[["gse8987_1"]] <- my_f[["gse8987_1"]] %>% 
  mutate(study_acc=paste(study_acc, "1", sep="_"))
my_f[["gse8987_2"]] <- my_f[["gse8987_2"]] %>% 
  mutate(study_acc=paste(study_acc, "2", sep="_"))
# manual edit: GSE16149: all female
my_f[["gse16149"]] <- my_f[["gse16149"]] %>% 
  mutate(sex_lab="female")

# 1. remove duplicates
# GSE21862, GSE19027, GSE46699
my_f[["gse21862"]]  <- my_f[["gse21862"]] %>% 
  group_by(`subject id`) %>%
  sample_n(1) %>% 
  ungroup()

my_f[["gse19027"]] <- my_f[["gse19027"]] %>% 
  group_by(patient_id) %>%
  sample_n(1) %>% 
  ungroup()

my_f[["gse46699"]] <- my_f[["gse46699"]] %>% 
  group_by(patient)  %>%
  sample_n(1) %>% 
  ungroup()


# create the combined tables with counts

my_df <- do.call(rbind, lapply(my_f, function(x)
  x %>% select(study_acc, sample_acc, sex_lab, smok, tissue)))

count_table <- my_df %>%   
  filter(!is.na(sex_lab),sex_lab != "unlabeled", !is.na(smok), smok != "FS" ) %>%
  group_by(study_acc, tissue, sex_lab, smok) %>% 
  count() %>%
  unite(grp, c(sex_lab, smok)) %>%
  pivot_wider(names_from="grp", values_from="n", values_fill=0) %>%
  mutate(num_S=female_S+male_S, 
         num_NS=female_NS+male_NS,
         num_f=female_S+female_NS,
         num_m=male_S+male_NS,
         tot=num_S+num_NS)

smok_tab <- count_table %>% 
  select(contains("num"), tot)

smok_sex <- count_table %>% 
  filter(across(contains("S"), ~.>=3)) %>%
  select(-contains("num"), -tot) # 9 studies

# grab platform
library('GEOmetadb')
unique(my_df$study_acc)
list_studies <- paste(c(unique(my_df$study_acc), "GSE8987" ), collapse="','")
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
plat <- dbGetQuery(con, sprintf("SELECT gse,gpl FROM gse_gpl WHERE gse IN ('%s');",
                     list_studies))

dbDisconnect(con)

library('readxl')
study_annot <- read_excel("table1_study_info.xlsx") %>%
  select(-remove, -`duplicate?`, -notes)

smok_tab2 <- smok_tab %>% left_join(
  plat %>%
    group_by(gse) %>% 
    summarize(gpl=paste(gpl, collapse=";")) %>%
    mutate(gse=ifelse(gse=="GSE8987", "GSE8987_1;GSE8987_2", gse)) %>%
    separate_rows(gse, sep=";"),
  by=c("study_acc"="gse")
) %>%
  rename(platform=gpl) %>%
  select(study_acc, tissue, platform, everything()) %>% 
  left_join(study_annot %>% select(-tissue), by="study_acc")

smok_tab3 <- smok_tab2 %>%
  filter(num_S >=5, num_NS >= 5) %>%
  ungroup() %>%
  mutate(study_acc=ifelse(study_acc %in% c("GSE8987_1", "GSE8987_2"), 
                          "GSE8987", study_acc)) %>%
  left_join(incl_studies %>% select(study_acc, title)) %>%
  select(study_acc, title, tissue, platform, type, everything()) %>%
  rename(num_samples=tot) %>%
  mutate(covars=str_replace_all(covars, ",", ";")) 
smok_tab3 %>% View()

count_table %>% 
  mutate(frac_f_s=female_S/(female_S+male_S),
         frac_f_ns=female_NS/(female_NS+male_NS)) %>%
  mutate(across(contains("frac"), ~signif(., 2))) %>%
  select(study_acc, contains("frac")) %>%
  write_csv("data/f_frac.csv")
smok_tab3 %>%
  write_csv("data/supp_tables/s4_list_smok_studies.csv")

smok_tab3 %>% group_by(tissue) %>% count() %>% arrange(desc(n))
# 3 - lung, PBMC, 2 - alveolar macrophages, BE, BM, whole blood 

smok_sex %>% left_join(
  plat %>%
    group_by(gse) %>% 
    summarize(gpl=paste(gpl, collapse=";")),
  by=c("study_acc"="gse")
) %>%
  rename(platform=gpl) %>%
  select(study_acc, tissue, platform, everything()) %>%
           write_csv("data/supp_tables/s4b_breakdown_smok_sex.csv")

samples_kept <- my_df %>% 
  mutate(study_acc=ifelse(str_detect(study_acc,"GSE8987"),"GSE8987", 
                          study_acc ) )%>%
  filter(study_acc %in% smok_tab3$study_acc)
samples_kept %>%
  write_csv("data/filt_samples.csv")

# fraction missing initial labels?

# sex breakdown by study
ggplot(smok_sex %>%
         pivot_longer(
    c("female_NS","female_S", "male_NS", "male_S"), 
    names_to="group", 
    values_to="count") %>%
      mutate(group=str_replace_all(group, "_NS", " nonsmoker")) %>%
      mutate(group=str_replace_all(group, "_S", " smoker")), 
  aes(x=study_acc, y=count, fill=group)) +
  ylab("number of samples")+
  xlab("")+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

ggsave("figures/paper_figs/grp_breakdown.png")
  