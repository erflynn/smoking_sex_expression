library(tidyverse)

all_pdat <- lapply(list.files("data/pdata_filt2_v2/"), 
                   function(x) read_csv(sprintf("data/pdata_filt2_v2/%s", x)))
studies <- toupper(str_replace_all(list.files("data/pdata_filt2_v2/"), ".csv", ""))
names(all_pdat) <- studies

lapply(all_pdat, colnames)
new_study_df <- do.call(rbind, lapply(studies , function(study) {
  x <- all_pdat[[study]]
  if (!"metadata_sex" %in% colnames(x)){
    x <- x %>% mutate(metadata_sex="")
  }
  if (!"tissue" %in% colnames(x)){
    x <- x %>% mutate(tissue="")
  }
  if ("gsm" %in% colnames(x)){
    x <- x %>% rename(sample_acc=gsm)
  }
  x %>% 
    mutate(study_acc=study) %>%
    select(study_acc, sample_acc, smok, metadata_sex, sex, tissue)
})
)

df <- data.frame(new_study_df)
df2 <- df %>% 
  mutate(tissue=case_when(
    study_acc=="GSE20681" | tissue=="whole blood" ~ "blood",
    tissue=="bronchial epithelium" ~ "airway epithelium",
    str_detect(tissue, "cortex") ~ "brain",
    TRUE ~ tissue)) %>%
  group_by(tissue, study_acc, smok, sex)  %>%
  count() %>%
  unite(grp, c("smok", "sex")) %>%
  pivot_wider(names_from=grp, values_from=n) %>%
  mutate(NS=NS_female + NS_male,
         S=S_female+S_male,
         total=NS+S) 
df2 %>% write_csv("table1_part2.csv")

# get title, platform
post_annot <- read_csv("data/manual_annot_extra_050221_v2.csv")

df3 <- df2 %>% inner_join(post_annot %>% filter(study_acc %in% studies) %>%
                            select(study_acc, title) %>%
                            left_join(gse_gpl, by=c("study_acc"="gse")) %>%
                            mutate(gpl=ifelse(is.na(gpl), "GPL570", gpl))) %>%
  select(tissue, study_acc, title, gpl, total, S, NS, everything())


df3 %>% write_csv("table_part2_w_gpl.csv")
df3 <- read_csv("table_part2_w_gpl.csv")
annot_part2 <- read_csv("table1_part2_v2.csv")


# load the platform data and save it
# gse_to_gpl <- df3 %>% filter(gpl !="GPL570") %>% ungroup() %>%
#   select(study_acc, gpl) %>%
#   group_by(gpl) %>%
#   sample_n(1)
# lapply(gse_to_gpl$study_acc, function(gse){
#   load(sprintf("data/small_gses2/%s.RData", gse))
#   probe_gene <- tibble("probes"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
#   save(probe_gene, file=sprintf("ref/%s_probe_gene.RData", tolower(my_gse$platform)))
# })


library('readxl')
study_annot <- read_excel("table1_study_info.xlsx") %>%
  select(-remove, -`duplicate?`, -notes)
form_list <- read_csv("data/supp_tables/s4_list_smok_studies_formatted.csv")
head(form_list)
form_list %>% distinct(tissue)
form_list %>%
  filter(tissue %in% c("whole blood", "pbmc")) 
form_list %>%
  filter(tissue %in% c("alveolar macrophages"))

ae_be <- form_list %>%
  filter(tissue %in% c("bronchial epithelium", "airway epithelium")) %>%
  View()

part1_info <- form_list %>% select(study, title, tissue, platform, mismatch_sex_lab, `additional phenoty[es`) %>%
  rename(covars=`additional phenoty[es`, num_mismatch=mismatch_sex_lab)

part1_counts <- form_list %>% select(study, samples, smokers, `non-smokers`,
                                     contains("frac"))
annot2 <- annot_part2 %>% select(tissue, study_acc, num_mismatch, covars)
part2_info <- annot2 %>% left_join(df3 %>% select(study_acc, title, gpl)) %>% rename(study=study_acc, platform=gpl) 
part2_counts <- df3 %>% 
  mutate(frac_female_s=S_female/S,
         frac_female_ns=NS_female/NS) %>%
  rename(study=study_acc, samples=total, smokers=S, `non-smokers`=NS) %>%
  select(colnames(part1_counts))

part1.2 <- part1_info %>% mutate(
  tissue=case_when(
    tissue=="bronchial epithelium" ~ "airway epithelium",
    tissue=="whole blood" ~ "blood - whole",
    tissue=="pbmc" ~ "blood - pbmcs",
    tissue =="b cell" ~ "blood - b cells",
    tissue =="leukocyte" ~ "blood - leukocytes",
    tissue=="trachea" ~ "trachea epithelium",
    tissue=="oral" ~ "oral cavity",
    tissue=="hippocampus" ~ "brain - hippocampus",
    TRUE ~ tissue
  )
)

comb_info <- part1.2 %>%
  bind_rows(part2_info %>% select(colnames(part1.2))) %>%
  mutate(tissue=ifelse(tissue=="brain", "brain - prefrontal cortex", tissue)) %>%
  select(tissue, everything()) %>% 
  arrange(tissue, study)
counts <- part1_counts %>% bind_rows(part2_counts) %>% 
  filter(study!="GSE44456",study!="GSE55962") %>%
  mutate(female_NS=round(`non-smokers`*frac_female_ns), 
         female_S=round(smokers*frac_female_s),
         male_NS=`non-smokers`-female_NS,
         male_S=smokers-female_S) %>%
  dplyr::select(-contains("frac"))
counts %>% write_csv("data/supp_tables/counts_by_study.csv")
median(counts$smokers)
summary(counts)
counts2 <- counts  %>% group_by(study) %>%
  mutate(s_ns=chisq.test(c(smokers, `non-smokers`))$p.value,
         cat=chisq.test(matrix(c(female_NS, female_S, male_NS, male_S), nrow=2, byrow=T))$p.value)
counts2 %>% dplyr::select(study, smokers, `non-smokers`, s_ns, cat) %>%
  filter(s_ns < 0.05)

counts2 %>% dplyr::select(study, female_NS, female_S, male_NS, male_S, s_ns, cat) %>%
  filter(cat < 0.05)

comb_info2 <- comb_info %>% 
  inner_join(part1_counts %>% bind_rows(part2_counts), by="study") %>%
  select(-num_mismatch, -covars, everything()) %>%
  filter(study != "GSE55962") %>%
  mutate(smokers=sprintf("%s (%s)", smokers, round(frac_female_s*smokers), 0),
         `non-smokers`=sprintf("%s (%s)", `non-smokers`, round(frac_female_ns*`non-smokers`), 0)) %>%
  select(-contains("frac")) 

comb_info2 %>% write_csv("data/table1b_updated.csv")

comb_info3 <- comb_info2 %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex", "brain - hippocampus")))
