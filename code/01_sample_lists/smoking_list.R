# smoking_list.R 
# E Flynn
# 4/18/2019, updated 11/15/2019
#
# Get a list of studies that mention smoking or nicotine. 


# 0) all human studies that contain smoking or nicotine
# 1) do we have tissue or text sex labels? cell line? drug?
# 2) GPLs - which are GPLs we know how to map? gpls?
# 3) aggregate study-level + sample-level files
# 4) download all
# 5) use exprsex to sex label --> make the alluvial diagram
# 6) go thru manually and label more

require('tidyverse')
require('GEOmetadb')
source("code/utils/deduplication_utils.R")

# ---- 0. all human studies ---- #
gse_data <- read_csv("data/sample_lists/gse_all_geo_info.csv")
human_rows <- gse_data %>% 
  filter(str_detect(organism, "Homo sapiens"))
human_str <- human_rows %>% 
  group_by(gse) %>% 
  mutate(str=paste(c(title, overall_design, summary), collapse=" ")) 

human_str <- human_str %>% 
  ungroup(gse) 
smok <- human_str %>%
  filter(grepl( "smok|nicotine|tobacco|cigarette", str, ignore.case=TRUE))  # 220


smok %>% write_csv("data/smok_dat/smoking_data.csv")

# ---- 1. tissue or text sex labels? cell line? drug? ---- #
ale_dat <- read_csv("data/ale_combined_data.csv")
cell_line_dat <- read.table("data/labeled_data/cell_line_mapped_gse.txt", 
                            header=TRUE, stringsAsFactors=FALSE)
drug_dat <- read.table("data/labeled_data/drugbank_mapped_gse.txt", 
                       header=TRUE, stringsAsFactors=FALSE)

length(intersect(smok$gse, unique(ale_dat$gse))) # 161
length(intersect(smok$gse, unique(cell_line_dat$gse))) # 220
length(intersect(smok$gse, unique(drug_dat$gse))) # 63


# ---- 2. GPL data? duplications? ---- #

con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")

gse_list_str <- paste(smok$gse, collapse=c("\',\'"))
smokDatGSM <- dbGetQuery(con, 
sprintf("SELECT gse.gse, gsm, gse.title,
organism, gse_gpl.gpl FROM gse JOIN gse_gpl ON gse.gse=gse_gpl.gse JOIN 
                         gpl ON gse_gpl.gpl=gpl.gpl JOIN gse_gsm ON gse_gsm.gse = gse.gse WHERE gse.gse IN (\'%s\')",
        gse_list_str))

smokDatGSMH <- filter(smokDatGSM, organism=="Homo sapiens")
smokDatGSMH <- smokDatGSMH[!duplicated(smokDatGSMH[,c("gse", "gsm")]),] # 218 studies, 17349 samples

write.csv(smokDatGSMH, file="data/smok_dat/smokDat2.csv")
deduplicate("data/smok_dat/smokDat2.csv", "data/smok_dat/smokDat_dedup2.csv")

dedup_smok <- read.csv("data/smok_dat/smokDat_dedup2.csv", stringsAsFactors = FALSE) # 199 studies, 15455 samples
dedup_smok$X <- NULL

list.gpls <- unique(dedup_smok$gpl) # 53
data.frame(list.gpls) %>% write_tsv("data/smok_dat/list_gpls_smok.txt")
gpls_downloaded <- read_csv("data/smok_dat/gpl_download_meta.csv") # from server information

smok_gpl_meta <- gpls_downloaded %>% filter(gpl %in% list.gpls) # 52
setdiff(list.gpls, smok_gpl_meta$gpl) # 1 isn't present... hmm
#  GPL201  - should map, has GB_ACC column (now it does!)

smok_gpl_meta %>% filter(bad & !downloaded) %>% dplyr::select(gpl, title, description) # 2 possible bad
# GPL16956 - lncRNA, only probe seq
# GPL6955  - miRNA

smok_gpl_meta %>% filter(!downloaded & !bad) %>% dplyr::select(gpl, title) # 4 failed to download that we do not expect
# GPL23159 - this should work
# GPL2829 - this should work
# GPL2986 - this should work
# GPL4133 - this should work

smok_gpl_meta %>% filter(downloaded & (missing | err)) %>% dplyr::select(gpl, title) # 1 downloaded but missing?
# GPL18044 - miRNA

smok_gpl_meta %>% filter(downloaded) %>% nrow() # 46 work

# so should be removing 44 GSMs, 3 GSEs for miR/LNC
remove_mir_lnc <- dedup_smok %>% filter(gpl %in% c("GPL16956", "GPL18044", "GPL6955" ))
dedup_smok2 <- dedup_smok %>% filter(!gpl %in% c("GPL16956", "GPL18044", "GPL6955" ))

# // TODO - figure out how to map these 5 missing
#  this accounts for 588 GSMs and 10 studies
missing_gpl <- dedup_smok %>% filter(gpl %in% c("GPL201", "GPL2986", "GPL2829", "GPL23159",  "GPL4133")) %>% 
  dplyr::select(gse, gsm)
length(unique(missing_gpl$gse))
length(unique(missing_gpl$gsm))

# // TODO need to filter to remove less than 10k probes that map and > 1/3 NAs

dedup_smok2 %>% 
  dplyr::select(gse, gpl) %>% 
  unique() %>%  # 196
  write_csv("data/smok_dat/gse_to_download2.csv")

# ---- 3. aggregate to study level! ---- #

smok_ale <- left_join(dedup_smok2 %>% dplyr::select(-organism), 
                       ale_dat %>% dplyr::select(gse, gsm, text_sex, text_tissue_name), 
                       by=c("gse", "gsm")) # 12418

gse_level_annot <- 
  dedup_smok2 %>% 
  dplyr::select(gse) %>% unique() %>%
  left_join(cell_line_dat) %>% 
  left_join(drug_dat) %>%
  dplyr::select(-ATC, -dbID, -accession) %>%
  mutate(cl=str_replace_all(cl, "-| ","")) %>%
  unique() %>%
  dplyr::rename(drug_name=name) %>% 
  group_by(gse) %>%
  summarize(cl=paste(unique(cl[!is.na(cl)]), collapse=";"), 
            cell_line=any(cell_line), 
            drug_name=paste(unique(drug_name[!is.na(drug_name)]), collapse=";")) %>%
  mutate(cell_line=ifelse(cl!="", TRUE, cell_line)) %>%
  ungroup() 



# now summarize the rest
summary_dat <- smok_ale %>% 
  group_by(gse) %>%
  summarize(gpl=paste(unique(gpl), collapse=";"),
            tissue=paste(unique(text_tissue_name[!is.na(text_tissue_name)]), collapse=";"),
            num_samples=n(),
            num_f=sum(text_sex=="F"),
            num_m=sum(text_sex=="M"))

comb_summary_dat <- inner_join(summary_dat, gse_level_annot )

# add in the title
smok_metadat <- smok %>% dplyr::select(gse, title, pubmed_id, submission_date,overall_design, summary)
all_summary_dat <- left_join(comb_summary_dat, smok_metadat)

summarized <- all_summary_dat %>% 
  dplyr::rename(cell_name=cl) %>%
  mutate(study_type= case_when(
    (is.na(num_f) & is.na(num_m)) ~ "unlabeled",
    (!is.na(num_f) & !is.na(num_m) & num_f > 0 & num_m > 0 ) ~ "mixed",
    (!is.na(num_f) & num_f > 0 ) ~ "female-only",
     (!is.na(num_m) & num_m > 0 )~ "male-only")) %>%
   dplyr::select(gse, title, pubmed_id, submission_date, overall_design, gpl, 
         cell_line, cell_name, tissue, drug_name, num_samples, num_f, num_m,
         study_type, summary)

# grab annotations from before and add to the data 
annot_dat <- read_csv("data/smok_dat/smoking_data_summary_annot_1218.csv")

annot_stud <- annot_dat %>% filter(annot_tiss!="") %>% dplyr::select(-"type")

no_annot <- summarized %>% anti_join(annot_stud, by=c("gse"))
no_annot$annot_tiss <- ""
no_annot$exposure <- ""
no_annot$keep <- ""
no_annot$reason <- ""
no_annot$design <- ""

summarized %>% write_csv("data/smok_dat/smoking_data_summary_0109.csv")

comb <- rbind(annot_stud, no_annot)
comb %>% write_csv("data/smok_dat/smoking_data_summary_0109-annot.csv")


# ----- some descriptive analysis ------ #
head(summarized)
ggplot(summarized, aes(x=1))+geom_bar(aes(fill=study_type))+
  ylab("Number of studies")+theme(axis.title.x = element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_blank())
ggplot(summarized, aes(x=1))+geom_bar(aes(fill=cell_line))+
  ylab("Number of studies")+theme(axis.title.x = element_blank(), 
                                  axis.ticks.x=element_blank(), axis.text.x = element_blank())

tissue_counts <- summarized %>% 
  separate_rows(tissue, sep=";") %>% 
  mutate(tissue = ifelse(tissue %in% c("gravid adult", "tissues", "culture medium", "cell culture"), "", tissue)) %>%
  mutate(tissue=case_when(
    tissue %in% c("bronchial epithelium",	"bronchial epithelial cell") ~ "bronchial epithelium",
    tissue %in% c("epithelium", "epithelial cell line", "epithelial cell") ~ "epithelium",
    str_detect(tissue, "umbilical") ~ "umbilical cord",
    str_detect(tissue, "blood") ~ "blood",
    tissue %in% c("keratinocyte", "skin") ~ "skin",
    str_detect(tissue, "oral squamous cell carcinoma") ~ "oral squamous cell carcinoma",
    tissue == "" ~ "unlabeled",
    TRUE ~ tissue)) %>%
  group_by(tissue) %>%
  count() %>% 
  ungroup(tissue) %>%
  mutate(tissue=ifelse(n > 2, tissue, "(<2 counts)")) %>%
  group_by(tissue) %>%
  mutate(n=ifelse(tissue=="(<2 counts)", n(), n)) %>%
  unique()

tissue_counts$tissue <- factor(tissue_counts$tissue, levels=tissue_counts$tissue[order(tissue_counts$n)])
ggplot(tissue_counts, aes(x=tissue, y=n))+geom_histogram(stat="identity")+
  theme(axis.text.x = element_text(angle=90))+
  ylab("Number of studies")

gpl_counts <- summarized %>% group_by(gpl) %>% count() %>% 
  ungroup(gpl) %>%
  mutate(gpl=ifelse(n > 2, gpl, "(< 2 counts)")) %>%
    group_by(gpl) %>%
  mutate(n=ifelse(gpl=="(< 2 counts)", n(), n)) %>%
  unique() %>%
  arrange(desc(n))

gpl_counts$gpl <- factor(gpl_counts$gpl, levels=gpl_counts$gpl[order(gpl_counts$n)])
ggplot(gpl_counts, aes(x=gpl, y=n))+geom_histogram(stat="identity")+
  theme(axis.text.x = element_text(angle=90))+
  ylab("Number of studies")

# oral squamous cell carcinoma
# peripheral blood cell, blood
# umbilical cord, umbilical cord blood
# 
# cancer vs non-cancer
# 
# mouth
#   - oral squamous cell carcinoma
# blood
# lung
#  - bronchial alveolar lavage
#  - respiratory epithelium

# ----- some labels ----- #

all_drug_labels <- read_csv("data/tmp/all_labels_drug.csv")
overlap_gses <- intersect(unique(all_drug_labels$gse), summarized$gse)
smok_lab <- all_drug_labels %>% filter(gse %in% overlap_gses)
lab2 <- smok_lab %>% mutate(text_sex = 
                                     case_when(text_sex=="M" ~ "male",
                                               text_sex=="F" ~"female")) %>%
  filter(!is.na(gse))

text_labels <- lab2 %>% 
  group_by(gse) %>%
  mutate(n_m=sum(text_sex=="male"),n_f=sum(text_sex=="female")) %>%
  select(gse, n_m, n_f) %>% 
  unique() %>% 
  mutate(study_type=case_when(n_m==0 & n_f==0  ~ "unlabeled",
                              n_m==0 & n_f> 0  ~ "female",
                              n_m>0 & n_f== 0  ~ "male",
                              n_m>0 & n_f> 0  ~ "mixed-sex",
                              TRUE ~ "unlabeled"
  )) %>% 
  rename(text_sex=study_type) %>% 
  select(-n_m, -n_f)

table(text_labels$text_sex)


exprsex_labels <- lab2 %>% 
  group_by(gse) %>% 
  mutate(n_m=sum(rank_sex=="male"),n_f=sum(rank_sex=="female")) %>%
  select(gse, n_m, n_f) %>% 
  unique() %>% 
  mutate(study_type=case_when(n_m==0 & n_f==0  ~ "unlabeled",
                              n_m==0 & n_f> 0  ~ "female",
                              n_m>0 & n_f== 0  ~ "male",
                              n_m>0 & n_f> 0  ~ "mixed-sex",
                              TRUE ~ "unlabeled")) %>% 
  rename(expr_sex=study_type)  %>% 
  select(-n_m, -n_f)

table(exprsex_labels$expr_sex)



comb_labels <- inner_join(text_labels, exprsex_labels)
flow_freq_counts <- comb_labels %>% ungroup(gse) %>% group_by(text_sex, expr_sex) %>% mutate(Freq=n()) %>% select(-gse) %>% unique()
comb_labels_long <- comb_labels %>% rename(metadata=text_sex, exprsex=expr_sex) %>% gather(key="labeling_method", value="sex", -gse) %>% mutate(freq=1) %>% ungroup(gse) %>%  mutate(sex=as.factor(sex), labeling_method=factor(labeling_method, levels=c("metadata", "exprsex")), gse=as.factor(gse)) %>% rename(study=gse)

ggplot(comb_labels_long,
       aes(x = labeling_method, stratum = sex, alluvium = study,
           y = freq,
           fill = sex, label = sex)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  xlab("Label source")+ylab("Number of studies")+
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 

sample_lab <- lab2 %>% select(gsm, gse, text_sex, rank_sex) %>%
  unique() %>%
  rename(text=text_sex, exprsex=rank_sex, sample=gsm, study=gse) %>%
  mutate(exprsex=ifelse(is.na(exprsex), "unlabeled", exprsex),
         text=ifelse(is.na(text), "unlabeled", text)
         ) %>%
  pivot_longer(c("text", "exprsex"), 
               names_to="labeling_method", values_to ="sex") %>%
  mutate(freq=1)

sample_lab$labeling_method <- 
  factor(sample_lab$labeling_method, levels=c("text", "exprsex"))

ggplot(sample_lab,
       aes(x = labeling_method, 
           stratum = sex, alluvium = sample,
           y = freq,
           fill = sex, label = sex)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  xlab("Label source")+ylab("Number of samples")+
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") #+ facet_wrap(. ~study)





# how do we want to collapse the tissue labels?

