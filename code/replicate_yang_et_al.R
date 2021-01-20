# code for replicating Yang et al
# Steps:
# 1. Make list of studies and samples
#   - do we get the same n?
# 2. Grab raw cel files --> RMA --> report probe level
# 3. Visualize
# 4. DE analysis w covars they used
# 5. Add date to covariate analysis

# 16 studies
# sex, smoking status, age, COPD status, ethnicity and pack-years were available for 211 subjects
# (never smokers n=68; current smokers n=143) after removing duplicate data and outliers
#
# validation GSE7895
library(tidyverse)
library(readxl)
supp_files1 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=1, skip=2, col_names=TRUE)

head(supp_files1)
incl_files <- supp_files1 %>%
  filter(`Microarray Platform`=="U133 Plus 2.0 Array")
list_geo_studies <- incl_files %>% pull(`GEO Accession`)

# get the lists of files associated with these data
library(GEOmetadb)
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
gse_gsm <- dbGetQuery(con, sprintf("SELECT gse, gsm FROM gse_gsm WHERE gse IN ('%s');", 
                                   paste(list_geo_studies, collapse="','")))

gsm2 <- dbGetQuery(con, sprintf("SELECT gsm, title, source_name_ch1, description, 
characteristics_ch1 FROM gsm 
                        WHERE gsm IN ('%s');", paste(unique(gse_gsm$gsm), collapse="','")))
dbDisconnect(con)

# clean up AE phenotype data
gsm2.1 <- gsm2 %>%
  separate_rows(characteristics_ch1, sep=";\t") %>%
  mutate(characteristics_ch1=tolower(characteristics_ch1)) %>%
  separate(characteristics_ch1, into=c("key", "value"), sep=": ") %>%
  dplyr::select(-title, -source_name_ch1, -description) %>%
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
    `copd status`=="yes" ~ "yes",
    smoking == "copd" ~ "yes",
    smoking == "early-copd" ~ "early",
    TRUE  ~ "no"
  )) %>%
  mutate(smoking=case_when(
    smoking %in% c("non-smoker", "nonsmoker", "ns") ~ "NS",
    smoking %in% c("smoker", "s") ~ "S"
  )) %>%
  dplyr::select(-`copd status`) %>%
  mutate(pack_years=as.numeric(str_replace_all(pack_years, " pack-years", ""))) %>%
  dplyr::select(gsm, age, sex, smoking, race_ethnicity, copd, pack_years, everything())

gsm2.2 %>% fct_summ()
# none of the DGM IDs are replicated

gsm2.3 <- gsm2.2 %>% filter(!is.na(smoking), copd=="no", 
                  !is.na(race_ethnicity), !is.na(age)) %>%
  dplyr::select(gsm, age, sex, smoking, race_ethnicity, pack_years) %>%
  filter(smoking=="NS" | (smoking=="S" & !is.na(pack_years)))
table(gsm2.3$smoking) # 212 S, 153 NS
length(unique(gsm2.3$gsm))

# get download information
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
download_info = dbGetQuery(con, sprintf("SELECT gsm, gpl, submission_date, supplementary_file FROM gsm WHERE gsm IN ('%s')", 
                        paste(gsm2.3$gsm, collapse="','")))
dbDisconnect(con)
download_info2 <- download_info %>% s
eparate_rows(supplementary_file, sep=";\t") %>%
  filter(str_detect(supplementary_file, "CEL"))
gsm2.4 <- gsm2.3 %>% left_join(download_info %>% 
                                 dplyr::select(gsm, submission_date))

download_info2 %>% write_csv("data/list_ae_to_download_replication.csv")




# ------- READ IN significant genes tables ------ #

fill_empty_cells <- function(df){
  # fill empty cells w previous values
  df2 <- df %>% 
    mutate(gene=case_when(
      is.na(gene) ~ lag(gene),
      TRUE ~ gene),
      chromosome=case_when(
        is.na(chromosome) ~ lag(chromosome),
        TRUE ~ chromosome
      )) 
  num_nas <- length(which(is.na(df2$gene)))
  if(num_nas==0){
    return(df2)
  }
  return(fill_empty_cells(df2))
}

format_sig_genes <- function(df){
  disc <- df[,1:6]
  rep <- df[,8:13] 
  colnames(disc) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
  colnames(rep) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
  disc1 <- fill_empty_cells(disc %>% filter(!is.na(pval)))
  rep1 <- fill_empty_cells(rep %>% filter(!is.na(pval)))
  return(list("disc"=disc1, "rep"=rep1))
}
supp_files2 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=2, skip=4, col_names=TRUE)
supp_files3 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=3, skip=4, col_names=TRUE)
supp_files4 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=4, skip=4, col_names=TRUE)

smok_dr <- format_sig_genes(supp_files2)
sex_dr <- format_sig_genes(supp_files3)
int_dr <- format_sig_genes(supp_files4)

disc1 <- int_dr$disc
rep1 <- int_dr$rep

disc_sex1 <- smok_dr$disc
disc_smok1 <- smok_dr$disc
rep_sex1 <- sex_dr$rep
rep_smok1 <- smok_dr$rep


stopifnot(length(which(is.na(disc1$gene)))==0)
stopifnot(length(which(is.na(rep1$gene)))==0)

# fix date genes
which(sapply(rep1$gene, function(x) !is.na(as.numeric(x))))
# "43350" == SEPTIN7   (ch7)
# "43164" ==  MARCH5 (ch10) MARCHF5
# "43160" == Mar1 or March1  (ch1) MARC1 --> MTARC1
rep2 <- rep1 %>%
  mutate(gene=case_when(
    gene=="43350" ~ "SEPTIN7",
    gene=="43164" ~ "MARCHF5",
    gene=="43160" ~ "MTARC1",
    TRUE ~ gene
  )) %>%
  filter(!is.na(pval))
stopifnot(length(which(sapply(rep2$gene, function(x) !is.na(as.numeric(x)))))==0)


which(sapply(disc1$gene, function(x) !is.na(as.numeric(x))))
# "43167"== Mar-8 (ch10) MARCHF8
# "43355" = 12-Sep (ch16) SEPTIN12
# "43168" = Mar 9 (ch12) MARCHF12
# "43349" = 6-sep (X) SEPTIN6
# "43165" = 6-Mar (5) MARCHF6
# "43354" = 11-Sep (4) SEPTIN11
# "43350" = 7-Sep (ch7) SEPTIN7
# "43352" = 9-Sep (17) SEPTIN9
disc2 <- disc1 %>%
  mutate(gene=case_when(
    gene=="43167" ~ "MARCHF8",
    gene=="43355" ~ "SEPTIN12",
    gene=="43168" ~ "MARCHF12",
    gene=="43349" ~ "SEPTIN6",
    gene=="43165" ~ "MARCHF6",
    gene=="43354" ~ "SEPTIN11",
    gene=="43350" ~ "SEPTIN7",
    gene=="43352" ~ "SEPTIN9",
    TRUE ~ gene
  ))
stopifnot(length(which(sapply(disc2$gene, function(x) !is.na(as.numeric(x)))))==0)


disc2 %>% write_csv("ref/yang_int_disc.csv")
rep2 %>% write_csv("ref/yang_int_rep.csv")


# look at overlap btw discovery and validation
overlapping <- intersect(disc2 %>% distinct(gene) %>% pull(gene),
          rep2 %>% distinct(gene) %>% pull(gene))
# 333 genes


# look at directionality
both_g <- disc2 %>% inner_join(rep2 %>% dplyr::select(-chromosome), by="gene") 
ggplot(both_g, aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.7)+
  theme_bw()+
  ylab("logFC in yang replication")+
  xlab("logFC in yang discovery")
# --> 184 (149 did not)
cor.test(both_g$logFC.x, both_g$logFC.y, method="kendall")
ggsave("figures/yang_disc_valid_compare.png")

# simulate random and view
#  TODO: sample with replacement!
both_g %>% filter(logFC.x*logFC.y>0) %>% nrow() # 302
both_g2 <- both_g 
sim_counts <- sapply(1:1000, function(x){
  both_g2$logFC.y <- sample(both_g2$logFC.y, nrow(both_g2))
  both_g2 %>% filter(logFC.x*logFC.y>0) %>% nrow()
})

ggplot(both_g2, aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.7)+
  theme_bw()+
  ylab("random logFC")+
  xlab("logFC in yang discovery")
ggsave("figures/random_logFC.png")

ggplot(tibble("num_probes"=sim_counts), aes(x=num_probes))+
  geom_histogram()+theme_bw()+
  geom_vline(xintercept=302, col="red")+
  xlab("Number of same direction probes")+
  ylab("Number of randomized runs")
ggsave("figures/same_dir_probes.png")

ggplot(both_g2, aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.7)+
  theme_bw()+
  ylab("logFC in yang replication")+
  xlab("shuffled logFC")

# look at sex coefficient disc/rep directionality
both_sex <- disc_sex1 %>% 
  inner_join(rep_sex1 %>% dplyr::select(-chromosome), by="gene") 
both_sex %>%
  ggplot(aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.4)+
  theme_bw()+
  ylab("logFC in yang replication")+
  xlab("logFC in yang discovery")
ggsave("figures/rep_fc_sex.png")

both_sex %>% filter(logFC.x*logFC.y > 0) # 70 / 79
both_sex2 <- both_sex 

sim_counts_s <- sapply(1:1000, function(x){
  both_sex2$logFC.y <- sample(both_sex$logFC.y, nrow(both_sex))
  both_sex2 %>% filter(logFC.x*logFC.y>0) %>% nrow()
})

ggplot(tibble("num_probes"=sim_counts_s), aes(x=num_probes))+
  geom_histogram()+theme_bw()+
  geom_vline(xintercept=70, col="red")+
  xlab("Number of same direction probes")+
  ylab("Number of randomized runs")
ggsave("figures/random_sex.png")

# look at directionality of multiple probes across genes
mult_g <- disc2 %>% semi_join(rep2, by="gene") %>%
  arrange(gene) %>%
  group_by(gene) %>%
  mutate(n=n()) %>%
  filter(n>1)

# neg*neg --> pos
# neg*neg*neg --> neg
mult_g2 <- mult_g %>% 
  group_by(gene) %>% 
  summarize(num_neg=sum(logFC<0), num_pos=sum(logFC>0), n=unique(n)) %>%
  arrange(desc(n))

# 32 of the overlapping genes have probes in opposite directionsin disc
disc_opp <- mult_g2 %>% 
  filter(num_neg!= 0 & num_pos!=0) %>%
  pull(gene)

rep_g2 <- rep2 %>% 
  semi_join(disc2, by="gene") %>%
  arrange(gene) %>%
  group_by(gene) %>%
  mutate(n=n()) %>%
  filter(n>1) %>% 
  group_by(gene) %>% 
  summarize(num_neg=sum(logFC<0), 
            num_pos=sum(logFC>0), 
            n=unique(n)) %>%
  arrange(desc(n)) 

# 28 of the overlapping gnes have probes in opp dir in rep
rep_opp <- rep_g2 %>% 
  filter(num_neg!= 0 & num_pos!=0) %>%
  pull(gene)

# filter and get the list
both_g_filt <- both_g %>% 
  filter(!gene %in% rep_opp,
         !gene %in% disc_opp) %>%
  arrange(gene) %>%
  filter(logFC.x*logFC.y>0) 
length(unique(both_g_filt$gene)) # --> 144

both_g_filt %>% 
  select(gene, chromosome, everything()) %>%
  write_csv("ref/disc_rep_filt.csv")


