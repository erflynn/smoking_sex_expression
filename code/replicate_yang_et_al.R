
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
# replication GSE7895
library(tidyverse)
library(readxl)
supp_files1 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=1, skip=2, col_names=TRUE)

head(supp_files1)
incl_files <- supp_files1 %>%
  filter(`Microarray Platform`=="U133 Plus 2.0 Array")
list_geo_studies <- incl_files %>% pull(`GEO Accession`)


supp_files4 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=4, skip=4, col_names=TRUE)
colnames(supp_files4)
disc <- supp_files4[,1:6]
rep <- supp_files4[,8:13] 
colnames(disc) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
colnames(rep) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
fill_empty_cells <- function(df){
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
disc1 <- fill_empty_cells(disc)

rep1 <- fill_empty_cells(rep %>% filter(!is.na(pval)))
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

overlapping <- intersect(disc2 %>% distinct(gene) %>% pull(gene),
          rep2 %>% distinct(gene) %>% pull(gene))
# 333 genes

# look at directionality
both_g <- disc2 %>% inner_join(rep2 %>% select(-chromosome), by="gene") 
ggplot(both_g, aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.7)+
  theme_bw()+
  ylab("logFC in replication")+
  xlab("logFC in discovery")
# --> 184 (149 did not)
ggsave("figures/yang_disc_valid_compare.png")
# ... this does not impress me... hmm.

# what is the probability of each direction?
# correlation of betas

# what abt grouping by gene
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

# 32 of the overlapping genes have probes in opposite directions  in disc
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
  summarize(num_neg=sum(logFC<0), num_pos=sum(logFC>0), n=unique(n)) %>%
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
