
library('tidyverse')

load("data/eset_comb.RData") # --> eset_comb, dup_split

condense_cols <- function(x){
  paste(sort(unique(x[!is.na(x)])), collapse=";")
}

comb_metadata <- read_csv("data/comb_sample_epithelium_metadata.csv")
dup_dat <- dup_split %>% 
  mutate(grp=1:n()) %>%
  select(-n) %>%
  rename(geo_accession=sample) %>%
  separate_rows(geo_accession) %>% 
  left_join(comb_metadata, by=c("geo_accession")) %>%
  group_by(grp) %>%
  summarize(across(everything(), ~condense_cols(.))) 

dup_dat %>% filter(str_detect(metadata_sex,";") |
                     str_detect(race_ethnicity,";") |
                     str_detect(smok,";") | # one
                     str_detect(age,";")  | #one
                     str_detect(expr_sex, ";"))

# these are fine
dup_dat %>% filter(str_detect(tissue, ";") ) %>% group_by(tissue) %>% count()

dup_dat2 <- dup_dat %>%
  mutate(
    tissue=case_when(
      ! str_detect(tissue, ";")  ~ tissue,
      tissue=="airway epithelium;small airway epithelium" ~ "small airway epithelium",
      tissue=="trachea;trachea epithelium" ~ "trachea epithelium" 
      ), # take the more specific tissue
    age=case_when(
      ! str_detect(age, ";") ~ age,
      TRUE ~ as.character(mean(as.numeric(str_split(age, ";")[[1]])))
      ) # take the mean of the ages
    ) %>%
  filter(!str_detect(smok, ";")) # remove conflicting

# add the phenotype data
present_cols <- colnames(eset_comb)


dup_dat2_long <- dup_dat2 %>% 
  separate_rows(geo_accession, sep=";") %>% select(-grp) %>%
  mutate(across(c(age, pack_years, pred), as.numeric))
comb_metadata2 <- dup_dat2_long %>%
  bind_rows(comb_metadata %>% anti_join(dup_dat2_long, by=c("geo_accession"))  )

comb_metadata3 <- comb_metadata2 %>% 
  filter(geo_accession %in% present_cols) %>%
  arrange(geo_accession) %>%
  mutate(across(everything(), ~ifelse(.=="", NA, .))) %>%
  mutate(tissue=ifelse(tissue=="trachea", "trachea epithelium", tissue)) %>%
  mutate(across(c(smok, tissue, metadata_sex, race_ethnicity,
                  expr_sex), as.factor))

summary(comb_metadata3 %>% select(-id, -geo_accession, -study))
eset_comb2 <- eset_comb[,comb_metadata3$geo_accession]
# remove missing data
row_miss <- apply(eset_comb2, 1, function(x) sum(is.na(x)))
eset_comb3 <- eset_comb2[row_miss==0,]

save(eset_comb3, comb_metadata3, file="data/ae_eset_pheno.RData")



#

pcs <- prcomp(t(eset_comb3), rank.=4, scale=TRUE)
pcs2 <- data.frame(pcs$x)
pcs2$smoking <- comb_metadata3$smok
pcs2$sex <- comb_metadata3$expr_sex
pcs2$tissue <- comb_metadata3$tissue

ggplot(pcs2, 
         aes(x=PC1, y=PC2, col=sex, shape=smoking))+
    geom_point(alpha=0.5)+
    theme_bw()
ggsave("figures/ae_pcs_disc.pdf")
ggplot(pcs2, 
       aes(x=PC1, y=PC2, col=tissue))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/ae_pcs_disc_tissue.pdf")


# separate out: bronchial epithelium, alveolar macrophages
pcs2 %>% filter(!tissue %in% c("bronchial epithelium",
                               "alveolar macrophages")) %>%
  ggplot(aes(x=PC1, y=PC2, col=tissue, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/ae_only_pcs_disc.pdf")
# // TODO should we remove

# TODO - demographic information checks!
# age
# race-ethnicity
# missingness of data


# ... rerun expression analysis? ... #
head(eset_comb3[,1:5])
ae_only_meta <- comb_metadata3 %>% filter(!tissue %in% c("bronchial epithelium",
                                                                "alveolar macrophages")) 
ae_only <- eset_comb3[, ae_only_meta$geo_accession] # 539 vs 671

library(limma)
ae_only[,1:5]

smok <- factor(ae_only_meta$smok)
sex <- factor(ae_only_meta$expr_sex)

#batch <- factor(blood_pDat2$batch)
design<- model.matrix(~smok+sex+sex*smok)
fit <- lmFit(ae_only, design)
fit <- eBayes(fit)

# add in HGNC symbol
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
convert_genes <- getBM(attributes = c("hgnc_symbol", "entrezgene_id", "chromosome_name"), 
        filters = "entrezgene_id", values = rownames(ae_only),
        mart=mart)

add_gene_info <- function(df){
  df <- data.frame(df)
  df$entrezgene_id <- rownames(df)
  df %>% mutate(entrezgene_id=as.numeric(entrezgene_id)) %>%
    left_join(convert_genes) %>% 
    dplyr::rename(gene=hgnc_symbol) %>%
    dplyr::select(-entrezgene_id) %>%
    as_tibble()
}

# adj.P < 0.05 

res_smok <- topTable(fit, coef="smokS", number=nrow(ae_only)) %>% add_gene_info()
# 7,607 

res_smok2 <- res_smok %>% filter(adj.P.Val < 0.01 & 
                                 chromosome_name %in% c(1:22, "X", "Y")) %>%
  mutate(chromosome_name=factor(chromosome_name, 
                                levels= c(1:22, "X", "Y")))
res_smok2 %>%
  ggplot(aes(x=chromosome_name))+
  geom_histogram(stat="count")+
  theme_bw()+
  ylab("number of genes")+ 
  ggtitle("smoking-related genes (P<0.01)")
ggsave("figures/ae_chr_dist_smok.png")

res_sex <- topTable(fit, coef="sexmale", number=nrow(ae_only)) %>% add_gene_info()
# 767
res_sex2 <- res_sex %>% filter(adj.P.Val < 0.01 & 
                                 chromosome_name %in% c(1:22, "X", "Y")) %>%
                        mutate(chromosome_name=factor(chromosome_name, 
                                                      levels= c(1:22, "X", "Y")))

res_sex2 %>%
  ggplot(aes(x=chromosome_name))+
  geom_histogram(stat="count")+
  theme_bw()+
  ylab("number of genes")+
  ggtitle("sex-related genes (P<0.01)")
ggsave("figures/ae_chr_dist_sex.png")

res_int <- topTable(fit, coef="smokS:sexmale", number=nrow(ae_only)) %>% add_gene_info()

res_smok %>% filter(adj.P.Val < 0.05 ) %>% nrow()
res_int %>% filter(gene %in% c("RALGDS", "DNAI2", "DGKZ"))

res_int2 <- res_int %>% filter(adj.P.Val < 0.01 & 
                                   chromosome_name %in% c(1:22, "X", "Y")) %>%
  mutate(chromosome_name=factor(chromosome_name, 
                                levels= c(1:22, "X", "Y")))
res_int2 %>%
  ggplot(aes(x=chromosome_name))+
  geom_histogram(stat="count")+
  theme_bw()+
  ylab("number of genes")+
  ggtitle("smoking*sex genes (P<0.01)")
ggsave("figures/ae_chr_dist_smok_sex.png")

# volcano plot
library(ggrepel)
res_int2 <- res_int %>% 
  mutate(gene_grp=case_when(
  adj.P.Val > 0.05 ~ "not sig",
  logFC < 0 ~ "down",
  logFC > 0 ~ "up")) %>%
  mutate(gene_grp=factor(gene_grp, levels=c("down", "up", "not sig")))
ggplot(res_int2, 
       aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(aes(col=gene_grp), alpha=0.5)+
  geom_label_repel(data=res_int2 %>% head(20), aes(label=gene), size=3)+
  scale_color_manual(values=c("red", "blue", "gray"))+
  theme_bw()+
  ggtitle("smoking x sex interaction effects")
ggsave("figures/volcano_ae_sex_smoking.pdf")


res_sex2 <- res_sex %>% 
  mutate(gene_grp=case_when(
    adj.P.Val > 0.05 ~ "not sig",
    logFC < 0 ~ "down",
    logFC > 0 ~ "up")) %>%
  mutate(gene_grp=factor(gene_grp, levels=c("down", "up", "not sig")))
ggplot(res_sex2, 
       aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(aes(col=gene_grp), alpha=0.5)+
  geom_label_repel(data=res_sex2 %>% head(20), aes(label=gene), size=3)+
  scale_color_manual(values=c("red", "blue", "gray"))+
  theme_bw()+
  ggtitle("sex effects")
ggsave("figures/volcano_ae_sex.pdf")

res_smok2 <- res_smok %>% 
  mutate(gene_grp=case_when(
    adj.P.Val > 0.05 ~ "not sig",
    logFC < 0 ~ "down",
    logFC > 0 ~ "up")) %>%
  mutate(gene_grp=factor(gene_grp, levels=c("down", "up", "not sig")))
ggplot(res_smok2, 
       aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(aes(col=gene_grp), alpha=0.5)+
  geom_label_repel(data=res_smok2 %>% head(20), aes(label=gene), size=3)+
  scale_color_manual(values=c("red", "blue", "gray"))+
  theme_bw()+
  ggtitle("smoking effects")
ggsave("figures/volcano_ae_smoking.pdf")


# GO analysis - GSVA
library("GSVA")

# how does this match the previous analysis?


# LIONESS COEXPR NETWORKS



