
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

# get the lists of files associated with these data
library(GEOmetadb)
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
gse_gsm <- dbGetQuery(con, sprintf("SELECT gse, gsm FROM gse_gsm WHERE gse IN ('%s');", 
                                   paste(list_geo_studies, collapse="','")))



load("data/ae_only_eset.RData")
missing_gsm <- gse_gsm %>% anti_join(ae_only_meta, by=c("gsm"="geo_accession"))

gsm2 <- dbGetQuery(con, sprintf("SELECT gsm, title, source_name_ch1, description, 
characteristics_ch1 FROM gsm 
                        WHERE gsm IN ('%s');", paste(unique(gse_gsm$gsm), collapse="','")))

gsm2.1 <- gsm2 %>%
  separate_rows(characteristics_ch1, sep=";\t") %>%
  mutate(characteristics_ch1=tolower(characteristics_ch1)) %>%
  separate(characteristics_ch1, into=c("key", "value"), sep=": ") %>%
  dplyr::select(-title, -source_name_ch1, -description) %>%
  pivot_wider(names_from=key, values_from=value)


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
  separate(`smoking status`, into=c("smoking", "pack_years"), sep=", ", extra="merge") %>%
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
download_info = dbGetQuery(con, sprintf("SELECT gsm, gpl, submission_date, supplementary_file FROM gsm WHERE gsm IN ('%s')", 
                        paste(gsm2.3$gsm, collapse="','")))
download_info2 <- download_info %>% separate_rows(supplementary_file, sep=";\t") %>%
  filter(str_detect(supplementary_file, "CEL"))
gsm2.4 <- gsm2.3 %>% left_join(download_info %>% dplyr::select(gsm, submission_date))

download_info2 %>% write_csv("data/list_ae_to_download_replication.csv")
# ae_only_meta2 <- ae_only_meta %>% 
#   dplyr::select(geo_accession, study) %>%
#   separate_rows(study, sep=";")
# ae_only_meta3 <- ae_only_meta2 %>% filter(study %in% list_geo_studies)
# # ... hmm only 11 studies? #


# ---> download server side
# while IFS=',' read -r gse gpl submission_date fpath
# do
# wget "$fpath";
# done < "../list_ae_to_download_replication.csv"
# <--- read in
library(affy)

setwd("data/cel_files_rep/")
data1 <- ReadAffy(filenames=list.celfiles()) 
save(data1, file="../data_rep.RData")



eset1 <- affy::rma(data1)
save(eset1, file="../eset_rep.RData")

# remove duplicated? or lots of NAs
load("data/replication/eset_rep.RData")
expDat <- exprs(eset1)
dups <- duplicated(t(expDat))
expDat2 <- expDat[,!dups] # removes 21
pDat <- pData(eset1)
pDat2 <- data.frame(pDat[!dups,])
rownames(pDat2) <- rownames(pDat)[!dups]
colnames(pDat2) <- "sample"
pDat2$gsm <- str_replace_all(rownames(pDat2), ".CEL.gz", "")
colnames(expDat2) <- str_replace_all(colnames(expDat2), ".CEL.gz", "")
pDat3 <- pDat2 %>% left_join(gsm2.4)
rownames(pDat3) <- pDat3$gsm

# correlation matrix
pcor2 <- coop::pcor(expDat2)
list_crs <- pcor2[upper.tri(pcor2)]

ggplot(tibble("cor"=list_crs), aes(x=cor))+
  geom_histogram(bins=60)+
  theme_bw()+
  xlab("pairwise correlation")+
  ylab("")
ggsave("figures/ae_pairwise_correlation.png")

which(pcor2<1 & pcor2 > 0.99, arr.ind=T)
# should we combine these?
pDat3 %>% filter(gsm %in% c("GSM101106", "GSM101105"))
pDat3 %>% filter(gsm %in% c("GSM114089", "GSM101115"))
pDat3 %>% filter(gsm %in% c("GSM114090", "GSM101114"))
pDat3 %>% filter(gsm %in% c("GSM219258", "GSM219252")) 
pDat3 %>% filter(gsm %in% c("GSM219263", "GSM219260")) 
pDat3 %>% filter(gsm %in% c("GSM469996", "GSM469994"))

pDat4 <- pDat3 %>% mutate(pack_years=ifelse(is.na(pack_years), 0, pack_years)) %>%
  mutate(across(c(age, pack_years), as.numeric)) %>%
  mutate(across(c(sex, smoking, race_ethnicity), as.factor))
rownames(pDat4) <- pDat4$gsm
expDat3 <- expDat2[,rownames(pDat4)]
# ... skip for now
#dups2 <- duplicated(pDat3[,c("age", "sex", "smoking", "race_ethnicity", "pack_years", "submission_date")])

#pDat3 %>% 
#  group_by(age, sex, smoking,race_ethnicity, pack_years, submission_date) %>% 
#  count() %>%
#  filter(n>1)


# --- check SL  --- #


# --- covariate viz --- #
library(bigpca)
start_time = Sys.time(); 
pcs <- big.PCA(as.matrix(expDat3), pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time
pcs2 <- data.frame(pcs$PCs)
pcs2$gsm <- rownames(pcs2) 
pcs2.2 <- pcs2 %>% left_join(pDat4, by=c("gsm"))

to_exclude <- pcs2.2 %>% filter(PC1 > 0.1)

# no dups - ok with this
pcs2.2 %>% mutate(across(c(PC1, PC2, PC3), ~round(.,digits=4))) %>%
  group_by(PC1, PC2, PC3) %>%
  count() %>%
  filter(n>1)


plotPC3(pcs2.2, sex)
plotPC3(pcs2.2, smoking)
plotPC3(pcs2.2, race_ethnicity)
plotPC3(pcs2.2, submission_date)





pDat5 <- pDat4 %>% anti_join(to_exclude, by="gsm")

pDat5 %>% 
  group_by( smoking, sex, race_ethnicity) %>%
  count() %>%
  pivot_wider(names_from=race_ethnicity, values_from=n)



chisq.test(t(matrix(c(17, 37, 31, 71, 12, 29, 14, 25, 8, 10, 7, 23), nrow=3, byrow=T)))

pDat5 %>%
  mutate(sex=ifelse(sex=="f", "female", "male")) %>%
  unite(group, c(sex, smoking), sep="-") %>%
  ggplot(aes(x=group, y=age))+
           geom_violin()+
  geom_boxplot(width=0.1)+
  theme_bw()+
  ylab("age (years)")+
  xlab("")

summary(aov(age ~ sex + smoking + sex*smoking, data=pDat5))
t.test(pDat5 %>% filter(sex=="f", smoking=="S") %>% pull(age),
       pDat5 %>% filter(sex=="m", smoking=="S") %>% pull(age))
t.test(pDat5 %>% filter(sex=="f", smoking=="NS") %>% pull(age),
       pDat5 %>% filter(sex=="m", smoking=="NS") %>% pull(age))

pDat5 %>% 
  group_by( smoking, sex) %>%
  summarize(n=n(),
            m_age=mean(age),
            sd_age=sd(age, na.rm = T),
            m_pack_years=mean(pack_years),
            sd_pack_years=sd(pack_years, na.rm=T))
rownames(pDat5) <- pDat5$gsm
t.test(pDat5 %>% filter(sex=="f", smoking=="S") %>% pull(pack_years),
       pDat5 %>% filter(sex=="m", smoking=="S") %>% pull(pack_years))

start_time = Sys.time(); 
pcs_v2 <- big.PCA(as.matrix(expDat3[,pDat5$gsm]), pcs.to.keep=3)

end_time=Sys.time()
end_time-start_time
pcs_v2_alt <- prcomp(t(as.matrix(expDat3[,pDat5$gsm])))
                     
pcs22 <- data.frame(pcs_v2$PCs)
pcs22$gsm <- rownames(pcs22) 
pcs22.2 <- pcs22 %>% left_join(pDat5, by=c("gsm"))

plotPC3(pcs22.2, sex)
plotPC3(pcs22.2, smoking)
plotPC3(pcs22.2, race_ethnicity)
plotPC3(pcs22.2, submission_date)

summary(aov(PC1 ~ sex, data=pcs22.2))
summary(aov(PC2 ~ sex, data=pcs22.2))
summary(aov(PC3 ~ sex, data=pcs22.2))

summary(aov(PC1 ~ smoking, data=pcs22.2))
summary(aov(PC2 ~ smoking, data=pcs22.2))
summary(aov(PC3 ~ smoking, data=pcs22.2))

summary(aov(PC1 ~ race_ethnicity, data=pcs22.2))
summary(aov(PC2 ~ race_ethnicity, data=pcs22.2))
summary(aov(PC3 ~ race_ethnicity, data=pcs22.2))

summary(aov(PC1 ~ submission_date, data=pcs22.2))
summary(aov(PC2 ~ submission_date, data=pcs22.2))
summary(aov(PC3 ~ submission_date, data=pcs22.2))

cor.test(pcs22.2$PC1, pcs22.2$age)
cor.test(pcs22.2$PC2, pcs22.2$age)
cor.test(pcs22.2$PC3, pcs22.2$age)

pcs_s <- pcs22.2 %>% filter(smoking=="S")
cor.test(pcs_s$PC1, pcs_s$pack_years)
cor.test(pcs_s$PC2, pcs_s$pack_years)
cor.test(pcs_s$PC3, pcs_s$pack_years)

# --- run models --- #
design_p <- model.matrix(~sex + smoking + sex:smoking + age + pack_years + race_ethnicity,data=pDat5)
head(expDat2[,1:5])
library('limma')
fit_p <- lmFit(expDat3[,pDat5$gsm], design_p)
fit_p <- eBayes(fit_p)
colnames(fit_p$coefficients)
topTable(fit_p, coef="sexm:smokingS") %>% add_gene() %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))
tt_smok <- topTable(fit_p, coef="smokingS", n=nrow(expDat2)) %>%  add_gene() %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))
tt_sex <- topTable(fit_p, coef="sexm", n=nrow(expDat2)) %>%  add_gene() %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))
tt_int <- topTable(fit_p, coef="sexm:smokingS", n=nrow(expDat2)) %>%  add_gene() %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))
library(ggrepel)
volcano_plot_de(tt_int %>% dplyr::select(-gene) %>% 
                  mutate(gene=hgnc_symbol) %>%
                  dplyr::select(-adj.P.Val) %>%
                  dplyr::rename(adj.P.Val=P.Value), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_int2.png")

volcano_plot_de(tt_sex %>% dplyr::select(-gene) %>% 
                  mutate(gene=hgnc_symbol) %>%
                  dplyr::select(-adj.P.Val) %>%
                  dplyr::rename(adj.P.Val=P.Value), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_sex2.png")

volcano_plot_de(tt_smok %>% dplyr::select(-gene) %>% 
                  mutate(gene=hgnc_symbol) %>%
                  dplyr::select(-adj.P.Val) %>%
                  dplyr::rename(adj.P.Val=P.Value), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_smok2.png")

disc2 <- read_csv("ref/yang_int_disc.csv")
tt_int <- data.frame(tt_int)
tt_int$probe <- rownames(tt_int)

tt_w_other <- tt_int %>% left_join(disc2, by="probe") %>% head(50)
head(tt_w_other)
tt_int %>% inner_join(disc2, by="probe") %>% ggplot(aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.5)





# 
design_date <- model.matrix(~sex + smoking + sex:smoking + age + pack_years + race_ethnicity+
                              submission_date,data=pDat5)
fit_d <- lmFit(expDat3[,pDat5$gsm], design_date)
fit_d <- eBayes(fit_d)
colnames(fit_d$coefficients)

tt_int_d <- topTable(fit_d, coef="sexm:smokingS", n=nrow(expDat2))  %>%  add_gene() %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))
tt_sex_d <- topTable(fit_d, coef="sexm", n=nrow(expDat2))  %>%  add_gene() %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))
tt_smok_d <- topTable(fit_d, coef="smokingS", n=nrow(expDat2))  %>%  
  add_gene() %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))

volcano_plot_de(tt_int_d %>% dplyr::select(-gene) %>% 
                  mutate(gene=hgnc_symbol) %>%
                  dplyr::select(-adj.P.Val) %>%
                  dplyr::rename(adj.P.Val=P.Value), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_d_int2.png")

volcano_plot_de(tt_sex_d %>% dplyr::select(-gene) %>% 
                  mutate(gene=hgnc_symbol) %>%
                  dplyr::select(-adj.P.Val) %>%
                  dplyr::rename(adj.P.Val=P.Value), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_d_sex2.png")

volcano_plot_de(tt_smok_d %>% dplyr::select(-gene) %>% 
                  mutate(gene=hgnc_symbol) %>%
                  dplyr::select(-adj.P.Val) %>%
                  dplyr::rename(adj.P.Val=P.Value), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_d_smok2.png")


# --- meta-analyze probes to genes
tt_w_other_d<- tt_int_d %>% left_join(disc2, by="probe") %>% head(50)
tt_int_d %>% inner_join(disc2, by=c("probes"="probe")) %>% 
  ggplot(aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.5)

# comparison plots
int_both1 <- tt_int %>% inner_join(disc2, by=c("probes"="probe")) 
int_both1 %>% 
  ggplot(aes(x=-log10(P.Value), y=-log10(pval)))+
  geom_point(alpha=0.5)+
  ylim(c(0, 7))+
  ylab("-log10(P) for Yang et al discovery")+
  xlab("-log10(P) for repeated analysis")+
  geom_vline(xintercept=1.94)+
  geom_hline(yintercept=1.94)+
  theme_bw()
tt_int2 <- tt_int %>% filter(P.Value < 0.01) # 1531

ggsave("figures/pval_comparison.png")

int_both <- tt_int_d %>% inner_join(disc2, by=c("probes"="probe")) 
int_both %>% 
  ggplot(aes(x=-log10(P.Value), y=-log10(pval)))+
  geom_point(alpha=0.5)+
  ylim(c(0, 7))+
  ylab("-log10(P) for Yang et al discovery")+
  xlab("-log10(P) for repeated analysis with date")+
  geom_vline(xintercept=1.94)+
  geom_hline(yintercept=1.94)+
  theme_bw()
tt_int_d2 <- tt_int_d %>% filter(P.Value < 0.01) # 538
ggsave("figures/pval_comparison_date.png")
length(intersect(tt_int_d2$probes, tt_int2$probes)) # 124

compare_d <- tt_int %>% full_join(tt_int_d, by="probes")
compare_d %>% 
  filter(P.Value.x < 0.01 | P.Value.y < 0.01) %>%
  ggplot(aes(x=-log10(P.Value.x), y=-log10(P.Value.y)))+
  geom_point(alpha=0.5)+
  ylim(c(0, 6))+
  xlab("-log10(P) for initial analysis")+
  ylab("-log10(P) for analysis with date")+
  geom_vline(xintercept=1.94)+
  geom_hline(yintercept=1.94)+
  theme_bw()
ggsave("figures/pval_comparison_ours.png")
compare_d %>% 
  filter(P.Value.x < 0.01 | P.Value.y < 0.01) %>%
  ggplot(aes(x=logFC.x, y=logFC.y))+
  geom_label_repel(data=compare_d %>% 
              filter(P.Value.x < 0.01 | P.Value.y < 0.01) %>%
              filter(logFC.x > abs(0.8) | logFC.y > abs(0.8)), 
              aes(label=hgnc_symbol.x))+
  geom_point(alpha=0.5)+
  xlab("logFC for initial analysis")+
  ylab("logFC for analysis with date")+
  theme_bw()
ggsave("figures/fc_comparison_ours.png")
compare_d %>% 
  filter(P.Value.x < 0.01 | P.Value.y < 0.01) %>%
  mutate(logFC.d=logFC.x-logFC.y,
         logP.d=-log10(P.Value.x)--log10(P.Value.y)) %>%
  ggplot(aes(x=logFC.d, y=logP.d))+
  geom_point(alpha=0.5)+
  xlab("logFC (orig-date)")+
  ylab("-log10 P (orig-date)")+
  theme_bw()
ggsave("figures/fc_p_comparison_ours.png")

f_comp <- tibble("m1"=fit_p$F, "m2"=fit_d$F)
ggplot(f_comp, aes(x=m1, y=m2))+geom_point(alpha=0.5)+ylim(c(0, 500000))+xlim(c(0, 500000))


compare_d %>% arrange(P.Value.x) %>% head(10)

smok_ci_int_p <- topTable(fit_p, coef="sexm:smokingS", number=nrow(expDat3), 
                        confint = TRUE)
smok_ci_int_p <- data.frame(smok_ci_int_p)
smok_ci_int_p$probes <- rownames(smok_ci_int_p)

smok_ci_int_d <- topTable(fit_d, coef="sexm:smokingS", number=nrow(expDat3), 
                          confint = TRUE)
smok_ci_int_d <- data.frame(smok_ci_int_d)
smok_ci_int_d$probes <- rownames(smok_ci_int_d)
both_ci_int <- smok_ci_int_p %>% inner_join(smok_ci_int_d, by="probes") %>%  
  dplyr::select(-probe) %>%
  left_join(probe_gene, by="probes") %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))
both_ci_int %>% arrange(P.Value.x) %>% head(10) %>%
  ggplot(aes(x=logFC.x, y=logFC.y))+
  geom_abline(slope=1, intercept=0, col="gray")+
  geom_point()+
  geom_errorbar(aes(ymin=CI.L.y, ymax=CI.R.y))+
  geom_errorbarh(aes(xmin=CI.L.x, xmax=CI.R.x))+
  ylim(-0.7, 1.2)+
  theme_bw()+
  geom_label_repel(aes(label=hgnc_symbol), alpha=0.5)+
  ylab("logFC with date")+
  xlab("logFC top genes")
ggsave("figures/topFC_genes_ci.png")

both_ci_int %>% arrange(P.Value.y) %>% head(10) %>%
  ggplot(aes(x=logFC.x, y=logFC.y))+
  geom_abline(slope=1, intercept=0, col="gray")+
  geom_point()+
  geom_errorbar(aes(ymin=CI.L.y, ymax=CI.R.y))+
  geom_errorbarh(aes(xmin=CI.L.x, xmax=CI.R.x))+
  ylim(-0.7, 0.7)+
  xlim(-0.7, 0.7)+
  theme_bw()+
  geom_label_repel(aes(label=hgnc_symbol), alpha=0.5)+
  ylab("logFC top genes with date")+
  xlab("logFC without date")
ggsave("figures/topFC_genes_ci_date.png")
smok_ci_int_p %>% mutate(model="baseline") %>% 
  bind_rows(smok_ci_int_d %>% mutate(model="date")) %>%  
  filter(probes %in% head(smok_ci_int_d$probes, 10) |
           probes %in% head(smok_ci_int_p$probes, 10)) %>%
  dplyr::select(-probe) %>%
  left_join(probe_gene, by="probes") %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)),
            by=c("gene"="entrezgene_id")) %>%
  filter(!is.na(hgnc_symbol)) %>%
  ggplot(  aes(x=hgnc_symbol, y=logFC))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=CI.L, ymax=CI.R))+
  coord_flip()+
  theme_bw()+
  facet_grid(.~factor(model))+
  ylab("")
ggsave("figures/comparison_ci_top_genes.png")
  
ggplot(smok_ci_int %>% head(40) %>%
         mutate(gene=factor(gene, levels=rev(c(smok_ci_int %>% head(40))$gene))),
       aes(x=gene, y=logFC, fill=chromosome))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=CI.L, ymax=CI.R))+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(values=c("gray", "blue", "yellow"))

# 0. COMPARE CI for top 10 probes?
# 1. plot example gene
expDat4 <- data.frame(expDat3)
expDat4$probes <- rownames(expDat3)
exp_long <- expDat4 %>% 
  pivot_longer(-probes,
               names_to="gsm", values_to="expr") %>%
  left_join(probe_gene, by="probes") %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)),
            by=c("gene"="entrezgene_id")) %>%
  left_join(pDat5, by=c("gsm"))

# TODO - figure out how to add lines!
exp_long %>%
  filter(!is.na(hgnc_symbol), !is.na(sex), !is.na(smoking)) %>%
  filter(probes %in% head(tt_int$probes, 5)) %>%
  unite(grp, c("smoking", "sex"), remove=FALSE) %>%
  ggplot(aes(y=expr, x=smoking, fill=sex, group=grp))+
  geom_boxplot()+
  facet_wrap(~hgnc_symbol, scales="free")+
  theme_bw()+
  ylab("")+
  xlab("")
ggsave("figures/example_genes.png")



#

# model comparison
probe_list <- exp_long %>% group_split(probes)
run_mod_comp <- function(df){
  mod1 <- lm(expr ~ age + sex+smoking+smoking*sex+race_ethnicity+pack_years, 
             data=df)
  mod2 <- lm(expr ~ age + sex+smoking+smoking*sex+race_ethnicity+pack_years+submission_date, 
             data=df)
  res <- anova(mod1, mod2)
  return(list("probe"=unique(df$probes), "F"=res$`F`[[2]], "p"=res$`Pr(>F)`[[2]], 
              "SS"=res$`Sum of Sq`[[2]]))
}

names(probe_list) <- unique(exp_long$probes)
mod_comp2 <- lapply(probe_list[tt_int2$probes], run_mod_comp)
my_mods <- do.call(rbind,mod_comp2)
df <- data.frame(my_mods)
df2 <- df %>% mutate(probe=unlist(probe),
              F=unlist(F), 
              p=unlist(p),
              SS=unlist(SS)) %>%
  as_tibble()
df2 %>% filter(SS < 0)
df2 %>% distinct() %>% arrange(p) %>% filter(p<0.05/1423)
df2 %>% filter(p<0.05/1531)

# STOP

overlapping <- tt_int_d %>% inner_join(disc2, by=c("probes"="probe")) %>%  
  filter(-log10(P.Value) > 2)

compare_res <- disc2 %>% left_join(tt_int_d, by="probe")
cor.test(log10(compare_res$pval), log10(compare_res$P.Value))
cor.test(compare_res$logFC.x,compare_res$logFC.y)

compare_res_date <- tt_int_d %>% left_join(tt_int, by="probe")
cor.test(log10(compare_res_date$P.Value.x), log10(compare_res_date$P.Value.y))
cor.test(compare_res_date$logFC.x,compare_res_date$logFC.y)

# present_data <- ae_only_meta %>% 
#   semi_join(ae_only_meta3 %>% distinct(geo_accession) ) %>% 
#   filter(!is.na(race_ethnicity)) %>%
#   filter(!str_detect(tissue, "trachea")) # 203
# 
# present_data %>% 
#   select(geo_accession)
# # never smokers n=68; current smokers n=143


# filter for present pheno data

# deduplicate?


# ------- NEXT PART ------ #
supp_files2 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=2, skip=4, col_names=TRUE)

supp_files3 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=3, skip=4, col_names=TRUE)

supp_files4 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=4, skip=4, col_names=TRUE)
colnames(supp_files4)
disc_sex <- supp_files3[,1:6]
rep_sex <- supp_files3[,8:13] 
colnames(disc_sex) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
colnames(rep_sex) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
disc_smok<- supp_files2[,1:6]
rep_smok <- supp_files2[,8:13] 
colnames(disc_smok) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
colnames(rep_smok) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")

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

disc_sex1 <- fill_empty_cells(disc_sex)
disc_smok1 <- fill_empty_cells(disc_smok)
rep_sex1 <- fill_empty_cells(rep_sex)
rep_smok1 <- fill_empty_cells(rep_smok)


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
both_sex <- disc_sex1 %>% inner_join(rep_sex1 %>% dplyr::select(-chromosome), by="gene") 
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
# ... this does not impress me... hmm.
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

# PROBES TO GENES
library(dmetar)
library(meta)
library(metafor)

smok_int_p <- smok_ci_int_p %>% dplyr::select(-probe) %>% 
       left_join(probe_gene, by="probes") %>% 
       left_join(convert_genes %>% 
                   mutate(entrezgene_id=as.character(entrezgene_id)),
                 by=c("gene"="entrezgene_id")) 
# effect size
# SE for effect size - we have 95%  CI = 1.96*SD
smok_int_p2 <- smok_int_p %>%
  mutate(SD=(logFC-CI.L)/1.96) %>%
  dplyr::select(probes, hgnc_symbol, logFC, SD, P.Value) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol!="") %>%
  distinct() %>% 
  group_by(hgnc_symbol) %>%
  mutate(n=n()) 
multi_probe_gene <- smok_int_p2 %>% filter(n>1) %>% group_split(hgnc_symbol)
ma_probes_genes <- function(df){
  ma <- metagen(df %>% pull(logFC),
                df %>% pull(SD),
                studlab=df %>% pull(probes),
                comb.fixed = TRUE,
                comb.random = FALSE,
                method.tau = "DL",
                hakn = FALSE,
                prediction = FALSE,
                sm = "SMD")
  return(list("gene"=unique(df$hgnc_symbol),
              "logFC"=ma$TE.fixed,
              "logFC.l"=ma$lower.fixed,
              "logFC.u"=ma$upper.fixed,
              "p"=ma$pval.fixed))
  
}
ma_vals <- lapply(multi_probe_gene, ma_probes_genes)
mult_gene_df <- data.frame(apply(do.call(rbind, ma_vals) , c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), ~as.numeric(as.character(.)) ))
comb_ma_dat <- smok_int_p2 %>% filter(n==1) %>% 
  dplyr::rename(gene=hgnc_symbol, p=P.Value) %>%
  mutate(logFC.l=logFC-1.96*SD,
         logFC.u=logFC+1.96*SD) %>%
  dplyr::select(colnames(mult_gene_df)) %>%
  mutate(src="single") %>%
  bind_rows(mult_gene_df %>% mutate(src="mult")) %>%
  arrange(p)
comb_ma_dat$adj.p <- p.adjust(comb_ma_dat$p, method="fdr")
comb_ma_dat %>% filter(adj.p < 0.05) %>% nrow() # 179
comb_ma_dat %>% filter(adj.p < 0.05) %>% dplyr::select(-src) %>% head(10)

# TRY W D

smok_int_d <- smok_ci_int_d %>% 
  left_join(probe_gene, by="probes") %>% 
  left_join(convert_genes %>% 
              mutate(entrezgene_id=as.character(entrezgene_id)),
            by=c("gene"="entrezgene_id")) 

smok_int_d2 <- smok_int_d %>%
  mutate(SD=(logFC-CI.L)/1.96) %>%
  dplyr::select(probes, hgnc_symbol, logFC, SD, P.Value) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol!="") %>%
  distinct() %>% 
  group_by(hgnc_symbol) %>%
  mutate(n=n()) 
multi_probe_gene_d <- smok_int_d2 %>% filter(n>1) %>% group_split(hgnc_symbol)
ma_probes_genes <- function(df){
  ma <- metagen(df %>% pull(logFC),
                df %>% pull(SD),
                studlab=df %>% pull(probes),
                comb.fixed = TRUE,
                comb.random = FALSE,
                method.tau = "DL",
                hakn = FALSE,
                prediction = FALSE,
                sm = "SMD")
  return(list("gene"=unique(df$hgnc_symbol),
              "logFC"=ma$TE.fixed,
              "logFC.l"=ma$lower.fixed,
              "logFC.u"=ma$upper.fixed,
              "p"=ma$pval.fixed))
  
}
ma_vals_d <- lapply(multi_probe_gene_d, ma_probes_genes)
mult_gene_df_d <- data.frame(apply(do.call(rbind, ma_vals_d) , c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), ~as.numeric(as.character(.)) ))
comb_ma_dat_d <- smok_int_d2 %>% filter(n==1) %>% 
  dplyr::rename(gene=hgnc_symbol, p=P.Value) %>%
  mutate(logFC.l=logFC-1.96*SD,
         logFC.u=logFC+1.96*SD) %>%
  dplyr::select(colnames(mult_gene_df_d)) %>%
  mutate(src="single") %>%
  bind_rows(mult_gene_df_d %>% mutate(src="mult")) %>%
  arrange(p)
comb_ma_dat_d$adj.p <- p.adjust(comb_ma_dat_d$p, method="fdr")
comb_ma_dat_d %>% filter(adj.p < 0.05) %>% nrow() # 7
comb_ma_dat_d %>% filter(adj.p < 0.05) %>% dplyr::select(-src)
intersect(comb_ma_dat_d  %>% filter(adj.p < 0.05) %>% pull(gene), unique(disc2$gene))

exp_long %>%
  filter(!is.na(hgnc_symbol), !is.na(sex), !is.na(smoking)) %>%
  filter(hgnc_symbol %in% (comb_ma_dat_d %>% filter(adj.p < 0.05) %>% pull(gene))) %>%
  unite(grp, c("smoking", "sex"), remove=FALSE) %>%
  ggplot(aes(y=expr, x=smoking, fill=sex, group=grp))+
  geom_boxplot()+
  facet_wrap(~hgnc_symbol, scales="free", nrow=2)+
  theme_bw()+
  ylab("")+
  xlab("")
ggsave("figures/example_genes2.png")

# --- all AE --- #
