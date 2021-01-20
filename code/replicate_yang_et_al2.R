
# TODO:
#  - check SL
#  - phenotype dedup
#  - go back to prcomp!
#  - update MA src to n, set it up to allow 

library('tidyverse')
library('Biobase')
library('ggrepel')
library('bigpca')
library('limma')
library('meta')
library('metafor')
source("code/00_utils.R")

# ---> download server side
# while IFS=',' read -r gse gpl submission_date fpath
# do
# wget "$fpath";
# done < "../list_ae_to_download_replication.csv"
# <--- read in

# ---> load server side
# library(affy)
# 
# setwd("data/cel_files_rep/")
# data1 <- ReadAffy(filenames=list.celfiles()) 
# save(data1, file="../data_rep.RData")
# 
# eset1 <- affy::rma(data1)
# save(eset1, file="../eset_rep.RData")

# --- 0. load data and de-duplicate --- #
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

# which(pcor2<1 & pcor2 > 0.99, arr.ind=T)
# # should we combine these?
# pDat3 %>% filter(gsm %in% c("GSM101106", "GSM101105"))
# pDat3 %>% filter(gsm %in% c("GSM114089", "GSM101115"))
# pDat3 %>% filter(gsm %in% c("GSM114090", "GSM101114"))
# pDat3 %>% filter(gsm %in% c("GSM219258", "GSM219252")) 
# pDat3 %>% filter(gsm %in% c("GSM219263", "GSM219260")) 
# pDat3 %>% filter(gsm %in% c("GSM469996", "GSM469994"))

pDat4 <- pDat3 %>% mutate(pack_years=ifelse(is.na(pack_years), 0, pack_years)) %>%
  mutate(across(c(age, pack_years), as.numeric)) %>%
  mutate(across(c(sex, smoking, race_ethnicity), as.factor))
rownames(pDat4) <- pDat4$gsm
expDat3 <- expDat2[,rownames(pDat4)]

#>> TODO: phenotype info de-duplication??
#pDat3 %>% 
#  group_by(age, sex, smoking,race_ethnicity, pack_years, submission_date) %>% 
#  count() %>%
#  filter(n>1)

# --- check SL  --- #


# --- 1. covariate viz --- #
pcs <- big.PCA(as.matrix(expDat3), pcs.to.keep=3); 
pcs2 <- data.frame(pcs$PCs)
pcs2$gsm <- rownames(pcs2) 
pcs2.2 <- pcs2 %>% left_join(pDat4, by=c("gsm"))

# double check no dups in PC space - there are none
# pcs2.2 %>% mutate(across(c(PC1, PC2, PC3), ~round(.,digits=4))) %>%
#   group_by(PC1, PC2, PC3) %>%
#   count() %>%
#   filter(n>1)

plotPC3(pcs2.2, sex)
to_exclude <- pcs2.2 %>% filter(PC1 > 0.1)

# exclude "outliers" and re-run PCA
pDat5 <- pDat4 %>% anti_join(to_exclude, by="gsm")

pcs_v2 <- big.PCA(as.matrix(expDat3[,pDat5$gsm]), pcs.to.keep=3)
pcs_v2_alt <- prcomp(t(as.matrix(expDat3[,pDat5$gsm])))

pcs22 <- data.frame(pcs_v2$PCs)
pcs22$gsm <- rownames(pcs22) 
pcs22.2 <- pcs22 %>% left_join(pDat5, by=c("gsm"))

# visualize PC results
plotPC3(pcs22.2, sex)
plotPC3(pcs22.2, smoking)
plotPC3(pcs22.2, race_ethnicity)
plotPC3(pcs22.2, submission_date)

# look at association btw PCs and covariates
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


# --- 1b. test for different covariate breakdown --- #
# TODO - make sure not hardcoded!
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


# --- 2. run models --- #

# original model
design_p <- model.matrix(~sex + smoking + sex:smoking + age + pack_years + race_ethnicity,data=pDat5)
head(expDat2[,1:5])
fit_p <- lmFit(expDat3[,pDat5$gsm], design_p)
fit_p <- eBayes(fit_p)
colnames(fit_p$coefficients)
add_gene_ids <- function(df){
  convert_genes <- load_gene_convert(df, "affy")
  df  %>% add_gene() %>% 
    left_join(convert_genes %>% 
                mutate(entrezgene_id=as.character(entrezgene_id)), by=c("gene"="entrezgene_id"))
}

tt_smok <- topTable(fit_p, coef="smokingS", n=nrow(expDat2)) %>%  add_gene_ids() 
tt_sex <- topTable(fit_p, coef="sexm", n=nrow(expDat2)) %>%  add_gene_ids() 
tt_int <- topTable(fit_p, coef="sexm:smokingS", n=nrow(expDat2)) %>%  add_gene_ids() 

plot_w_p <- function(df){
  df %>% dplyr::select(-gene) %>% 
    mutate(gene=hgnc_symbol) %>%
    dplyr::select(-adj.P.Val) %>%
    dplyr::rename(adj.P.Val=P.Value)
}

volcano_plot_de(tt_int %>% plot_w_p(), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_int2.png")

volcano_plot_de(tt_sex %>% plot_w_p(), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_sex2.png")

volcano_plot_de(tt_smok %>% plot_w_p(), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_smok2.png")

# model with date
design_date <- model.matrix(~sex + smoking + sex:smoking + age + pack_years + race_ethnicity+
                              submission_date,data=pDat5)
fit_d <- lmFit(expDat3[,pDat5$gsm], design_date)
fit_d <- eBayes(fit_d)
colnames(fit_d$coefficients)

tt_int_d <- topTable(fit_d, coef="sexm:smokingS", n=nrow(expDat2))  %>%  add_gene_ids()
tt_sex_d <- topTable(fit_d, coef="sexm", n=nrow(expDat2))  %>%  add_gene_ids()
tt_smok_d <- topTable(fit_d, coef="smokingS", n=nrow(expDat2))  %>%  add_gene_ids()

volcano_plot_de(tt_int_d %>% plot_w_p(), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_d_int2.png")

volcano_plot_de(tt_sex_d %>% plot_w_p(), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_d_sex2.png")

volcano_plot_de(tt_smok_d %>% plot_w_p(), pcut=0.01)+
  ylab("-log10 P value")
ggsave("figures/ae_rep_d_smok2.png")


# --- 3. compare results + Yang disc --- #
disc2 <- read_csv("ref/yang_int_disc.csv")
tt_int_d %>% inner_join(disc2, by=c("probes"="probe")) %>% 
  ggplot(aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.5)

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


# look at top genes CI
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



# plot example genes
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

# --- 4. model comparison --- #
probe_list <- exp_long %>% group_split(probes)
run_mod_comp <- function(df){
  mod1 <- lm(expr ~ age + sex+smoking+smoking*sex+race_ethnicity+pack_years, 
             data=df)
  mod2 <- lm(expr ~ age + sex+smoking+smoking*sex+race_ethnicity+pack_years+submission_date, 
             data=df)
  res <- anova(mod1, mod2)
  return(list("probe"=unique(df$probes), "Fstat"=res$`F`[[2]], "p"=res$`Pr(>F)`[[2]], 
              "SS"=res$`Sum of Sq`[[2]]))
}

names(probe_list) <- unique(exp_long$probes)
mod_comp2 <- lapply(probe_list[tt_int2$probes], run_mod_comp)
my_mods <- do.call(rbind,mod_comp2)
df <- data.frame(my_mods)
df2 <- df %>% mutate(probe=unlist(probe),
                     Fstat=unlist(Fstat), 
                     p=unlist(p),
                     SS=unlist(SS)) %>%
  as_tibble()
df2 %>% filter(SS < 0)
df2 %>% distinct() %>% arrange(p) %>% filter(p<0.05/1423)



# --- 5. probes to genes --- #

# run meta-analysis to convert probes to genes
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

# prepare data for MA
prep_data_ma <- function(df){
  df1 <- df %>% dplyr::select(-probe) %>% 
    left_join(probe_gene, by="probes") %>% 
    left_join(convert_genes %>% 
                mutate(entrezgene_id=as.character(entrezgene_id)),
              by=c("gene"="entrezgene_id")) 
  df2 <- df1 %>%
    mutate(SD=(logFC-CI.L)/1.96) %>% # SE for effect size - we have 95%  CI = 1.96*SD
    dplyr::select(probes, hgnc_symbol, logFC, SD, P.Value) %>%
    filter(!is.na(hgnc_symbol), hgnc_symbol!="") %>%
    distinct() %>% 
    group_by(hgnc_symbol) %>%
    mutate(n=n()) 
  return(df2)
}

# clean + run probe data thru meta-analysis, clean up output
# TODO update src to n
run_clean_ma <- function(df){
  df2 <- prep_data_ma(df)
  multi_probe_gene <- df2 %>% filter(n>1) %>% group_split(hgnc_symbol)
  ma_vals <- lapply(multi_probe_gene, ma_probes_genes)
  mult_gene_df <- data.frame(apply(do.call(rbind, ma_vals) , c(1,2), unlist)) %>%
    mutate(across(c(contains("logFC"), p), ~as.numeric(as.character(.)) ))
  comb_ma_dat <- df2 %>% filter(n==1) %>% 
    dplyr::rename(gene=hgnc_symbol, p=P.Value) %>%
    mutate(logFC.l=logFC-1.96*SD,
           logFC.u=logFC+1.96*SD) %>%
    dplyr::select(colnames(mult_gene_df)) %>%
    mutate(src="single") %>%
    bind_rows(mult_gene_df %>% mutate(src="mult")) %>% 
    arrange(p)
  comb_ma_dat$adj.p <- p.adjust(comb_ma_dat$p, method="fdr")
  return(comb_ma_dat)
}

# baseline model
comb_ma_dat <- run_clean_ma(smok_ci_int_p)
comb_ma_dat %>% filter(adj.p < 0.05) %>% nrow() # 179
comb_ma_dat %>% filter(adj.p < 0.05) %>% dplyr::select(-src) %>% head(10)

# date model
comb_ma_dat <- run_clean_ma(smok_ci_int_d)
comb_ma_dat_d %>% filter(adj.p < 0.05) %>% nrow() # 7
comb_ma_dat_d %>% filter(adj.p < 0.05) %>% dplyr::select(-src)
intersect(comb_ma_dat_d  %>% filter(adj.p < 0.05) %>% pull(gene), unique(disc2$gene))

# plot example genes!
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