# 00a_exploratory_analysis.R
#
# Exploratory analysis of expression data, examine PCs and look for duplicates
# This exploratory analysis focuses on the ae data (affy) -AND- does not
# identify expression duplicates.

library('tidyverse')
library('bigpca')
library('coop')
library('doppelgangR')

load("data/ae_eset_pheno.RData")

# --- plot PCs --- #

#pcs <- prcomp(t(eset_comb3), rank.=4, scale=TRUE)
#pcs2 <- data.frame(pcs$x)

# use bigpca b/c it is faster
start_time = Sys.time(); 
bp1 <- big.PCA(as.matrix(eset_comb3), pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time

# put together with the metadata
pcs2 <- data.frame(bp1$PCs)
pcs2$smoking <- comb_metadata3$smok
pcs2$sex <- comb_metadata3$expr_sex
pcs2$tissue <- comb_metadata3$tissue

pcs2$sample_acc <- rownames(pcs2) 

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



# ---- checking for expression duplicates ---- #

# check for duplicates in PC space by grouping
#  for AE processed the original way, this is NOT a problem
pcs2.3 <- pcs2 %>% 
  mutate(across(c(PC1, PC2, PC3), ~round(., digits=4))) %>%
  group_by(PC1, PC2, PC3) %>%
  summarize(sample_acc=paste(unique(sample_acc), collapse=";"),
    sex=paste(unique(sex), collapse=";"),
    smoking=paste(unique(smoking), collapse=";")) 


# look at pairwise correlations with coop (fast)
s_ae_only <- scale(ae_only)
row_vars <- apply(ae_only, 1, var)
quantile(row_vars, c(0.01, 0.05, 0.10)) # remove lowest percentile var
s2_ae_only <- s_ae_only[row_vars >0.005,]
pcc <- coop::pcor(s2_ae_only)
cors <- pcc[upper.tri(pcc, diag=FALSE)]

ggplot(tibble("cor"=cors), aes(x=cor))+
  geom_histogram(bins=1000)+
  theme_bw()

ggplot(tibble("cor"=cors), aes(x=cor))+
  geom_histogram(bins=1000)+
  theme_bw()+
  xlim(0.99, 1) # zoom in


# check for duplicates with doppelgangR
#  nothing pops out
ae_only_meta = data.frame(ae_only_meta)
rownames(ae_only_meta) <- ae_only_meta$geo_accession
ae_eset <- ExpressionSet(assayData=as.matrix(ae_only), 
                         phenoData=AnnotatedDataFrame(ae_only_meta))

dR <- doppelgangR(ae_eset, phenoFinder.args=NULL)
plot(dR)

# change the sens/spec (uses cached version so faster)
dR1 <- doppelgangR(ae_eset, phenoFinder.args=NULL, 
                   outlierFinder.expr.args = list(bonf.prob = 1.0, tail = "upper"))
plot(dR1)
