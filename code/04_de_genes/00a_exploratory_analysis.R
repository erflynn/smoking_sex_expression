# 00a_exploratory_analysis.R
#
# Exploratory analysis of expression data, examine Pcs, etc

library('tidyverse')
load("data/ae_eset_pheno.RData")

# --- plot PCs --- #

#pcs <- prcomp(t(eset_comb3), rank.=4, scale=TRUE)
#pcs2 <- data.frame(pcs$x)
library(bigpca)
start_time = Sys.time(); 
bp1 <- big.PCA(as.matrix(eset_comb3), pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time
pcs2 <- data.frame(bp1$PCs)
pcs2$smoking <- comb_metadata3$smok
pcs2$sex <- comb_metadata3$expr_sex
pcs2$tissue <- comb_metadata3$tissue

pcs2$sample_acc <- rownames(pcs2) 
pcs2.3 <- pcs2 %>% 
  #arrange(PC1, PC2, PC3) %>%
  mutate(across(c(PC1, PC2, PC3), ~round(., digits=4))) %>%
  group_by(PC1, PC2, PC3) %>%
  summarize(#study=paste(unique(study), collapse=";"),
            sample_acc=paste(unique(sample_acc), collapse=";"),
            sex=paste(unique(sex), collapse=";"),
            smoking=paste(unique(smoking), collapse=";")) 
# ... ok mb there isn't a problem with these data? ... #
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
