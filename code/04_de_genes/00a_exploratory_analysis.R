# 00a_exploratory_analysis.R
#
# Exploratory analysis of expression data, examine Pcs, etc

library('tidyverse')
load("data/ae_eset_pheno.RData")

# --- plot PCs --- #

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
